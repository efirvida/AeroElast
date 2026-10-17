"""Uniform shell boxes and independent beam references for T6 (SI units).

Geometry is explicit: generic ShellProperty has no reference-surface offset.
Moving the four wall planes therefore tests a MODEL effect, not BladeMesh.
All assembly/centrifugal recovery/spin softening is delegated to S-4 unchanged.
"""

import numpy as np
from numpy.polynomial import Legendre
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.laminate import create_laminate_from_angles
from aeroelast.core.material import IsotropicMaterial, Material
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import CompositeShellProperty, ShellProperty
from aeroelast.elements import ElementFamily
from tools.run_s4_rotating_modal import assemble_rotating
from tests.validation.parity.test_composite_beam_parity import _hand_clt_abd

LENGTH, WIDTH, HEIGHT, THICKNESS = 30.0, 0.6, 1.0, 0.05
E, NU, RHO = 70e9, 0.3, 2700.0
ANGLES = [0, 45, -45, 90, 90, -45, 45, 0]
GLASS = (45e9, 12e9, 5e9, 0.28)


def section(width=WIDTH, height=HEIGHT, thickness=THICKNESS):
    """Integrate four wall strips, including their local t^3/12 bending terms.

    Corner overlap is intentionally absent: four shell strips, NOT the exact
    area of a solid rectangle minus a smaller rectangle. Mass is rho*p*t.
    """
    w, h, t = width, height, thickness
    area = 2 * (w + h) * t
    ix = t * (w * h**2 / 2 + h**3 / 6) + w * t**3 / 6
    iy = t * (h * w**2 / 2 + w**3 / 6) + h * t**3 / 6
    j = 2 * w**2 * h**2 * t / (w + h)
    return area, ix, iy, j


def build_box(nspan, nwall, width=WIDTH, height=HEIGHT, thickness=THICKNESS,
              composite=False, triangles=False):
    """Shared corners, outward normals, span +z, axis +y, root at r0=0."""
    mesh = MeshModel()
    corners = np.array([[-width / 2, -height / 2], [width / 2, -height / 2],
                        [width / 2, height / 2], [-width / 2, height / 2]])
    ring = np.array([a + (b - a) * i / nwall
                     for a, b in zip(corners, np.roll(corners, -1, axis=0), strict=True)
                     for i in range(nwall)])
    nodes = []
    for z in np.linspace(0, LENGTH, nspan + 1):
        row = []
        for x, y in ring:
            node = Node([float(x), float(y), float(z)], geometric_node=False)
            mesh.add_node(node)
            row.append(node)
        nodes.append(row)
    for lo, hi in zip(nodes[:-1], nodes[1:], strict=True):
        for i in range(len(ring)):
            j = (i + 1) % len(ring)
            if triangles:
                mesh.add_element(MeshElement(nodes=[lo[i], lo[j], hi[j]],
                                             element_type=ElementType.triangle))
                mesh.add_element(MeshElement(nodes=[lo[i], hi[j], hi[i]],
                                             element_type=ElementType.triangle))
            else:
                mesh.add_element(MeshElement(nodes=[lo[i], lo[j], hi[j], hi[i]],
                                             element_type=ElementType.quad))
    mesh.add_node_set(NodeSet("RootNodes", set(nodes[0])))
    mesh.add_element_set(ElementSet("walls", set(mesh.elements)))
    if composite:
        e1, e2, g, nu = GLASS
        material = Material(name="glass", E=(e1, e2, e2), G=(g, g, g),
                            nu=(nu, nu, 0), rho=RHO)
        prop = CompositeShellProperty(create_laminate_from_angles(
            material, thickness / len(ANGLES), ANGLES))
    else:
        prop = ShellProperty(IsotropicMaterial(name="aluminium", E=E, nu=NU, rho=RHO),
                             thickness)
    asm = MeshAssembler(mesh, {"elements": {"element_family": ElementFamily.SHELL,
                        "span_direction": (0., 0., 1.), "properties": {"walls": prop}}})
    # MeshModel owns all IDs and ordering; do not reset process-global counters.
    indices = np.array([[mesh.node_id_to_index[n.id] for n in row] for row in nodes])
    return asm, mesh, indices, ring


def modal_families(case, omega=0):
    """Classify section-average rigid motions, not one arbitrary corner node.

    Torsion uses the least-squares rigid twist and its transverse displacement
    norm. Reject local/distortional modes with <80% rigid-section content.
    """
    asm, mesh, indices, ring = case
    # PETSc's built-in LU default graph ordering hangs in MatGetOrdering on
    # this regular box (native trace, SLURM11617667). RCM changes only the
    # factorization permutation; the unchanged S-4 matrices and eigenproblem
    # give the identical spectrum (dense/SLEPc diagnostic cross-check).
    from petsc4py import PETSc

    options = PETSc.Options()
    key = "st_pc_factor_mat_ordering_type"
    previous = options.getString(key) if options.hasName(key) else None
    options[key] = "rcm"
    try:
        frequencies, modes, free = assemble_rotating(asm, mesh, omega, 12)
    finally:
        if previous is None:
            del options[key]
        else:
            options[key] = previous
    families = {"flap": [], "edge": [], "torsion": [], "axial": []}
    for frequency, mode in sorted(zip(frequencies, modes, strict=True)):
        full = np.zeros(asm.dofs_count)
        full[free] = mode
        u = full.reshape(-1, 6)[indices, :3]
        mean = u.mean(axis=1)
        twist = (u[:, :, 1] @ ring[:, 0] - u[:, :, 0] @ ring[:, 1]) / (ring**2).sum()
        scores = [np.sum(mean[:, 1]**2), np.sum(mean[:, 0]**2),
                  np.sum(twist**2) * np.mean(np.sum(ring**2, axis=1)),
                  np.sum(mean[:, 2]**2)]
        if sum(scores) >= 0.8 * np.sum(u**2) / len(ring):
            families[list(families)[int(np.argmax(scores))]].append(float(frequency))
    return families


def assembled_mass(case):
    """Rigid-translation quadratic form, independently for all three axes."""
    asm = case[0]
    rows, cols, values = asm._rust.assemble_m()
    mass = coo_matrix((values, (rows, cols)), shape=(asm.dofs_count, asm.dofs_count)).tocsr()
    result = []
    for axis in range(3):
        v = np.zeros(asm.dofs_count)
        v[axis::6] = 1
        result.append(float(v @ mass @ v))
    return np.array(result)


def shear_areas():
    """Thin-wall shear energy: q'=-t*y/I, closed-cell integral(q/t ds)=0.

    Unit V, A_s=1/integral(q^2/t ds). This is a section-derived shear
    coefficient A_s/A, NOT the solid-rectangle 5/6. Local t^3 terms omitted
    consistently in this thin-wall shear-flow calculation (O(t^2/h^2)).
    """
    corners = np.array([[-WIDTH/2, -HEIGHT/2], [WIDTH/2, -HEIGHT/2],
                        [WIDTH/2, HEIGHT/2], [-WIDTH/2, HEIGHT/2]])
    points, weights = np.polynomial.legendre.leggauss(8)
    result = []
    for direction, inertia in [(1, section(thickness=1)[1]-WIDTH/6),
                                (0, section(thickness=1)[2]-HEIGHT/6)]:
        qs, ws, qend = [], [], 0.0
        for a, b in zip(corners, np.roll(corners, -1, axis=0), strict=True):
            length = np.linalg.norm(b-a)
            s = (points+1)*length/2
            # t/I = 1 / (I/t).
            q = qend-(a[direction]*s+(b[direction]-a[direction])*s*s/(2*length))/inertia
            qs.extend(q)
            ws.extend(weights*length/2)
            qend -= (a[direction]+b[direction])*length/(2*inertia)
        qs, ws = np.array(qs), np.array(ws)
        qs -= np.sum(qs*ws)/sum(ws)
        result.append(THICKNESS/np.sum(qs*qs*ws))
    return np.array(result)


def timoshenko_reference(degree=16):
    """Independent 1D Timoshenko energy Ritz, including rotary inertia.

    K=integral(EI theta'^2 + G A_s (w'-theta)^2); M=integral(m w^2+
    rho I theta^2). Both fields vanish at the clamp; natural free end.
    Integrated Legendre basis avoids monomial conditioning. No shell matrices.
    """
    points, weights = np.polynomial.legendre.leggauss(80)
    x, weights = (points+1)/2, weights*LENGTH/2
    q, dq = [], []
    for i in range(degree):
        p = Legendre.basis(i, domain=[0, 1]).integ()
        p -= p(0)
        q.append(p(x))
        dq.append(p.deriv()(x)/LENGTH)
    q, dq = np.array(q), np.array(dq)
    mm = (q*weights) @ q.T
    dd = (dq*weights) @ dq.T
    cross = (dq*weights) @ q.T
    area, ix, iy, _ = section()
    result = {}
    for family, inertia, shear_area in zip(["flap", "edge"], [ix, iy], shear_areas(), strict=True):
        ga = E/(2*(1+NU))*shear_area
        stiffness = np.block([[ga*dd, -ga*cross],
                              [-ga*cross.T, E*inertia*dd+ga*mm]])
        mass = np.block([[RHO*area*mm, np.zeros_like(mm)],
                         [np.zeros_like(mm), RHO*inertia*mm]])
        result[family] = np.sqrt(eigh(stiffness, mass, subset_by_index=[0, 1])[0])/(2*np.pi)
    return result


def parked_reference():
    """Euler-Bernoulli bending; Saint-Venant/Bredt torsional wave cantilever.

    Blevins, Formulas for Natural Frequency and Mode Shape: beta L = 1.875104,
    4.694091; torsion k L=pi/2. Rotary inertia/shear are O((h/L)^2):
    their effect is bounded separately by the independently assembled 1D model.
    """
    area, ix, iy, j = section()
    beta = np.array([1.875104068711961, 4.694091132974175])
    return {"flap": beta**2 * np.sqrt(E * ix / (RHO * area)) / (2*np.pi*LENGTH**2),
            "edge": beta**2 * np.sqrt(E * iy / (RHO * area)) / (2*np.pi*LENGTH**2),
            "torsion": np.array([np.sqrt(E / (2*(1+NU)) * j / (RHO*(ix+iy)))
                                  / (4*LENGTH)])}


def rotating_beam_ritz(speed_ratio=0.4, degree=12, order=0):
    """Independent weak rotating EB ODE, r0=0, dimensionless x in [0,1].

    y'''' - (T y')' = lambda*y, T=mu^2*(1-x^2)/2;
    clamped y(0)=y'(0)=0, natural free tip. Integrated Legendre basis enforces
    only essential BCs. mu=Omega*sqrt(m L^4/EI). K1 below is the finite-speed
    Southwell coefficient, not assumed speed independent. Low-speed limit
    ~1.193: Wright, Smith, Thresher & Wang (1982), J. Appl. Mech. 49, 197-202,
    'Vibration modes of centrifugally stiffened beams'. Lead-lag subtracts 1.
    """
    points, weights = np.polynomial.legendre.leggauss(80)
    x, w = (points + 1)/2, weights/2
    q, dq, ddq = [], [], []
    for i in range(degree):
        curvature = Legendre.basis(i, domain=[0, 1])
        slope = curvature.integ()
        slope = slope - slope(0)
        shape = slope.integ()
        shape = shape - shape(0)
        q.append(shape(x))
        dq.append(slope(x))
        ddq.append(curvature(x))
    q, dq, ddq = np.array(q), np.array(dq), np.array(ddq)
    mass = (q*w) @ q.T
    stiffness = (ddq*w) @ ddq.T
    geometric = (dq*(w*(1-x*x)/2)) @ dq.T
    lam0, vectors = eigh(stiffness, mass, subset_by_index=[order, order])
    lam0 = lam0[0]
    mu2 = speed_ratio**2 * lam0
    lam = eigh(stiffness + mu2*geometric, mass, subset_by_index=[order, order])[0][0]
    vector = vectors[:, 0]
    low_speed = float(vector @ geometric @ vector)
    return float((lam-lam0)/mu2), low_speed


def static_rigidities(case):
    """EA from interior extension slope; EI from cubic section-mean bending.

    Saint-Venant interior (0.2L..0.8L), avoiding tip wall distortion and root
    Poisson restraint; shear adds a linear term so it cannot bias curvature.
    Unit total nodal forces on the tip. No reimplementation of material CLT.
    """
    asm, mesh, indices, _ = case
    rows, cols, vals = asm._rust.assemble_k()
    k = coo_matrix((vals, (rows, cols)), shape=(asm.dofs_count, asm.dofs_count)).tocsr()
    fixed = (indices[0, :, None]*6 + np.arange(6)).ravel()
    free = np.setdiff1d(np.arange(asm.dofs_count), fixed)
    force = np.zeros((asm.dofs_count, 2))
    force[indices[-1]*6+2, 0] = 1/indices.shape[1]
    force[indices[-1]*6+1, 1] = 1/indices.shape[1]
    u = np.zeros_like(force)
    u[free] = spsolve(k[free][:, free].tocsc(), force[free])
    uz = u[indices*6+2, 0].mean(axis=1)
    uy = u[indices*6+1, 1].mean(axis=1)
    z = np.array([mesh.nodes[int(i)].coords[2] for i in indices[:, 0]])
    interior = (z >= .2*LENGTH) & (z <= .8*LENGTH)
    ea = 1 / np.polyfit(z[interior], uz[interior], 1)[0]
    curve = np.polyfit(z[interior], uy[interior], 3)
    curvature = np.polyval(np.polyder(curve, 2), LENGTH/2)
    return np.array([ea, (LENGTH/2)/abs(curvature)])


def composite_reference(thickness=0.01):
    """Jones CLT: A_eff=1/a11=A11-A12^2/A22 (balanced), NOT A11.

    Independent hand CLT already used by parity tests; no production ABD data
    enters the judge. Local wall D is negligible for t/h=.01.
    """
    a, _, _ = _hand_clt_abd(*GLASS, ANGLES, thickness/len(ANGLES))
    effective = 1/np.linalg.inv(a)[0, 0]
    return effective * np.array([2*(WIDTH+HEIGHT), WIDTH*HEIGHT**2/2+HEIGHT**3/6])
