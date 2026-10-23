from copy import deepcopy

import numpy as np
from numpy import ndarray

from aeroelast.models.blade.numad.objects.bom import BillOfMaterials
from aeroelast.models.blade.numad.objects.keypoints import KeyPoints
from aeroelast.models.blade.numad.utils.interpolation import interpolator_wrap


class StackDatabase:
    def __init__(self):
        self.stacks: ndarray | None = None
        self.swstacks: ndarray | None = None

    def __eq__(self, other):
        if self.stacks is None or self.swstacks is None:
            raise RuntimeError("StackDatabase.__eq__ requires a generated database")
        assert self.stacks.shape == other.stacks.shape

        assert self.swstacks.shape == other.swstacks.shape

        self_flat_stacks = self.stacks.flatten()
        other_flat_stacks = other.stacks.flatten()
        for i in range(self_flat_stacks.shape[0]):
            self_stack = self_flat_stacks[i]
            other_stack = other_flat_stacks[i]
            if self_stack != other_stack:
                return False
        self_flat_stacks = self.swstacks.flatten()
        other_flat_stacks = other.swstacks.flatten()
        for i in range(self_flat_stacks.shape[0]):
            self_stack = self_flat_stacks[i]
            other_stack = other_flat_stacks[i]
            if self_stack != other_stack:
                return False
        return True

    def generate(self, keypoints: KeyPoints, bom: BillOfMaterials, definition=None):
        # build the material stack for each area
        if (
            keypoints.key_areas is None
            or keypoints.web_points is None
            or keypoints.web_indices is None
        ):
            raise RuntimeError("StackDatabase.generate() requires generated keypoints")
        n_segments = keypoints.key_areas.shape[0]
        n_stations = keypoints.key_areas.shape[1]
        n_webs = len(keypoints.web_points)
        segment_labels = [
            "HP_TE_FLAT",
            "HP_TE_REINF",
            "HP_TE_PANEL",
            "HP_SPAR",
            "HP_LE_PANEL",
            "HP_LE",
            "LP_LE",
            "LP_LE_PANEL",
            "LP_SPAR",
            "LP_TE_PANEL",
            "LP_TE_REINF",
            "LP_TE_FLAT",
        ]

        # initialize stack array
        self.stacks = np.empty(shape=(n_segments, n_stations), dtype=object)
        for k_seg in range(n_segments):
            for k_stat in range(n_stations):
                self.stacks[k_seg, k_stat] = Stack()
                self.stacks[k_seg, k_stat].name = "{:02d}_{:02d}_{}".format(
                    k_seg, k_stat, segment_labels[k_seg]
                )
                self.stacks[k_seg, k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    k_seg,
                    k_seg + 1,
                ]

        for k in range(len(bom["hp"])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = bom["hp"][k].name
            cur_ply.materialid = bom["hp"][k].materialid
            cur_ply.thickness = bom["hp"][k].thickness
            cur_ply.angle = 0  # TODO, set to 0 for now, bom['lp'](k, );
            cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = bom.indices["hp"][k]
            for k_seg in range(ind[2], ind[3]):
                for k_stat in range(ind[0], ind[1]):
                    # deepcopy is important to keep make ply object unique in each stack
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))

        for k in range(len(bom["lp"])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = bom["lp"][k].name
            cur_ply.materialid = bom["lp"][k].materialid
            cur_ply.thickness = bom["lp"][k].thickness
            cur_ply.angle = 0  # TODO, set to 0 for now, bom['lp'](k, );
            cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = bom.indices["lp"][k]
            for k_seg in range(ind[2], ind[3]):
                for k_stat in range(ind[0], ind[1]):
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))

        self.swstacks = np.empty(shape=(n_webs, n_stations), dtype=object)
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                self.swstacks[k_web, k_stat] = Stack()
                # name the stacks <webnumber+1>_<stationnumber+1>_SW
                self.swstacks[k_web, k_stat].name = "{:02d}_{:02d}_SW".format(k_web, k_stat)
                # currently, the shearweb indices do not change down the span
                ind = keypoints.web_indices[k_web]
                self.swstacks[k_web][k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    ind[0],
                    ind[1],
                ]
        for k_web in range(n_webs):
            for k in range(len(bom["sw"][k_web])):
                # for each row in the BOM, get the ply definition ...
                cur_ply = Ply()
                cur_ply.component = bom["sw"][k_web][k].name
                cur_ply.materialid = bom["sw"][k_web][k].materialid
                cur_ply.thickness = bom["sw"][k_web][k].thickness
                cur_ply.angle = 0  # TODO, set to 0 for now, bom['lp'](k, );
                cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

                # ... and add the ply to every area that is part of the region
                ind = bom.indices["sw"][k_web][k]
                for k_stat in range(ind[0], ind[1]):
                    self.swstacks[k_web, k_stat].addply(deepcopy(cur_ply))

        if definition is not None and any(
            c.continuous_layers for c in definition.components.values()
        ):
            self._map_continuous_layers(definition, keypoints)

        self._order_skin_plies_for_outward_normals()

    def _order_skin_plies_for_outward_normals(self):
        """Reverse the skin stacks so the outermost layer is the last ply.

        Layers (windIO ``layers``, Excel BOM rows) are listed from the outer mould line inward and
        are added in that order, but the shell laminate stacks ply 0 on the ``-normal`` side and the
        skin elements' normals point outward. Emitting the stack inner-to-outer puts the gelcoat on
        the outer face and the innermost layer on the inner face (stress recovery TOP = outermost).
        Web stacks are symmetric about their core and are left alone.
        """
        for stack in self.stacks.flat:
            stack.plygroups.reverse()

    def _map_continuous_layers(self, definition, keypoints):
        """Map WindIO laminates without converting thickness to integer stock plies.

        Keep the existing NuMAD region topology. Integrate thickness across
        each span bay (two-point Gauss is exact for linear/cubic interpolation),
        rather than holding the inboard value over a thinning bay. Evaluate
        skin ND bounds at the bay midpoint. At a boundary, use the arc-overlap fraction of
        the region: this is a conservative, region-averaged laminate, NOT a
        new geometric edge. Thickness is unchanged in fully covered regions.
        The YAML skin bounds are symmetric about 0.5; reject asymmetric ones
        rather than silently guessing the WindIO/NuMAD side correspondence.
        Excel components never enter this path.
        """
        span = definition.ispan
        eta = (span-span[0])/(span[-1]-span[0])
        mid = .5*(eta[:-1]+eta[1:])
        offset = .5*np.diff(eta)/np.sqrt(3.)
        for stack in self.stacks.flat:
            stack.plygroups = []
        for stack in self.swstacks.flat:
            stack.plygroups = []
        for component in definition.components.values():
            if not component.continuous_layers:
                raise ValueError("Cannot mix continuous WindIO and integer-ply components")
            counts = .5*(np.maximum(component.get_num_layers(mid-offset), 0.)
                         + np.maximum(component.get_num_layers(mid+offset), 0.))
            material = definition.materials[component.materialid]
            if component.group:
                regions = [(self.swstacks[component.group-1], None)]
            else:
                hp, lp = component.find_region_extents()
                regions = [(self.stacks[j], j) for region in (hp, lp) if region
                           for j in range(*region)]
            start = np.zeros_like(mid)
            end = np.ones_like(mid)
            if component.nd_arc_start is not None and component.nd_arc_end is not None:
                start = interpolator_wrap(component.nd_arc_start["grid"],
                                          component.nd_arc_start["values"], mid, "pchip")
                end = interpolator_wrap(component.nd_arc_end["grid"],
                                        component.nd_arc_end["values"], mid, "pchip")
                if not np.allclose(start+end, 1., rtol=0., atol=1e-12):
                    raise ValueError("Asymmetric WindIO skin ND bounds are not supported")
                if np.any(start < 0.) or np.any(end > 1.) or np.any(start > end):
                    raise ValueError("Invalid WindIO skin ND bounds")
            for stacks, segment in regions:
                for station, stack in enumerate(stacks):
                    fraction = 1.
                    if segment is not None:
                        arcs = keypoints.key_arcs[:, station]
                        perimeter = arcs[-1]-arcs[0]
                        a, b = (arcs[segment:segment+2]-arcs[0])/perimeter
                        if b <= a:
                            # A sharp TE has coincident key arcs before TE expansion.
                            fraction = float(start[station] <= a <= end[station])
                        else:
                            fraction = max(0., min(b, end[station])-max(a, start[station]))/(b-a)
                    count = counts[station]*fraction
                    if count <= 0.:
                        continue
                    ply = Ply()
                    ply.component = component.name
                    ply.materialid = component.materialid
                    ply.thickness = material.layerthickness
                    ply.angle = 0  # Preserve NuMAD's existing orientation convention.
                    ply.nPlies = float(count)
                    stack.addply(ply)


class Stack:
    """A class definition for a stack of composite layers.

    Parameters
    ----------

    Attributes
    ----------
    name : string
        Name of the stack or composite material used by NuMAD, e.g. '000000_HP_LE_PANEL'
    indices : list
        Indices of stack, ``[in board station, out board station,
        1st kepoint, 2nd keypoint]``, e.g. ``[ibSta,obSta,keypt1,keypt2]``
    plygroups : list
        List of ``ply`` objects

    """

    def __init__(self):
        self.name: str = ""
        self.indices = []
        self.plygroups: list = []

    def __eq__(self, other):
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int, float, str, list, dict)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True

    def addply(self, ply):
        """This method adds a Ply object to stack

        Parameters
        ----------
        ply : Ply object
            Ply object to be added
        Returns
        -------
        None

        """
        # modifies
        plygroups = self.plygroups
        if (
            plygroups
            and ply.component == plygroups[-1].component
            and ply.angle == plygroups[-1].angle
        ):
            plygroups[-1].nPlies += ply.nPlies
        else:
            plygroups.append(ply)

        return self


class Ply:
    """A simple class to organize the attributes of a ply

    Attributes
    ----------

    component : str
        parent component
    materialid : str
        Material id of ply
    thickness : float
        thickness of single ply (mm)
    angle : float
        ply angle
    nPlies : float
        number of plies (fractional for continuous WindIO laminate thickness)
    """

    def __init__(self):
        self.component: str | None = None
        self.materialid: str | None = None
        self.thickness: float | None = None
        self.angle: float | None = None
        self.nPlies: float | None = None

    def __eq__(self, other):
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int, float, str)):
                if self_attr != other_attr:
                    return False
        return True
