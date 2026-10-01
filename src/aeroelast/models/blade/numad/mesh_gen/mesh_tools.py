import numpy as np

from aeroelast.models.blade.numad.mesh_gen.element_utils import (
    cross_prod,
)
from aeroelast.models.blade.numad.mesh_gen.spatial_grid_list2d import spatial_grid_list2d
from aeroelast.models.blade.numad.mesh_gen.spatial_grid_list3d import spatial_grid_list3d


def rotate_vector(vec, axis, angle):
    if angle < 0.0000000001:
        return vec.copy()
    else:
        axAr = np.array(axis)
        mag = np.linalg.norm(axis)
        unitAxis = (1.0 / mag) * axAr
        alp1 = np.zeros((3, 3), dtype=float)
        alp1[0] = unitAxis
        i1 = 0
        if abs(unitAxis[1]) < abs(unitAxis[0]):
            i1 = 1
        if abs(unitAxis[2]) < abs(unitAxis[i1]):
            i1 = 2
        alp1[1, i1] = np.sqrt(1.0 - alp1[0, i1] * alp1[0, i1])
        for i2 in range(0, 3):
            if i2 != i1:
                alp1[1, i2] = -alp1[0, i1] * alp1[0, i2] / alp1[1, i1]
        alp1[2] = cross_prod(alp1[0], alp1[1])
        theta = angle * np.pi / 180.0
        cs = np.cos(theta)
        sn = np.sin(theta)
        alp2 = np.array([[1.0, 0.0, 0.0], [0.0, cs, -sn], [0.0, sn, cs]])
        rV = np.matmul(alp1, vec)
        rV = np.matmul(alp2, rV)
        rV = np.matmul(rV, alp1)
        return rV


def translate_mesh(meshData, tVec):
    tAr = np.array(tVec)
    nLen = len(meshData["nodes"])
    newNds = np.zeros((nLen, 3), dtype=float)
    for i, nd in enumerate(meshData["nodes"]):
        newNds[i] = nd + tAr
    meshData["nodes"] = newNds
    return meshData


def rotate_mesh(meshData, pt, axis, angle):
    ptAr = np.array(pt)
    nLen = len(meshData["nodes"])
    newNds = np.zeros((nLen, 3), dtype=float)
    for i, nd in enumerate(meshData["nodes"]):
        tCrd = nd - ptAr
        rCrd = rotate_vector(tCrd, axis, angle)
        newNds[i] = ptAr + rCrd
    meshData["nodes"] = newNds
    return meshData


def get_direction_cosines(xDir, xyDir):
    mag = np.linalg.norm(xDir)
    a1 = (1.0 / mag) * xDir
    zDir = cross_prod(xDir, xyDir)
    mag = np.linalg.norm(zDir)
    a3 = (1.0 / mag) * zDir
    a2 = cross_prod(a3, a1)
    dirCos = np.array([a1, a2, a3])
    return dirCos


def get_average_node_spacing(nodes, elements) -> float:
    totDist = 0.0
    ct = 0
    for el in elements:
        for ndi in el:
            if ndi > -1:
                nd1 = nodes[ndi]
                for ndi2 in el:
                    if ndi2 > -1 and ndi2 != ndi:
                        nd2 = nodes[ndi2]
                        vec = nd1 - nd2
                        dist = np.linalg.norm(vec)
                        totDist = totDist + dist
                        ct = ct + 1
    return totDist / ct


def get_mesh_spatial_list(
    nodes, xSpacing: float = 0, ySpacing: float = 0, zSpacing: float = 0
):
    totNds = len(nodes)
    spaceDim = len(nodes[0])

    maxX = np.amax(nodes[:, 0])
    minX = np.amin(nodes[:, 0])
    maxY = np.amax(nodes[:, 1])
    minY = np.amin(nodes[:, 1])
    nto1_2 = np.power(totNds, 0.5)
    nto1_3 = np.power(totNds, 0.333333333333)
    if spaceDim == 3:
        maxZ = np.amax(nodes[:, 2])
        minZ = np.amin(nodes[:, 2])
        dimVec = np.array([(maxX - minX), (maxY - minY), (maxZ - minZ)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01 * meshDim
        minX = minX - 0.01 * meshDim
        maxY = maxY + 0.01 * meshDim
        minY = minY - 0.01 * meshDim
        maxZ = maxZ + 0.01 * meshDim
        minZ = minZ - 0.01 * meshDim
        if xSpacing == 0:
            xS = 0.5 * (maxX - minX) / nto1_3
        else:
            xS = xSpacing
        if ySpacing == 0:
            yS = 0.5 * (maxY - minY) / nto1_3
        else:
            yS = ySpacing
        if zSpacing == 0:
            zS = 0.5 * (maxZ - minZ) / nto1_3
        else:
            zS = zSpacing
        meshGL = spatial_grid_list3d(minX, maxX, minY, maxY, minZ, maxZ, xS, yS, zS)
        # tol = 1.0e-6*meshDim/nto1_3
    else:
        dimVec = np.array([(maxX - minX), (maxY - minY)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01 * meshDim
        minX = minX - 0.01 * meshDim
        maxY = maxY + 0.01 * meshDim
        minY = minY - 0.01 * meshDim
        if xSpacing == 0:
            xS = 0.5 * (maxX - minX) / nto1_2
        else:
            xS = xSpacing
        if ySpacing == 0:
            yS = 0.5 * (maxY - minY) / nto1_2
        else:
            yS = ySpacing
        meshGL = spatial_grid_list2d(minX, maxX, minY, maxY, xS, yS)
        # tol = 1.0e-6*meshDim/nto1_2
    return meshGL


## - Convert list of mesh objects into a single merged mesh, returning sets representing the elements/nodes from the original meshes
def mergeDuplicateNodes(meshData, tolerance=None):
    allNds = meshData["nodes"]
    allEls = meshData["elements"]
    totNds = len(allNds)
    totEls = len(allEls)
    elDim = len(allEls[0])

    avgSp = get_average_node_spacing(meshData["nodes"], meshData["elements"])
    sp = 2 * avgSp
    nodeGL = get_mesh_spatial_list(allNds, xSpacing=sp, ySpacing=sp, zSpacing=sp)
    if tolerance is None:
        tol = 1.0e-4 * avgSp
    else:
        tol = tolerance

    # glDim = nodeGL.getDim()
    # mag = np.linalg.norm(glDim)
    # nto1_3 = np.power(len(allNds),0.3333333333)
    # tol = 1.0e-6*mag/nto1_3

    i = 0
    for nd in allNds:
        nodeGL.addEntry(i, nd)
        i = i + 1

    ndElim = -np.ones(totNds, dtype=int)
    ndNewInd = -np.ones(totNds, dtype=int)
    for n1i in range(0, totNds):
        if ndElim[n1i] == -1:
            nearNds = nodeGL.findInRadius(allNds[n1i], tol)
            for n2i in nearNds:
                if n2i > n1i and ndElim[n2i] == -1:
                    proj = allNds[n2i] - allNds[n1i]
                    dist = np.linalg.norm(proj)
                    if dist < tol:
                        ndElim[n2i] = n1i
    ndi = 0
    nodesFinal = []
    for n1i in range(0, totNds):
        if ndElim[n1i] == -1:
            nodesFinal.append(allNds[n1i])
            ndNewInd[n1i] = ndi
            ndi = ndi + 1
    nodesFinal = np.array(nodesFinal)
    for eli in range(0, totEls):
        for j in range(0, elDim):
            nd = allEls[eli, j]
            if nd != -1:
                if ndElim[nd] == -1:
                    allEls[eli, j] = ndNewInd[nd]
                else:
                    allEls[eli, j] = ndNewInd[ndElim[nd]]

    meshData["nodes"] = nodesFinal
    meshData["elements"] = allEls

    return meshData


def add_node_set(meshData, newSet):
    try:
        meshData["sets"]["node"].append(newSet)
    except Exception:
        nSets = []
        nSets.append(newSet)
        try:
            meshData["sets"]["node"] = nSets
        except Exception:
            sets = {}
            sets["node"] = nSets
            meshData["sets"] = sets
    return meshData


def add_element_set(meshData, newSet):
    try:
        meshData["sets"]["element"].append(newSet)
    except Exception:
        elSets = []
        elSets.append(newSet)
        try:
            meshData["sets"]["element"] = elSets
        except Exception:
            sets = {}
            sets["element"] = elSets
            meshData["sets"] = sets
    return meshData
