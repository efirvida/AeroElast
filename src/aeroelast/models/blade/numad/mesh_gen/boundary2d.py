import numpy as np

import aeroelast.models.blade.numad.mesh_gen.mesh_tools as mt
from aeroelast.models.blade.numad.mesh_gen.segment2d import Segment2D


class Boundary2D:
    def __init__(self, segList=None):
        if segList is None:
            segList = []
        self.segList = []
        self.segList.extend(segList)

    def addSegment(self, segType, keyPts, numEls):
        self.segList.append(Segment2D(segType, keyPts, numEls))

    def getBoundaryMesh(self):
        allNds = []
        allEds = []
        totNds = 0
        for seg in self.segList:
            segMesh = seg.getNodesEdges()
            allNds.extend(segMesh["nodes"])
            allEds.extend(segMesh["edges"] + totNds)
            totNds = len(allNds)
        allNds = np.array(allNds)
        allEds = np.array(allEds)

        meshData = {}
        meshData["nodes"] = allNds
        meshData["elements"] = allEds

        output = mt.mergeDuplicateNodes(meshData)

        return output
