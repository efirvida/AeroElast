import numpy as np

import aeroelast.models.blade.numad.mesh_gen.mesh_tools as mt
from aeroelast.models.blade.numad.mesh_gen.shell_region import ShellRegion


class Surface:
    def __init__(self, regionList=None, regionNames=None, meshList=None, meshNames=None):
        if meshNames is None:
            meshNames = []
        if meshList is None:
            meshList = []
        if regionNames is None:
            regionNames = []
        if regionList is None:
            regionList = []
        self.shellRegions = []
        self.shellRegions.extend(regionList)
        self.regionNames = []
        self.regionNames.extend(regionNames)
        self.meshes = []
        self.meshes.extend(meshList)
        self.meshNames = []
        self.meshNames.extend(meshNames)

    def addShellRegion(
        self,
        regType,
        keyPts,
        numEls,
        name=None,
        natSpaceCrds=None,
        elType="quad",
        meshMethod="free",
    ):
        if natSpaceCrds is None:
            natSpaceCrds = []
        self.shellRegions.append(
            ShellRegion(regType, keyPts, numEls, natSpaceCrds, elType, meshMethod)
        )
        if name is None:
            numReg = len(self.shellRegions)
            regName = "Sub-Region_" + str(numReg)
            self.regionNames.append(regName)
        else:
            self.regionNames.append(name)


    def getSurfaceMesh(self):
        allNds = []
        allEls = []
        elSetList = []
        numNds = 0
        numEls = 0
        regi = 0
        for reg in self.shellRegions:
            regMesh = reg.createShellMesh()
            allNds.extend(regMesh["nodes"])
            setList = []
            eli = 0
            for el in regMesh["elements"]:
                for i in range(0, 4):
                    if el[i] != -1:
                        el[i] = el[i] + numNds
                allEls.append(el)
                setList.append((eli + numEls))
                eli = eli + 1
            thisSet = {}
            thisSet["name"] = self.regionNames[regi]
            thisSet["labels"] = setList
            elSetList.append(thisSet)
            numNds = len(allNds)
            numEls = len(allEls)
            regi = regi + 1
        mshi = 0
        for msh in self.meshes:
            setList = []
            eli = 0
            for el in msh["elements"]:
                newEl = -1 * np.ones(4, dtype=int)
                for i in range(0, 4):
                    if el[i] != -1:
                        newEl[i] = el[i] + numNds
                allEls.append(newEl)
                setList.append((eli + numEls))
                eli = eli + 1
            thisSet = {}
            thisSet["name"] = self.meshNames[mshi]
            thisSet["labels"] = setList
            elSetList.append(thisSet)
            allNds.extend(msh["nodes"])
            numNds = len(allNds)
            numEls = len(allEls)
            mshi = mshi + 1
        mData = {}
        mData["nodes"] = np.array(allNds)
        mData["elements"] = np.array(allEls)
        mData = mt.mergeDuplicateNodes(mData)
        mData["sets"] = {}
        mData["sets"]["element"] = elSetList
        return mData
