
from aeroelast.models.blade.numad.objects.stackdb import StackDatabase


class MaterialDatabase(dict):
    def __init__(self):
        pass

    def generate(self, materials: dict, stackdb: StackDatabase):
        """Adds material and composites information to MatDB

        Parameters
        ----------
        materials : dict
        stackdb : StackDatabase

        Returns
        -------
        Self
        """
        mm_to_m = 0.001
        stacks = stackdb.stacks
        swstacks = stackdb.swstacks
        if stacks is None or swstacks is None:
            raise RuntimeError(
                "MaterialDatabase.generate() requires a generated StackDatabase"
            )
        n_webs, n_stations = swstacks.shape

        # prepare material database ==========================================

        # add base materials
        for mat_name in materials:
            cur_entry = MaterialDatabaseEntry()
            cur_material = materials[mat_name]
            cur_entry.name = cur_material.name
            cur_entry.type = cur_material.type
            cur_entry.ex = cur_material.ex
            cur_entry.ey = cur_material.ey
            cur_entry.ez = cur_material.ez
            cur_entry.gxy = cur_material.gxy
            cur_entry.gyz = cur_material.gyz
            cur_entry.gxz = cur_material.gxz
            if cur_entry.type == "isotropic":
                cur_entry.nuxy = cur_material.prxy
            else:
                cur_entry.prxy = cur_material.prxy
                cur_entry.pryz = cur_material.pryz
                cur_entry.prxz = cur_material.prxz
            cur_entry.dens = cur_material.density
            cur_entry.reference = cur_material.reference
            self[mat_name] = cur_entry

        # add component stacks
        flat_stacks = stacks.flatten("F")
        for k_stack in range(flat_stacks.size):
            cur_stack = flat_stacks[k_stack]
            cur_entry = MaterialDatabaseEntry()
            cur_entry.name = cur_stack.name
            cur_entry.type = "composite"
            cur_entry.reference = "Reference text"
            cur_entry.thicknessType = "Constant"
            cur_entry.uniqueLayers = len(cur_stack.plygroups)
            cur_entry.symmetryType = "none"
            cur_entry.layer = [None] * cur_entry.uniqueLayers
            for j in range(cur_entry.uniqueLayers):
                cur_layer = Layer()
                plygroup = cur_stack.plygroups[j]
                matid = plygroup.materialid
                cur_layer.layer_name = self[matid].name
                cur_layer.thicknessA = mm_to_m * plygroup.thickness
                cur_layer.thicknessB = cur_layer.thicknessA
                cur_layer.quantity = plygroup.nPlies
                cur_layer.theta = plygroup.angle
                cur_entry.layer[j] = cur_layer
            self[cur_entry.name] = cur_entry

        # add shearweb stacks
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                cur_stack = swstacks[k_web, k_stat]
                cur_entry = MaterialDatabaseEntry()
                cur_entry.name = cur_stack.name
                cur_entry.type = "composite"
                cur_entry.reference = "Reference text"
                cur_entry.thicknessType = "Constant"
                try:
                    cur_entry.uniqueLayers = len(cur_stack.plygroups)
                except TypeError:
                    cur_entry.uniqueLayers = 0
                cur_entry.symmetryType = "none"
                cur_entry.layer = [None] * cur_entry.uniqueLayers
                for j in range(cur_entry.uniqueLayers):
                    cur_plygroup = cur_stack.plygroups[j]
                    cur_layer = Layer()
                    matid = cur_plygroup.materialid
                    cur_layer.layer_name = self[matid].name
                    cur_layer.thicknessA = mm_to_m * cur_plygroup.thickness
                    cur_layer.thicknessB = cur_layer.thicknessA
                    cur_layer.quantity = cur_plygroup.nPlies
                    cur_layer.theta = cur_plygroup.angle
                    cur_entry.layer[j] = cur_layer
                self[cur_entry.name] = cur_entry

        # shearweb information from NuMAD v1 is formatted in a specific
        # way, recreating that here
        # recreating data.shearweb ====================================
        # NOTE: do we care about this anymore? -kb
        # ctr = 0
        # self.shearweb = []
        # for k_web in range(n_webs):
        #     ind = keypoints.web_indices[k_web]
        #     for k_stat in range(n_stations):
        #         if swstacks[k_web, k_stat].plygroups:
        #             cur_sw = ShearWeb()
        #             cur_sw.Material = swstacks[k_web, k_stat].name
        #             cur_sw.BeginStation = swstacks[k_web, k_stat].indices[0]  # =k
        #             cur_sw.EndStation = swstacks[k_web, k_stat].indices[1]  # =k+1
        #             cur_sw.Corner = [
        #                 ind[1] - 1,
        #                 ind[0] - 1,
        #                 ind[0] - 1,
        #                 ind[1] - 1,
        #             ]  # dp number is offset by 1 in NuMAD v1
        #             self.shearweb.append(cur_sw)
        #             ctr += 1
        return self


class MaterialDatabaseEntry:
    """A simple class to organize the attributes of a material"""

    def __init__(self):
        self.type: str | None = None
        self.name: str | None = None
        self.reference: str | None = None
        self.dens: float | None = None
        self.nuxy: float | None = None
        self.ex: float | None = None
        self.ey: float | None = None
        self.ez: float | None = None
        self.gxy: float | None = None
        self.gyz: float | None = None
        self.gxz: float | None = None
        self.prxy: float | None = None
        self.pryz: float | None = None
        self.prxz: float | None = None
        self.xten: float | None = None
        self.xcmp: float | None = None
        self.yten: float | None = None
        self.ycmp: float | None = None
        self.zten: float | None = None
        self.zcmp: float | None = None
        self.xy: float | None = None
        self.yz: float | None = None
        self.xz: float | None = None
        self.xycp: float | None = None
        self.yzcp: float | None = None
        self.xzcp: float | None = None
        self.xzit: float | None = None
        self.xzic: float | None = None
        self.yzit: float | None = None
        self.yzic: float | None = None
        self.g1g2: float | None = None
        self.etal: float | None = None
        self.etat: float | None = None
        self.alp0: float | None = None
        self.thicknessType: str | None = None
        self.uniqueLayers: int | None = None
        self.symmetryType: str | None = None
        self.layer: list | None = None


class Layer:
    """A simple class to organize the attributes of a material layer.

    Attributes
    ----------
    layer_name : str
    thicknessA : float
    thicknessB :float
    quantity : int
    theta : float
    """

    def __init__(self):
        self.layer_name: str | None = None
        self.thicknessA: float | None = None
        self.thicknessB: float | None = None
        self.quantity: int | None = None
        self.theta: float | None = None
