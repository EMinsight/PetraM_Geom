
import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomModel')

from petram.model import Model

from petram.namespace_mixin import NS_mixin

class MFEM_GeomRoot(Model, NS_mixin):
    can_delete = False
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(MFEM_GeomRoot, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
    def get_possible_child(self):
        from .gmsh_geom_model import GmshGeom
        return [GmshGeom]

