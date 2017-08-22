
import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomModel')

from petram.model import Model

from petram.namespace_mixin import NS_mixin

class GeomBase(Model, NS_mixin):
    def __init__(self, *args, **kwargs):
        super(GeomBase, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
    
    def onItemSelChanged(self, evt):
        '''
        GUI response when model object is selected in
        the dlg_edit_model
        '''
        viewer = evt.GetEventObject().GetTopLevelParent().GetParent()
        viewer.use_toolbar_palette('petram_geom', mode = '3D')
        viewer._view_mode = 'geom'
    
class MFEM_GeomRoot(GeomBase):
    can_delete = False
    has_2nd_panel = False
    def get_possible_child(self):
        from .gmsh_geom_model import GmshGeom
        return [GmshGeom]
    

        
    

