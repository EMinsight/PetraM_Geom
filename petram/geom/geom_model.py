'''

    Model Tree to stroe MFEM model parameters

'''
import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomModel')

from petram.model import Model

from petram.namespace_mixin import NS_mixin
from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

class Geom(Model, Vtable_mixin, NS_mixin):
    hide_ns_menu = True
    has_2nd_panel = False    
    def __init__(self, *args, **kwargs):
        super(Geom, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    def attribute_set(self, v):
        v = super(Geom, self).attribute_set(v)
        self.vt.attribute_set(v)
        return v
    
    def panel1_param(self):
        return self.vt.panel_param(self)
        
    def get_panel1_value(self):
        return self.vt.get_panel_value(self)

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        return self.vt.import_panel_value(self, v)

    def panel1_tip(self):
        return self.vt.panel_tip()


class GeomGroup(Model, NS_mixin):
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(GeomGroup, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
    
    def get_possible_child(self):
        from .primitives import Circle, Rect, Polygon, Extrude
        return [Circle, Rect, Polygon, Extrude]

    
class MFEM_GeomRoot(Model, NS_mixin):
    can_delete = False
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(MFEM_GeomRoot, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
    def get_possible_child(self):
        return [GeomGroup]
