import numpy as np

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshMeshActions')

from petram.phys.vtable import VtableElement, Vtable
from petram.mesh.gmsh_mesh_model import GmshMeshActionBase



data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Line#',
                                   default = "remaining", 
                                   tip = "Line ID" )),
        ('num_seg', VtableElement('radius', type='float',
                                   guilabel = 'Number of segments',
                                   default = 5, 
                                   tip = "Number of segments" )),
        ('progression', VtableElement('progression', type='float',
                                   guilabel = 'Progression',
                                   default = 0, 
                                   tip = "Progression" )),
        ('bump', VtableElement('bump', type='float',
                                   guilabel = 'Bump',
                                   default = 0, 
                                   tip = "Bump" )),)

        
    
class TransfiniteLine(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, nseg, p, b = self.vt.make_value_or_expression(self)
        mesher.add('transfinite', gid, mode = 'Line', nseg=nseg,
                   progression = p,  bump = b)

    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['edge'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'edge'
    
data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Point#',
                                   default = "remaining", 
                                   tip = "Point ID" )),
        ('cl', VtableElement('cl', type='float',
                                guilabel = 'Size)',
                                default = 1.0, 
                                tip = "CharacteristicLength" )))

class CharacteristicLength(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, cl= self.vt.make_value_or_expression(self)
        mesher.add('characteristiclength', gid, cl = cl)

    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['point'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'point'


data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Element#',
                                   default = "remaining", 
                                   tip = "Element ID" )),
        ('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size)',
                                default = 1.0, 
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size',
                                default = 1.0, 
                                tip = "CharacteristicLengthMin" )),
        ('embed_s', VtableElement('embed_s', type='string',
                                   guilabel = 'Surface#',
                                   default = "", 
                                   tip = "Surface number" )),)


class FreeVolume(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, clmax, clmin, embed_s = self.vt.make_value_or_expression(self)
        mesher.add('freemesh', gid, clmax = clmax, clmin = clmin,
                   mode = 'Volume', embed_s = embed_s)

    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['volume'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'volume'
    
data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Element#',
                                   default = "remaining", 
                                   tip = "Element ID" )),
        ('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size)',
                                default = 0.0, 
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size',
                                default = 0.0, 
                                tip = "CharacteristicLengthMin" )),
        ('embed_l', VtableElement('embed_l', type='string',
                                   guilabel = 'Line#',
                                   default = "", 
                                   tip = "Line number" )),
        ('embed_p', VtableElement('embed_p', type='string',
                                   guilabel = 'Point#',
                                   default = "", 
                                   tip = "Point number" )),)

class FreeFace(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, clmax, clmin, embed_l, embed_p= self.vt.make_value_or_expression(self)
        mesher.add('freemesh', gid, clmax = clmax, clmin = clmin,
                   mode = 'Surface', embed_l=embed_l, embed_p=embed_p)
        
    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['face'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'face'

data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Line#',
                                   default = "remaining", 
                                   tip = "Line number" )),
        ('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size)',
                                default = 0.0, 
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size',
                                default = 0.0, 
                                tip = "CharacteristicLengthMin" )),)
class FreeEdge(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, clmax, clmin = self.vt.make_value_or_expression(self)
        mesher.add('freemesh', gid, clmax = clmax, clmin = clmin,
                   mode = 'Line')
        
    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['edge'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'edge'
    
data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Surface# (To)',
                                   default = "", 
                                   tip = "Surface number" )),
        ('src_id', VtableElement('src_id', type='string',
                                  guilabel = 'Source# (From)',
                                  default = "", 
                                  tip = "Surface number" )),)
class Rotate(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, src_id = self.vt.make_value_or_expression(self)
        mesher.add('rotate', gid, src=src_id, transform=self.name())
        
data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Point# (To)',
                                   default = "", 
                                   tip = "Point number" )),
        ('src_id', VtableElement('src_id', type='string',
                                  guilabel = 'Point# (From)',
                                  default = "", 
                                  tip = "Point number" )),)

class Translate(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, src_id = self.vt.make_value_or_expression(self)
        mesher.add('translate', gid, src=src_id, transform=self.name())
        
data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Surface# (To)',
                                   default = "", 
                                   tip = "Surface number" )),
        ('src_id', VtableElement('src_id', type='string',
                                  guilabel = 'Source # (From)',
                                  default = "", 
                                  tip = "Surface number" )),
        ('mapper', VtableElement('mapper', type='string',
                                  guilabel = 'Transform',
                                  default = "", 
                                  tip = "Coordinate transformatin " )),)

class CopyFace(GmshMeshActionBase):
    vt = Vtable(data)        
    def add_meshcommand(self, mesher):
        gid, src_id, transform = self.vt.make_value_or_expression(self)
        etg1 = mesher.new_etg()
        etg2 = mesher.new_etg()
        mesher.add('copymesh', gid, src_id, mode = 'Surface',
                   transform = transform, etg = [etg1, etg2])
        
    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            dest = [int(x) for x in self.geom_id.split(',')]
            src  = [int(x) for x in self.src_id.split(',')]
            ret['face'] = dest + src
        except:
            pass
        return ret, 'face'
