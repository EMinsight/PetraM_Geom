import numpy as np

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshMeshActions')

from petram.phys.vtable import VtableElement, Vtable
from petram.mesh.gmsh_mesh_model import GmshMeshActionBase



data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Line#',
                                   default = "remaining", 
                                   tip = "Line ID" )),
        ('num_seg', VtableElement('num_seg', type='int',
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
        mesher.add('transfinite_line', gid, nseg=nseg,
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
                                   guilabel = 'Surface#',
                                   default = "", 
                                   tip = "Surface ID" )),
        ('edge1', VtableElement('edge1', type='string',
                                 guilabel = '1st edge',
                                 default = "", 
                                 tip = "1st Edge" )),
        ('edge2', VtableElement('edge2', type='string',
                                 guilabel = '2nd edge',
                                 default = "", 
                                 tip = "2ndp Edge" )),
        ('edge3', VtableElement('edge3', type='string',
                                 guilabel = '3rd edge',
                                 default = "", 
                                 tip = "3rd Edge" )),
        ('edge4', VtableElement('edge4', type='string',
                                 guilabel = '4th edge',
                                 default = "", 
                                 tip = "4th Edge" )),)
class TransfiniteSurface(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, e1, e2, e3, e4 = self.vt.make_value_or_expression(self)
        mesher.add('transfinite_surface', gid, edges = (e1,e2,e3,e4))

    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['face'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'face'
    
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
                                default_txt = '',                                 
                                default = 0.,
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size',
                                default_txt = '',                                
                                default = 0.0, 
                                tip = "CharacteristicLengthMin" )),
        ('resolution', VtableElement('resolution', type='int',
                                guilabel = 'Resolution',
                                default = 5., 
                                tip = "Edge Resolution" )),
        ('embed_s', VtableElement('embed_s', type='string',
                                   guilabel = 'Surface#',
                                   default = "", 
                                   tip = "Surface number" )),
        ('embed_l', VtableElement('embed_l', type='string',
                                   guilabel = 'Line#',
                                   default = "", 
                                   tip = "Line number" )),
        ('embed_p', VtableElement('embed_p', type='string',
                                   guilabel = 'Point#',
                                   default = "", 
                                   tip = "Point number" )),)


class FreeVolume(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        values = self.vt.make_value_or_expression(self)
        gid, clmax, clmin, res, embed_s, embed_l, embed_p = values
        mesher.add('freevolume', gid,
                   maxsize = clmax,
                   minsize = clmin,
                   resolution = res,
                   embed_s = embed_s,                   
                   embed_l=embed_l, embed_p=embed_p)
        
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
                                default_txt = '',                                 
                                default = 0.,
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size',
                                default = 0., 
                                tip = "CharacteristicLengthMin" )),
        ('resolution', VtableElement('resolution', type='int',
                                guilabel = 'Resolution',
                                default = 5., 
                                tip = "Edge Resolution" )),
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
        gid, clmax, clmin, res, embed_l, embed_p= self.vt.make_value_or_expression(self)
        mesher.add('freeface', gid,
                   maxsize = clmax,
                   minsize = clmin,
                   resolution = res,
                   embed_l=embed_l, embed_p=embed_p)
        
    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            ret['face'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'face'
    
    def get_embed(self):
        gid, clmax, clmin, embed_l, embed_p= self.vt.make_value_or_expression(self)
        ll = [str(x) for x in embed_l.split(',')]
        pp = [str(x) for x in embed_p.split(',')]                
        return [], ll, pp
            

data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Line#',
                                   default = "remaining", 
                                   tip = "Line number" )),
        ('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size',
                                default_txt = '',                                 
                                default = 0.0,
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size',
                                default_txt = '', default = 0.0,
                                tip = "CharacteristicLengthMin" )),
        ('resolution', VtableElement('resolution', type='int',
                                guilabel = 'Resolution',
                                default_txt = '5',                                     
                                default = 5., 
                                tip = "Edge Resolution" )),)

class FreeEdge(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, clmax, clmin, res = self.vt.make_value_or_expression(self)
        mesher.add('freeedge', gid,
                   maxsize = clmax,
                   minsize = clmin,
                   resolution = res)
        
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
                                  tip = "Coordinate transformatin " )),
        ('cp_cl', VtableElement('cp_cl', type='bool',
                                 guilabel = 'Copy CL',
                                 default = True,
                                 tip = "Copy Characteristic Length")), )



class CopyFace(GmshMeshActionBase):
    vt = Vtable(data)        
    def add_meshcommand(self, mesher):
        gid, src_id, transform, cp_cl = self.vt.make_value_or_expression(self)
        mesher.add('copymesh', gid, src_id, self.copyface_params,
                   mode = 'Surface')

        
    def check_master_slave(self, mesher):
        gid, src_id, transform, cp_cl = self.vt.make_value_or_expression(self)
        self.copyface_params = mesher.copymesh_face_trans_txt(gid, src_id, cp_cl)
        
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
    


rsdata =  (('geom_id', VtableElement('geom_id', type='string',
                                    guilabel = 'Surfaces',
                                    default = "",
                                    tip = "surfacess to be recombined")), 
           ('max_angle', VtableElement('max_angle', type='float',
                                guilabel = 'Max size)',
                                default = 45, 
                                tip = "Maximum differend of angle" )),)

class RecombineSurface(GmshMeshActionBase):
    vt = Vtable(rsdata)
    def add_meshcommand(self, mesher):
        gid, max_angle  = self.vt.make_value_or_expression(self)
        mesher.add('recombine_surface', gid, max_angle=max_angle)

edata =  (('ex_target', VtableElement('ex_target', type='string',
                                      guilabel = 'Volume',
                                      default = "",
                                      tip = "extrusion target")),
          ('dst_id', VtableElement('dst_id', type='string',
                                   guilabel = 'Surface# (To)',
                                   default = "", 
                                   tip = "Surface number" )),
          ('src_id', VtableElement('src_id', type='string',
                                  guilabel = 'Source # (From)',
                                  default = "", 
                                  tip = "Surface number" )),
          ('nlayer', VtableElement('nlayer', type='int',
                             guilabel = 'Num. Layers',
                             default = 5,
                             tip = "Number of Layers" )), )         

class MeshExtrude(GmshMeshActionBase):
    vt = Vtable(edata)
    def add_meshcommand(self, mesher):
        gid, dst_id, src_id, n_layers = self.vt.make_value_or_expression(self)
        tag = mesher.new_etg()
        mesher.add('meshextrude', tag, gid, dst_id, src_id, n_layers, self.extrude_params)
        
    def check_master_slave(self, mesher):
        gid, dst_id, src_id, n_layers = self.vt.make_value_or_expression(self)        
        self.extrude_params = mesher.extrude_trans_txt(gid, dst_id, src_id, True)
        
    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = self.element_selection_empty()
        try:
            dest = [int(x) for x in self.dst_id.split(',')]
            src  = [int(x) for x in self.src_id.split(',')]
            ret['face'] = dest + src
        except:
            pass
        return ret, 'face'
    
edata =  (('ex_target', VtableElement('ex_target', type='string',
                                      guilabel = 'Volume',
                                      default = "",
                                      tip = "extrusion target")),
          ('dst_id', VtableElement('dst_id', type='string',
                                   guilabel = 'Surface# (To)',
                                   default = "", 
                                   tip = "Surface number" )),
          ('src_id', VtableElement('src_id', type='string',
                                  guilabel = 'Source # (From)',
                                  default = "", 
                                  tip = "Surface number" )),
          ('nlayer', VtableElement('nlayer', type='int',
                             guilabel = 'Num. Layers',
                             default = 5,
                             tip = "Number of Layers" )), )         
    
MeshRevolve = MeshExtrude

    
    
