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
        mesher.add('transfinite_edge', gid, nseg=nseg,
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
        mesher.add('cl', gid, cl = cl)

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
        mesher.add('copyface',
                   src_id,                   
                   gid,
                   transfomr = transform,
                   copy_cl = cp_cl)
        
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
                             tip = "Number of Layers" )),
          ('mapper', VtableElement('mapper', type='string',
                                   guilabel = 'Transform Hint',
                                   default = "", 
                                   tip = "Coordinate transformatin (ax, an), (dx,dy,dz), ")),)
# Transform Hint (extrude)
#    d     : dx, dy, dz : 3 float
#    l1, l2,,,,: set of edges : end points determines d (N.I.)
def process_hint(text):
    try:
        values = [float(x) for x in text.split(',')]
        if len(values) == 3:
            return {'d': values}
    except:
        pass
    return {}
    
class ExtrudeMesh(GmshMeshActionBase):
    vt = Vtable(edata)
    def add_meshcommand(self, mesher):
        gid, dst_id, src_id, nlayers, hint = self.vt.make_value_or_expression(self)
        kwargs = process_hint(hint)
        mesher.add('extrude_face', gid, src_id, dst_id, nlayers=nlayers, **kwargs)

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

# Transform Hint (revolve)
#    ax an : ax_x, ax_y, ax_z, angle(deg): 4 float
#    l1, angle  : l1 direction of axis, angle (deg) (N.I.)
#    s1, angle  : normal to face s1, angle (deg) (N.I.)
def process_hint(text):
    try:
        values = [float(x) for x in text.split(',')]
        if len(values) == 4:
            return {'axan': (values[0:3], values[-1])}
    except:
        pass
    return {}
    
class RevolveMesh(GmshMeshActionBase):
    vt = Vtable(edata)
    def add_meshcommand(self, mesher):
        gid, dst_id, src_id, nlayers, hint = self.vt.make_value_or_expression(self)
        kwargs = process_hint(hint)
        mesher.add('revolve_face', gid, src_id, dst_id, nlayers=nlayers, **kwargs)        

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
    
    


    
    
