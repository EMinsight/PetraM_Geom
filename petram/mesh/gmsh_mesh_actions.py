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
                                   default = 1.0, 
                                   tip = "Progression" )),
        ('bump', VtableElement('bump', type='float',
                                   guilabel = 'Bump',
                                   default = 1.0, 
                                   tip = "Bump" )),)

        
    
class TransfiniteLine(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, nseg, p, b = self.vt.make_value_or_expression(self)
        gid = self.eval_enitity_id(gid)
        
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
                                 guilabel = '1st corner',
                                 default = "", 
                                 tip = "1st Edge" )),
        ('edge2', VtableElement('edge2', type='string',
                                 guilabel = '2nd corner',
                                 default = "", 
                                 tip = "2ndp Edge" )),
        ('edge3', VtableElement('edge3', type='string',
                                 guilabel = '3rd corner',
                                 default = "", 
                                 tip = "3rd Edge" )),
        ('edge4', VtableElement('edge4', type='string',
                                 guilabel = '4th corner',
                                 default = "", 
                                 tip = "4th Edge" )),)
class TransfiniteSurface(GmshMeshActionBase):
    vt = Vtable(data)    
    def add_meshcommand(self, mesher):
        gid, e1, e2, e3, e4 = self.vt.make_value_or_expression(self)
        gid = self.eval_enitity_id(gid)
        
        c = [int(x)  for x in (e1,e2,e3,e4) if x.strip()!= '']
        mesher.add('transfinite_surface', gid, corner = c)

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
        gid = self.eval_enitity_id(gid)        
        mesher.add('cl', gid, cl)

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
        gid, embed_s, embed_l, embed_p = self.eval_enitity_id(gid, embed_s, embed_l, embed_p)
        
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
        gid, embed_l, embed_p = self.eval_enitity_id(gid, embed_l, embed_p)        
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
        gid  = self.eval_enitity_id(gid)
        
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

# Transform Hint (extrude)
#    d     : dx, dy, dz : 3 float
#    l1, l2,,,,: set of edges : end points determines d (N.I.)
def process_hint_ex(text):
    try:
        text = str(text)
        values = [x.strip() for x in text.split(',') if len(x.strip()) != 0]        
        values = [float(x) for x in values]
        if len(values) == 3:
            return {'axan': (values, 0.0)}
        else:
            assert False, "enter direction of translation"
    except:
        pass
    return {}
    
# Transform Hint (revolve)
#    ax an : ax_x, ax_y, ax_z, angle(deg): 4 float
#    l1, angle  : l1 direction of axis, angle (deg) (N.I.)
#    s1, angle  : normal to face s1, angle (deg) (N.I.)
def process_hint_rv(text):
    try:
        text = str(text)
        values = [x.strip() for x in text.split(',') if len(x.strip()) != 0]        
        values = [float(x) for x in values]
        if len(values) == 4:
            return {'axan': (values[0:3], values[-1])}
    except:
        pass
    return {}

data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Surface# (To)',
                                   default = "", 
                                   tip = "Surface number" )),
        ('src_id', VtableElement('src_id', type='string',
                                  guilabel = 'Source # (From)',
                                  default = "", 
                                  tip = "Surface number" )),
        ('mapper', VtableElement('mapper', type='string',
                                  guilabel = 'Transform Hint',                                 
                                  default = "", 
                                  tip = "Coordinate transformatin " )),
        ('cp_cl', VtableElement('cp_cl', type='bool',
                                 guilabel = 'Copy CL',
                                 default = True,
                                 tip = "Copy Characteristic Length")), )



class CopyFace(GmshMeshActionBase):
    vt = Vtable(data)        
    def add_meshcommand(self, mesher):
        gid, src_id, hint, cp_cl = self.vt.make_value_or_expression(self)
        gid  = self.eval_enitity_id(gid)
       
        kwargs = process_hint_ex(hint)
        kwargs['copy_cl'] = cp_cl
        print('adding here', src_id, gid,)
        mesher.add('copyface', src_id, gid,
                   **kwargs)
        
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
    
class CopyFaceRotate(GmshMeshActionBase):
    vt = Vtable(data)        
    def add_meshcommand(self, mesher):
        gid, src_id, hint, cp_cl = self.vt.make_value_or_expression(self)
        gid  = self.eval_enitity_id(gid)
                    
        kwargs = process_hint_rv(hint)
        kwargs['copy_cl'] = cp_cl
        kwargs['revolve'] = True
                    
        mesher.add('copyface', src_id, gid,
                   **kwargs)
        
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
    
merge_loc = [None, None, 36, {"col": 4, 
                              "labels":["0D","1D","2D","3D"]}]
class MergeText(GmshMeshActionBase):
    vt = Vtable(tuple())
    def panel1_param(self):
        ll = super(MergeText, self).panel1_param()
        
        ll.extend([["Text", "", 235, {"nlines":15}],
                    merge_loc,])

        return ll
      
    def attribute_set(self, v):
        v = super(MergeText, self).attribute_set(v)
        v["merge_txt"] = ""
        v["merge_dim"] = [True]+[False]*3
        return v
        
    def get_panel1_value(self):
        v = super(MergeText, self).get_panel1_value()
        v.extend([self.merge_txt, self.merge_dim])
        return  v

    def preprocess_params(self, engine):
        return

    def import_panel1_value(self, v):
        super(MergeText, self).import_panel1_value(v[:-2])
        self.merge_txt = str(v[-2])
        self.merge_dim = [x[1] for x in v[-1]]

    def panel1_tip(self):
        return [None, None, None, None]
    
    def add_meshcommand(self, mesher):
        self.vt.preprocess_params(self)                        
        mesher.add('mergetxt', text = self.merge_txt, dim = self.merge_dim)

data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Surfaces',
                                   default = "", 
                                   tip = "Entity number" )),)
class CompoundSurface(GmshMeshActionBase):
    vt = Vtable(data)
    def add_meshcommand(self, mesher):
        gid = self.vt.make_value_or_expression(self)[0]        
        gid  = self.eval_enitity_id(gid)

        # generate something like... Compound Surface{1, 5, 10};
        text = "Compound Surface{ " + gid + "};"        
        mesher.add('mergetxt', text = text , dim = [True, False, False, False])

data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Curves',
                                   default = "", 
                                   tip = "Entity number" )),)
        
class CompoundCurve(GmshMeshActionBase):
    vt = Vtable(data)
    def add_meshcommand(self, mesher):
        gid = self.vt.make_value_or_expression(self)[0]                
        gid  = self.eval_enitity_id(gid)

        # generate something like... Compound Curve{1, 5, 10};        
        text = "Compound Curve{ " + gid + "};"
        mesher.add('mergetxt', text = text , dim = [True, False, False, False])
        
        
rsdata =  (('geom_id', VtableElement('geom_id', type='string',
                                    guilabel = 'Surfaces',
                                    default = "",
                                    tip = "surfacess to be recombined")), )
#           ('max_angle', VtableElement('max_angle', type='float',
#                                guilabel = 'Max size)',
#                                default = 45, 
#                                tip = "Maximum differend of angle" )),)

class RecombineSurface(GmshMeshActionBase):
    vt = Vtable(rsdata)
    def add_meshcommand(self, mesher):
        gid = self.vt.make_value_or_expression(self)[0]
        gid  = self.eval_enitity_id(gid)
        
        mesher.add('recombine_surface', gid)

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
    
class ExtrudeMesh(GmshMeshActionBase):
    vt = Vtable(edata)
    def add_meshcommand(self, mesher):
        gid, dst_id, src_id, nlayers, hint = self.vt.make_value_or_expression(self)
        gid, dst_id, src_id  = self.eval_enitity_id(gid, dst_id, src_id)
        
        kwargs = process_hint_ex(hint)
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

    
class RevolveMesh(GmshMeshActionBase):
    vt = Vtable(edata)
    def add_meshcommand(self, mesher):
        gid, dst_id, src_id, nlayers, hint = self.vt.make_value_or_expression(self)
        gid, dst_id, src_id  = self.eval_enitity_id(gid, dst_id, src_id)
        
        kwargs = process_hint_rv(hint)
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
    
    


    
    
