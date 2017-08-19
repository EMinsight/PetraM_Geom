import numpy as np

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshPrimitives')

from petram.phys.vtable import VtableElement, Vtable
from petram.geom.gmsh_geom_model import GmshPrimitiveBase as GeomPB
from petram.geom.gmsh_geom_model import get_geom_key


cdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('ax1', VtableElement('ax1', type='float',
                                   guilabel = 'axis1',
                                   suffix =('x', 'y', 'z'),
                                   default = [1, 0, 0], 
                                   tip = "axis 1" )),
          ('ax2', VtableElement('ax2', type='float',
                                 guilabel = 'axis2',
                                 suffix =('x', 'y', 'z'),
                                 default = [0, 1, 0], 
                                 tip = "axis 2" )),
          ('radius', VtableElement('radius', type='float',
                                   guilabel = 'r',
                                   default = 1.0, 
                                   tip = "radius" )),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                                   tip = "characteristc length from point" )),)

class Circle(GeomPB):
    vt = Vtable(cdata)
    def build_geom(self, geom, objs):
        center, ax1, ax2, radius, lcar = self.vt.make_value_or_expression(self)
        a1 = np.array(ax1);  a2 = np.array(ax2)
        a2 = np.cross(np.cross(a1, a2), a1)
        a1 = a1/np.sqrt(np.sum(a1**2))*radius
        a2 = a2/np.sqrt(np.sum(a2**2))*radius                      

        c =np.array(center)
        p1 = geom.add_point(c+a1, lcar)
        p2 = geom.add_point(c+a2, lcar)
        p3 = geom.add_point(c-a1, lcar)
        p4 = geom.add_point(c-a2, lcar)                      
        pc = geom.add_point(c, lcar)
        ca1 = geom.add_circle_arc(p1, pc, p2)
        ca2 = geom.add_circle_arc(p2, pc, p3)
        ca3 = geom.add_circle_arc(p3, pc, p4)
        ca4 = geom.add_circle_arc(p4, pc, p1)
        ll1 = geom.add_line_loop([ca1, ca2, ca3, ca4])
        ps1 = geom.add_plane_surface(ll1)
        newkey = objs.addobj(ps1, 'ps')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]

rdata =  (('corner', VtableElement('corner', type='float',
                             guilabel = 'Corner',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('edge1', VtableElement('edge1', type='float',
                             guilabel = 'Edge(1)',
                             suffix =('x', 'y', 'z'),
                             default = [1,0,0],
                             tip = "Edge of rectangle" )),
          ('edge2', VtableElement('edge2', type='float',
                             guilabel = 'Edge(2)',
                             suffix =('x', 'y', 'z'),
                             default = [0,1,0],
                             tip = "Edge of rectangle" )),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.1, 
                              tip = "characteristc length from point" )),)

class Rect(GeomPB):
    vt = Vtable(rdata)
    
    def build_geom(self, geom, objs):
        c1,  e1,  e2,  lcar = self.vt.make_value_or_expression(self)
        c1 = np.array(c1);
        e1 = np.array(e1);        e2 = np.array(e2);
        p1 = geom.add_point(c1, lcar)
        p2 = geom.add_point(c1+e1, lcar)
        p3 = geom.add_point(c1+e1+e2, lcar)
        p4 = geom.add_point(c1+e2, lcar)
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p3)
        l3 = geom.add_line(p3, p4)
        l4 = geom.add_line(p4, p1)        
        ll1 = geom.add_line_loop([l1, l2, l3, l4])
        rec1 = geom.add_plane_surface(ll1)
        newkey = objs.addobj(rec1, 'rec')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]

    
pdata =  (('xarr', VtableElement('xarr', type='array',
                              guilabel = 'X',
                              default = '0.0',
                              tip = "X" )),
          ('yarr', VtableElement('yarr', type='array',
                              guilabel = 'Y',
                              default = '0.0',
                              tip = "Y" )),
          ('zarr', VtableElement('zarr', type='array',
                              guilabel = 'Z',
                              default = '0.0',
                              tip = "Z" )),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.1, 
                              tip = "characteristc length from point" )),)
class Polygon(GeomPB):
    vt = Vtable(pdata)
    def build_geom(self, geom, objs):
        xarr, yarr, zarr,  lcar = self.vt.make_value_or_expression(self)
        if len(xarr) < 2: return
        print type(xarr)
        try:
           pos = np.vstack((xarr, yarr, zarr)).transpose()
        except:
           print("can not make proper input array")
           return
        # check if data is already closed...
        if np.abs(np.sum((pos[0] - pos[-1])**2)) < 1e-17:
            pos = pos[:-1]
        poly = geom.add_polygon(pos, lcar = lcar)

        # apparently I should use this object (poly.surface)...?
        newkey = objs.addobj(poly.surface, 'pol')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]

edata =  (('ex_target', VtableElement('ex_target', type='string',
                                      guilabel = 'Target',
                                      default = "",
                                      tip = "extrusion target")),
          ('paxis', VtableElement('paxis', type='float',
                             guilabel = 'Point on Axis',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "point on axis" )),
          ('taxis', VtableElement('taxis', type='float',
                                   guilabel = 'Translation Axis',
                                   suffix =('x', 'y', 'z'),
                                   default = [0, 0, 1],
                                   tip = "translation axis" )),
          ('ex_len', VtableElement('exlen', type='float',
                             guilabel = 'Length',
                             default = 1.0,
                             tip = "Extrusion length" )), )         

class Extrude(GeomPB):    
    vt = Vtable(edata)

    def build_geom(self, geom, objs):
        targets, pax, tax, len = self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]

        tax = tax/np.sqrt(np.sum(np.array(tax)**2))*len          
        newkeys = []
        for t in targets:
             if not t in objs:
                 assert False, t + " does not exist"
             ret = geom.extrude(objs[t],
                          translation_axis=tax,
                          #rotation_axis=rax,
                          point_on_axis=pax)

             newkeys.append(objs.addobj(ret[0], t))
             newkeys.append(objs.addobj(ret[1], 'ex'))             
             #for o in ret[2:]:
             #   newkeys.append(objs.addobj(o,  get_geom_key(o))
                               
        self._objkeys = objs.keys()
        self._newobjs = newkeys

edata =  (('ex_target', VtableElement('ex_target', type='string',
                                      guilabel = 'Target',
                                      default = "",
                                      tip = "extrusion target")),
          ('paxis', VtableElement('paxis', type='float',
                             guilabel = 'Point on Axis',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "point on axis" )),
          ('raxis', VtableElement('raxis', type='float',
                                   guilabel = 'Translation Axis',
                                   suffix =('x', 'y', 'z'),
                                   default = [0, 0, 1],
                                   tip = "translation axis" )),
          ('angle', VtableElement('angel', type='float',
                             guilabel = 'angle',
                             default = 180.,
                             tip = "Extrusion length" )), )         
        
class Revolve(GeomPB):    
    vt = Vtable(edata)

    def build_geom(self, geom, objs):
        targets, pax, rax, angle = self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        for t in targets:
             if not t in objs:
                 assert False, t + " does not exist"
             ret = geom.extrude(objs[t],
                                rotation_axis=rax,
                                point_on_axis=pax,
                                angle = angle*np.pi/180.)

             newkeys.append(objs.addobj(ret[0], t))
             newkeys.append(objs.addobj(ret[1], 'ex'))             
             #for o in ret[2:]:
             #   newkeys.append(objs.addobj(o,  get_geom_key(o))
                               
        self._objkeys = objs.keys()
        self._newobjs = newkeys



class GeomPB_Bool(GeomPB):
    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        self.vt.attribute_set(v)
        v["delete_input"] = True
        return v
    
    def panel1_param(self):
        ll = self.vt.panel_param(self)
        ll.append(["Delete",
                    self.delete_input,  3, {"text":""}])
        return ll
        
    def get_panel1_value(self):
        return list(self.vt.get_panel_value(self)) + [self.delete_input]

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        self.vt.import_panel_value(self, v[:-1])
        self.delete_input = v[-1]

    def panel1_tip(self):
        return list(self.vt.panel_tip()) + ['delete input objects']


ddata =  (('objplus', VtableElement('objplus', type='string',
                                      guilabel = '+',
                                      default = "",
                                      tip = "added objects")),
          ('objminus', VtableElement('objminus', type='string',
                                      guilabel = '-',
                                      default = "",
                                      tip = "objects to be subtracted")),)
                
class Difference(GeomPB_Bool):    
    vt = Vtable(ddata)
    def build_geom(self, geom, objs):
        tp, tm  = self.vt.make_value_or_expression(self)
        tp = [x.strip() for x in tp.split(',')]
        tm = [x.strip() for x in tm.split(',')]          

        input_entity = [objs[x] for x in tp]
        tool_entity  = [objs[x] for x in tm]
        ret = geom.boolean_difference(
                          input_entity,
                          tool_entity,
                          delete = self.delete_input)
        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'diff'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o,  get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
udata =  (('objplus', VtableElement('obj1', type='string',
                                      guilabel = 'input',
                                      default = "",
                                      tip = "objects")),)

class Union(GeomPB_Bool):    
    vt = Vtable(udata)
    def build_geom(self, geom, objs):
        tp  = self.vt.make_value_or_expression(self)
        tp = [x.strip() for x in tp.split(',')]
        if len(tp) < 2: return

        input_entity = [objs[x] for x in tp[:1]]
        tool_entity  = [objs[x] for x in tp[1:]]
        ret = geom.boolean_union(
                          input_entity,
                          tool_entity,
                          delete = self.delete_input)
        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'diff'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o,  get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
class Intersection(GeomPB_Bool):    
    vt = Vtable(udata)
    def build_geom(self, geom, objs):
        tp  = self.vt.make_value_or_expression(self)
        tp = [x.strip() for x in tp.split(',')]
        if len(tp) < 2: return

        input_entity = [objs[x] for x in tp[:1]]
        tool_entity  = [objs[x] for x in tp[1:]]
        ret = geom.boolean_intersection(
                          input_entity,
                          tool_entity,
                          delete = self.delete_input)
        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'diff'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o,  get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys

class Fragments(GeomPB_Bool):    
    vt = Vtable(ddata)
    def build_geom(self, geom, objs):
        tp  = self.vt.make_value_or_expression(self)
        tp = [x.strip() for x in tp.split(',')]
        if len(tp) < 2: return

        input_entity = [objs[x] for x in tp[:1]]
        tool_entity  = [objs[x] for x in tp[1:]]
        ret = geom.boolean_fragments(
                          input_entity,
                          tool_entity,
                          delete = self.delete_input)
        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'diff'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o,  get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
          

        
    
