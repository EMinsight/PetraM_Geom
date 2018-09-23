import numpy as np
from collections import defaultdict

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshPrimitives')

from petram.phys.vtable import VtableElement, Vtable
from petram.geom.gmsh_geom_model import GmshPrimitiveBase as GeomPB
from petram.geom.gmsh_geom_model import get_geom_key

try:
  import pygmsh
  class Geometry(pygmsh.Geometry):
    def __init__(self, *args, **kwargs):
        self._point_loc = {}
        self._point_mask = []
        super(Geometry, self).__init__(*args, **kwargs)

    def add_point(self, *args, **kwargs):
        mask = kwargs.pop("mask", True)      
        pt = tuple(args[0])
        if not pt in self._point_loc:
            obj = super(Geometry, self).add_point(*args, **kwargs)
            self._point_loc[pt] = obj
        else:
            obj = self._point_loc[pt]

        if mask : self._point_mask.append(obj)
        return obj
      
    def add_surface_filling(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def add_sphere(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def add_wedge(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def add_torus(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def add_cone(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
        
    def boolean_union(self, *args, **kwargs):
        a = kwargs.pop("removeObject", False)
        b = kwargs.pop("removeTool", False)
        kwargs["delete"] = a or b
        return super(Geometry, self).boolean_union(*args, **kwargs)
      
    def boolean_difference(self, *args, **kwargs):
        a = kwargs.pop("removeObject", False)
        b = kwargs.pop("removeTool", False)
        kwargs["delete"] = a or b
        return super(Geometry, self).boolean_difference(*args, **kwargs)
      
    def boolean_intersection(self, *args, **kwargs):
        a = kwargs.pop("removeObject", False)
        b = kwargs.pop("removeTool", False)
        kwargs["delete"] = a or b
        return super(Geometry, self).boolean_intersection(*args, **kwargs)
      
    def boolean_fragments(self, *args, **kwargs):
        a = kwargs.pop("removeObject", False)
        b = kwargs.pop("removeTool", False)
        kwargs["delete"] = a or b
        return super(Geometry, self).boolean_fragments(*args, **kwargs)
      
    def copy(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def remove(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def rotate(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def translate(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def dilate(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
    def symmetrize(self, *args, **kwargs):
        assert False, "Not implemented for gmsh3. Use gmsh new API"
        
  has_gmsh = True
except:
  has_gmsh = False

'''
object created from givn coordinates
   Point
   Circle
   Rect
   Line
   Polygon

object created from given entities
   Spline -> Points to Spline
   CreateLine -> Points to multiple Lines

   LineLoop -> Lines to LineLoop
   CreateSurface -> Lines to PlaneSurface

   SurfaceLoop -> Surfacess to SurfaceLoop
   CreateVolume ->  Surfacess to Volume

object manipulation
   Extrude
   Revolve
   Recombine Surface

boolean operation
   Union
   Intersection
   Difference
   Fragments

'''
invalid_pdata = (('NotUsedValue', VtableElement('NotUsedValue', type='array',
                                  guilabel = 'not_implemented',
                                  default = '0.0',
                                  tip = "This panel is not implemented" )),)

pdata = (('xarr', VtableElement('xarr', type='array',
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
                                   default = 0.0, 
                              tip = "characteristc length from point" )),)
class Point(GeomPB):
    vt = Vtable(pdata)
    def build_geom(self, geom, objs):
        xarr, yarr, zarr,  lcar = self.vt.make_value_or_expression(self)
        try:
           pos = np.vstack((xarr, yarr, zarr)).transpose()
        except:
           print("can not make proper input array")
           return
        PTs = [geom.add_point(p, lcar=lcar) for p in pos]
        # apparently I should use this object (poly.surface)...?
        self._newobjs = []        
        for p in PTs:
           newkey = objs.addobj(p, 'pt')
           self._newobjs.append(newkey)
        self._objkeys = objs.keys()

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
                                   default = 0.0, 
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
          ('edge3', VtableElement('edge3', type='float',
                                  guilabel = 'Edge(3)',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,1],
                             tip = "Edge of rectangle" )),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                              tip = "characteristc length from point" )),)

class Box(GeomPB):
    vt = Vtable(rdata)
    
    def build_geom(self, geom, objs):
        c1,  e1,  e2,  e3, lcar = self.vt.make_value_or_expression(self)
        c1 = np.array(c1);
        e1 = np.array(e1);        e2 = np.array(e2);
        p1 = geom.add_point(c1, lcar)
        p2 = geom.add_point(c1+e1, lcar)
        p3 = geom.add_point(c1+e2, lcar)
        p4 = geom.add_point(c1+e3, lcar)
        p5 = geom.add_point(c1+e1+e2, lcar)        
        p6 = geom.add_point(c1+e2+e3, lcar)
        p7 = geom.add_point(c1+e3+e1, lcar)
        p8 = geom.add_point(c1+e3+e2+e1, lcar)
        
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p5)
        l3 = geom.add_line(p5, p3)
        l4 = geom.add_line(p3, p1)
        
        l5 = geom.add_line(p1, p4)
        l6 = geom.add_line(p2, p7)
        l7 = geom.add_line(p5, p8)
        l8 = geom.add_line(p3, p6)
        
        l9  = geom.add_line(p4, p7)
        l10 = geom.add_line(p7, p8)
        l11 = geom.add_line(p8, p6)
        l12 = geom.add_line(p6, p4)        
        
        ll1 = geom.add_line_loop([l1, l2, l3, l4])
        ll2 = geom.add_line_loop([l5, l9, l6, l1])
        ll3 = geom.add_line_loop([l6, l10, l7, l2])
        ll4 = geom.add_line_loop([l7, l11, l8, l3])
        ll5 = geom.add_line_loop([l8, l12, l5, l4])
        ll6 = geom.add_line_loop([l9, l10, l11, l12])
        
        rec1 = geom.add_plane_surface(ll1)
        rec2 = geom.add_plane_surface(ll2)
        rec3 = geom.add_plane_surface(ll3)
        rec4 = geom.add_plane_surface(ll4)
        rec5 = geom.add_plane_surface(ll5)
        rec5 = geom.add_plane_surface(ll6)

        sl = geom.add_surface_loop([ll1, ll2, ll3, ll4, ll5, ll6])
        v1 = geom.add_volume(sl)
        
        newkey = objs.addobj(v1, 'bx')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]


vtdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('x_radius', VtableElement('x_radius', type='float',
                                     guilabel = 'radius (X)',
                                     default = 1.0,
                                     tip = "radius in X direction" )),
          ('y_radius', VtableElement('y_radius', type='float',
                                     guilabel = 'radius (Y)',
                                     default = 1.0,
                                     tip = "radius in Y direction" )),
          ('z_rarius', VtableElement('z_radius', type='float',
                                     guilabel = 'radius (Z)',
                                     default = 1.0,                                
                                     tip = "radius in Z direction" )),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                                   tip = "characteristc length from point" )),)
                 
class Ball(GeomPB):
    vt = Vtable(vtdata)
    def build_geom(self, geom, objs):
        x0,  l1,  l2,  l3, lcar = self.vt.make_value_or_expression(self)

        # Add points.
        radii = [l1, l2, l3]
        rr = min(radii)
        v1 = geom.add_sphere(x0[0], x0[1], x0[2], rr)
        if (l1/rr != 1.0 or l2/rr != 1.0 or l3/rr != 1.0):
            geom.dilate([v1], x0[0], x0[1], x0[2], l1/rr, l2/rr, l3/rr)
        newkey = objs.addobj(v1, 'bl')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]

        '''
        print(radii)
        p = [ geom.add_point(x0, lcar=lcar),
              geom.add_point([x0[0]+radii[0], x0[1], x0[2]], lcar=lcar),
              geom.add_point([x0[0], x0[1]+radii[1], x0[2]], lcar=lcar),
              geom.add_point([x0[0], x0[1], x0[2]+radii[2]], lcar=lcar),
              geom.add_point([x0[0]-radii[0], x0[1], x0[2]], lcar=lcar),
              geom.add_point([x0[0], x0[1]-radii[1], x0[2]], lcar=lcar),
              geom.add_point([x0[0], x0[1], x0[2]-radii[2]], lcar=lcar), ]

        c = [geom.add_ellipse_arc(p[1], p[0], p[6]),
             geom.add_ellipse_arc(p[6], p[0], p[4]),
             geom.add_ellipse_arc(p[4], p[0], p[3]),
             geom.add_ellipse_arc(p[3], p[0], p[1]),
             geom.add_ellipse_arc(p[1], p[0], p[2]),
             geom.add_ellipse_arc(p[2], p[0], p[4]),
             geom.add_ellipse_arc(p[4], p[0], p[5]),
             geom.add_ellipse_arc(p[5], p[0], p[1]),
             geom.add_ellipse_arc(p[6], p[0], p[2]),
             geom.add_ellipse_arc(p[2], p[0], p[3]),
             geom.add_ellipse_arc(p[3], p[0], p[5]),
             geom.add_ellipse_arc(p[5], p[0], p[6]), ]
        ll = [geom.add_line_loop([c[4], c[9], c[3]]),
              geom.add_line_loop([c[8], -c[4], c[0]]),
              geom.add_line_loop([-c[9], c[5], c[2]]),
              geom.add_line_loop([-c[5], -c[8], c[1]]),
              geom.add_line_loop([c[7], -c[3], c[10]]),
              geom.add_line_loop([c[11], -c[7], -c[0]]),
              geom.add_line_loop([-c[10], -c[2], c[6]]),
              geom.add_line_loop([-c[1], -c[6], -c[11]]), ]
        s = [geom.add_surface_filling(l) for l in ll]
        # Combine the surfaces to avoid seams
        #new_surfs = [
        #        self.add_surface_filling(s[:4]),
        #        self.add_compound(s[4:])
        #        ]

        # Create the surface loop and volume
        sl = geom.add_surface_loop(s)
        v1= geom.add_volume(sl)
        
        newkey = objs.addobj(v1, 'bl')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]
        '''
vtdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
           ('axis', VtableElement('axis', type='float',
                             guilabel = 'Axis',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,1],
                             tip = "Center of Circle" )),           
          ('r1', VtableElement('r1', type='float',
                               guilabel = 'radius1',
                               default = 1.0,
                               tip = "r1")),
          ('r2', VtableElement('r2', type='float',
                               guilabel = 'radius2',
                               default = 0.0,
                               tip = "r2")),
          ('angle', VtableElement('angle', type='float',
                               guilabel = 'Angle',
                               default = 360,
                               tip = "angle")),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                                   tip = "characteristc length from point" )),)

class Cone(GeomPB):
    vt = Vtable(vtdata)
    def build_geom(self, geom, objs):
        x0,  d0,  r1, r2, angle, lcar = self.vt.make_value_or_expression(self)

        v1 = geom.add_cone(x0[0], x0[1], x0[2], d0[0], d0[1], d0[2],
                           r1, r2,  angle/180*np.pi)
        newkey = objs.addobj(v1, 'cn')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]

vtdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
           ('axis', VtableElement('axis', type='float',
                             guilabel = 'Axis',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,1],
                             tip = "Center of Circle" )),           
          ('radius', VtableElement('radius', type='float',
                               guilabel = 'Radius',
                               default = 1.0,
                               tip = "radius")),
          ('angle', VtableElement('angle', type='float',
                               guilabel = 'Angle',
                               default = 360,
                               tip = "angle")),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                                   tip = "characteristc length from point" )),)
class Cylinder(GeomPB):  
    vt = Vtable(vtdata)
    def build_geom(self, geom, objs):
        x0,  d0,  r1,  angle, lcar = self.vt.make_value_or_expression(self)

        v1 = geom.add_cylinder(x0[0], x0[1], x0[2], d0[0], d0[1], d0[2],
                           r1, angle/180*np.pi)
        newkey = objs.addobj(v1, 'cn')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]
        
vtdata =  (('corner', VtableElement('corner', type='float',
                             guilabel = 'Corner',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
           ('dxdydz', VtableElement('dxdydz', type='float',
                             guilabel = 'Size',
                             suffix =('x', 'y', 'z'),
                             default = [1,1,0.1],
                             tip = "Size of Wedge" )),           
          ('ltx', VtableElement('ltx', type='float',
                               guilabel = 'Top Extendradius1',
                               default = 1.0,
                               tip = "r1")),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                                   tip = "characteristc length from point" )),)

class Wedge(GeomPB):
    vt = Vtable(vtdata)
    def build_geom(self, geom, objs):
        x0,  d0,  ltx, lcar = self.vt.make_value_or_expression(self)
        print(x0, d0, ltx, lcar)
        v1 = geom.add_wedge(x0[0], x0[1], x0[2], d0[0], d0[1], d0[2], ltx)
        newkey = objs.addobj(v1, 'wdg')
        self._objkeys = objs.keys()
        self._newobjs = [newkey]

vtdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('r1', VtableElement('r1', type='float',
                               guilabel = 'radius1',
                               default = 1.0,
                               tip = "r1")),
          ('r2', VtableElement('r2', type='float',
                               guilabel = 'radius2',
                               default = 0.0,
                               tip = "r2")),
          ('angle', VtableElement('angle', type='float',
                               guilabel = 'Angle',
                               default = 360,
                               tip = "angle")),
          ('lcar', VtableElement('lcar', type='float',
                                   guilabel = 'lcar',
                                   default = 0.0, 
                                   tip = "characteristc length from point" )),)

class Torus(GeomPB):
    vt = Vtable(vtdata)
    def build_geom(self, geom, objs):
        x0,  r1,  r2, angle, lcar = self.vt.make_value_or_expression(self)
        v1 = geom.add_torus(x0[0], x0[1], x0[2], r1, r2, angle)
        newkey = objs.addobj(v1, 'trs')
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
                                   default = 0.0, 
                              tip = "characteristc length from point" )),)
class Line(GeomPB):
    vt = Vtable(pdata)
    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        self.vt.attribute_set(v)
        v["make_spline"] = False
        return v
    
    def panel1_param(self):
        ll = GeomPB.panel1_param(self)
        ll.append(["Spline",
                    self.make_spline,  3, {"text":""}])
        return ll
        
    def get_panel1_value(self):
        v = GeomPB.get_panel1_value(self)
        return v + [self.make_spline]

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        GeomPB.import_panel1_value(self, v[:-1])        
        self.make_spline = v[-1]

    def panel1_tip(self):
        tip = GeomPB.panel1_tip(self)
        return tip + ['make spline curve']

    def build_geom(self, geom, objs):
        xarr, yarr, zarr,  lcar = self.vt.make_value_or_expression(self)
        if len(xarr) < 2: return
        try:
           pos = np.vstack((xarr, yarr, zarr)).transpose()
        except:
           print("can not make proper input array")
           return

        pts = []
        for ii, p in enumerate(pos):
            pt = geom.add_point(p, lcar,
                                mask = (ii == 0 or ii == len(pos)-1))
            pts.append(pt)

        if not self.make_spline:
            pts1 = pts[:-1]
            pts2 = pts[1:]

            newkeys = []
            for p1, p2 in zip(pts1, pts2):
                ln = geom.add_line(p1, p2)
                newkeys.append(objs.addobj(ln, 'ln'))
            self._newobjs = newkeys
        else:     
            spline = geom.add_spline(pts)
            newobj = objs.addobj(spline, 'sp')
            self._newobjs = [newobj]
        newobj1 = objs.addobj(pts[0], 'pt')
        newobj2 = objs.addobj(pts[-1], 'pt')
        self._newobjs.append(newobj1)
        self._newobjs.append(newobj2)
        
        self._objkeys = objs.keys()

class Polygon(GeomPB):
    vt = Vtable(pdata)
    def build_geom(self, geom, objs):
        xarr, yarr, zarr,  lcar = self.vt.make_value_or_expression(self)
        if len(xarr) < 2: return
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

'''
object created from given entities
'''
ldata =  (('points', VtableElement('pts', type='string',
                                    guilabel = 'Points',
                                    default = "",
                                    tip = "points to be connected")), )
class Spline(GeomPB):
    vt = Vtable(ldata)

    def build_geom(self, geom, objs):
        pts = self.vt.make_value_or_expression(self)
        pts = [x.strip() for x in pts[0].split(',')]
        
        pts = [objs[x] for x in pts]
        
        spline = geom.add_spline(pts)
        newobj = objs.addobj(spline, 'sp')
        
        self._newobjs = [newobj]
        self._objkeys = objs.keys()
                           
class CreateLine(GeomPB):
    vt = Vtable(ldata)
    def build_geom(self, geom, objs):
        pts = self.vt.make_value_or_expression(self)
        pts = [x.strip() for x in pts[0].split(',')]        

        pts0 = pts[:-1]
        pts1 = pts[1:]
        self._newobjs = []
        
        for p0, p1 in zip(pts0, pts1):
             if not p0 in objs:
                 assert False, p0 + " does not exist"
             if not p1 in objs:
                 assert False, p1 + " does not exist"
             line = geom.add_line(objs[p0], objs[p1])
             self._newobjs.append(objs.addobj(line, 'ln'))

        self._objkeys = objs.keys()
        
ldata =  (('lines', VtableElement('lines', type='string',
                                    guilabel = 'Lines',
                                    default = "",
                                    tip = "lines to be connected")), )

class LineLoop(GeomPB):
    vt = Vtable(ldata)
    def build_geom(self, geom, objs):
        pts = self.vt.make_value_or_expression(self)
        pts = [x.strip() for x in pts[0].split(',')]
        
        pts = [(objs[x] if not x.startswith('-') else objs[x[1:]]) for x in pts]
        
        spline = geom.add_line_loop(pts)
        newobj = objs.addobj(spline, 'll')
        
        self._newobjs = [newobj]
        self._objkeys = objs.keys()

                           
class CreateSurface(GeomPB):        
    vt = Vtable(ldata)
    def build_geom(self, geom, objs):
        pts = self.vt.make_value_or_expression(self)
        pts = [x.strip() for x in pts[0].split(',')]
        
        pts = [(objs[x] if not x.startswith('-') else objs[x[1:]]) for x in pts]
        
        ll = geom.add_line_loop(pts)
        newobj1 = objs.addobj(ll, 'll')
        surface = geom.add_plane_surface(ll)
        newobj2 = objs.addobj(surface, 'ps')
        
        self._newobjs = [newobj1, newobj2]
        self._objkeys = objs.keys()
        
ldata =  (('surfs', VtableElement('surfs', type='string',
                                    guilabel = 'Surfaces',
                                    default = "",
                                    tip = "surfacess to be connected")), )       
class SurfaceLoop(GeomPB):
    vt = Vtable(ldata)
    def build_geom(self, geom, objs):
        pts = self.vt.make_value_or_expression(self)
        pts = [x.strip() for x in pts[0].split(',')]
        pts = [(objs[x] if not x.startswith('-') else objs[x[1:]]) for x in pts]
        
        sl = geom.add_surface_loop(pts)
        newobj = objs.addobj(sl, 'sl')
        
        self._newobjs = [newobj]
        self._objkeys = objs.keys()
                           
class CreateVolume(GeomPB):        
    vt = Vtable(ldata)
    def build_geom(self, geom, objs):
        pts = self.vt.make_value_or_expression(self)
        pts = [x.strip() for x in pts[0].split(',')]
        pts = [(objs[x] if not x.startswith('-') else objs[x[1:]]) for x in pts]
        
        sl = geom.add_surface_loop(pts)
        newobj1 = objs.addobj(sl, 'sl')
        vol = geom.add_volume(sl)
        newobj2 = objs.addobj(vol, 'vol')
        
        self._newobjs = [newobj1, newobj2]
        self._objkeys = objs.keys()

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


'''
 objection transformations
'''
data0 =  (('target_object', VtableElement('target_object', type='string',
                                          guilabel = 'Object',
                                          default = "",
                                          tip = "object to move")), 
          ('dx', VtableElement('dx', type='float',
                                   guilabel = 'dx',
                                   default = '0.0',
                                   tip = "x-displacement")),
          ('dy', VtableElement('dy', type='float',
                                   guilabel = 'dy',
                                   default = '0.0',                               
                                   tip = "x-displacement")),
          ('dz', VtableElement('dz', type='float',
                                   guilabel = 'dz',
                                   default = '0.0',
                                   tip = "z-displacement")),
          ('keep_org', VtableElement('kepp_org', type='bool',
                                      guilabel = 'Copy',
                                      default = True,
                                      tip = "Keep original")), )


class Move(GeomPB): # tanslate in gmsh
    vt = Vtable(data0)
    def build_geom(self, geom, objs):          
        targets, dx, dy, dz, keep  = self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]
          
        newkeys = []
        tt = [objs[t] for t in targets]
        if keep:
           tt = geom.copy(tt)          
        geom.translate(tt, dx, dy, dz)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'mv'))          
        self._objkeys = objs.keys()
        self._newobjs = newkeys

data0 =  (('target_object', VtableElement('target_object', type='string',
                                          guilabel = 'Object',
                                          default = "",
                                          tip = "object to move")), 
          ('cx', VtableElement('cx', type='float',
                                   guilabel = 'x(center)',
                                   default = '0.0',
                                   tip = "point on revolustion axis")),            
          ('cy', VtableElement('cy', type='float',
                                   guilabel = 'y(center)',
                                   default = '0.0',
                                   tip = "point on revolustion axis")),           
          ('cz', VtableElement('cz', type='float',
                                   guilabel = 'z(center)',
                                   default = '0.0',
                                   tip = "point on revolustion axis")),
          ('ax', VtableElement('ax', type='float',
                                   guilabel = 'ax',
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('ay', VtableElement('ay', type='float',
                                   guilabel = 'ay',
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('az', VtableElement('az', type='float',
                                   guilabel = 'az',
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('angle', VtableElement('angle', type='float',
                                   guilabel = 'angle',
                                   default = '0.0',
                                   tip = "angle of revoluiton")),
          ('keep_org', VtableElement('kepp_org', type='bool',
                                      guilabel = 'Copy',
                                      default = True,
                                      tip = "Keep original")), )
        
class Rotate(GeomPB):
    vt = Vtable(data0)
    def build_geom(self, geom, objs):          
        targets, cx, cy, cz, ax, ay, az, angle, keep  = self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]
          
        newkeys = []
        tt = [objs[t] for t in targets]        
        if keep:
           tt = geom.copy(tt)          
        geom.rotate(tt, cx, cy, cz, ax, ay, az, np.pi*angle/180.)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'rot'))          
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
data0 =  (('target_object', VtableElement('target_object', type='string',
                                          guilabel = 'Object',
                                          default = "",
                                          tip = "object to move")), 
          ('cx', VtableElement('cx', type='float',
                                   guilabel = 'x(center)',
                                   default = '0.0',
                                   tip = "point on revolustion axis")),            
          ('cy', VtableElement('cy', type='float',
                                   guilabel = 'y(center)',
                                   default = '0.0',
                                   tip = "point on revolustion axis")),           
          ('cz', VtableElement('cz', type='float',
                                   guilabel = 'z(center)',
                                   default = '0.0',
                                   tip = "point on revolustion axis")),
          ('scalex', VtableElement('scalex', type='float',
                                   guilabel = 'X scale',
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('scaley', VtableElement('scaley', type='float',
                                   guilabel = 'Y scale',
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('scalez', VtableElement('scalez', type='float',
                                   guilabel = 'Z scale',
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('keep_org', VtableElement('kepp_org', type='bool',
                                      guilabel = 'Copy',
                                      default = True,
                                      tip = "Keep original")), )

  
class Scale(GeomPB):  # Dilate in gmsh
    vt = Vtable(data0)
    def build_geom(self, geom, objs):          
        targets, cx, cy, cz, sx, sy, sz, keep  = self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]
          
        newkeys = []
        tt = [objs[t] for t in targets]
        if keep:
           tt = geom.copy(tt)          
        geom.dilate(tt, cx, cy, cz, sx, sy, sz)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'sc'))          

        self._objkeys = objs.keys()
        self._newobjs = newkeys

data0 =  (('target_object', VtableElement('target_object', type='string',
                                          guilabel = 'Object',
                                          default = "",
                                          tip = "object to move")), 
          ('flip_ax', VtableElement('flip_ax', type='float',
                                    guilabel = 'Flip Axis X',
                                    default = '0.0',
                                    tip = "direction on flip axis")),            
          ('flip_ay', VtableElement('flip_ay', type='float',
                                    guilabel = 'Flip Axis Y',
                                    default = '0.0',
                                    tip = "direction on flip axis")),           
          ('flip_az', VtableElement('flip_az', type='float',
                                    guilabel = 'Flip Axis Z',
                                   default = '0.0',
                                    tip = "direction on flip axis")), 
          ('flip_d', VtableElement('flip_d', type='float',
                                   guilabel = 'Offset', 
                                   default = '0.0',
                                   tip = "direction of revolustion axis")),
          ('keep_org', VtableElement('kepp_org', type='bool',
                                      guilabel = 'Copy',
                                      default = True,
                                      tip = "Keep original")), )

class Flip(GeomPB):
    vt = Vtable(data0)
    def build_geom(self, geom, objs):          
        targets, a, b, c, d,  keep  = self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]
          
        newkeys = []
        tt = [objs[t] for t in targets]
        if keep:
           tt = geom.copy(tt)          
        geom.symmetrize(tt, a, b, c, d)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'flp'))          

        self._objkeys = objs.keys()
        self._newobjs = newkeys
  
data0 =  (('target_object', VtableElement('target_object', type='string',
                                          guilabel = 'Object',
                                          default = "",
                                          tip = "object to move")), )

class Copy(GeomPB):
    vt = Vtable(data0)  
    def build_geom(self, geom, objs):
        targets  = self.vt.make_value_or_expression(self)[0]
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        tt = [objs[t] for t in targets]
        ret = geom.copy(tt)
        for r in ret:
            newkeys.append(objs.addobj(r, 'cp'))
        
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
data0 =  (('target_object', VtableElement('target_object', type='string',
                                          guilabel = 'Object',
                                          default = "",
                                          tip = "object to move")), 
          ('recursive', VtableElement('recursive', type='bool',
                                      guilabel = 'Recursive',
                                      default = True,
                                      tip = "delete recursively")), )

class Remove(GeomPB):
    vt = Vtable(data0)  
    def build_geom(self, geom, objs):
        targets, recursive= self.vt.make_value_or_expression(self)
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        tt = [objs[t] for t in targets]
        geom.remove(tt, recursive=recursive)
        for t in targets:
             del objs[t]
                               
        self._objkeys = objs.keys()
        self._newobjs = newkeys

class GeomPB_Bool(GeomPB):
    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        self.vt.attribute_set(v)
        v["delete_input"] = True
        v["delete_tool"] = True        
        return v
    
    def panel1_param(self):
        ll = GeomPB.panel1_param(self)
        ll.append(["Delete Input",
                    self.delete_input,  3, {"text":""}])
        ll.append(["Delete Tool",
                    self.delete_tool,  3, {"text":""}])
        return ll
        
    def get_panel1_value(self):
        v = GeomPB.get_panel1_value(self)
        return v + [self.delete_input, self.delete_tool]

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        GeomPB.import_panel1_value(self, v[:-2])        
        self.delete_input = v[-2]
        self.delete_tool = v[-1]

    def panel1_tip(self):
        tip = GeomPB.panel1_tip(self)
        return tip + ['delete input objects'] + ['delete tool objects']


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
                          removeObject = self.delete_input,
                          removeTool = self.delete_tool)

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
        tp = [x.strip() for x in tp[0].split(',')]
        if len(tp) < 2: return

        input_entity = [objs[x] for x in tp[:1]]
        tool_entity  = [objs[x] for x in tp[1:]]
        ret = geom.boolean_union(
                          input_entity,
                          tool_entity,
                          removeObject = self.delete_input,
                          removeTool = self.delete_tool)
        
        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'uni'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o,  get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
class Intersection(GeomPB_Bool):    
    vt = Vtable(udata)
    def build_geom(self, geom, objs):
        tp  = self.vt.make_value_or_expression(self)
        tp = [x.strip() for x in tp[0].split(',')]
        if len(tp) < 2: return

        input_entity = [objs[x] for x in tp[:1]]
        tool_entity  = [objs[x] for x in tp[1:]]
        ret = geom.boolean_intersection(
                          input_entity,
                          tool_entity,
                          removeObject = self.delete_input,
                          removeTool = self.delete_tool)
        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'its'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o,  get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys

class Fragments(GeomPB_Bool):    
    vt = Vtable(ddata)
    def build_geom(self, geom, objs):
        tp  = self.vt.make_value_or_expression(self)
        tp = [x.strip() for x in tp[0].split(',')]
        if len(tp) < 2: return

        input_entity = [objs[x] for x in tp[:1]]
        tool_entity  = [objs[x] for x in tp[1:]]
        ret = geom.boolean_fragments(
                          input_entity,
                          tool_entity,
                          removeObject = self.delete_input,
                          removeTool = self.delete_tool)

        newkeys = []
        newkeys.append(objs.addobj(ret[0], 'frag'))
        if len(ret) > 1:
            for o in ret[1:]:
                newkeys.append(objs.addobj(o, get_geom_key(o)))
            
        self._objkeys = objs.keys()
        self._newobjs = newkeys
        
          

        
    
