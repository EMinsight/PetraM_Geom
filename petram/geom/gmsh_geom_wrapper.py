from __future__ import print_function

import numpy as np

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshGeomWrapper')

from petram.phys.vtable import VtableElement, Vtable
from petram.geom.gmsh_geom_model import GmshPrimitiveBase as GeomPB
from petram.geom.gmsh_geom_model import get_geom_key



try:
   import gmsh
   has_gmsh = True   
except ImportError:
   has_gmsh = False
   assert False, "gmsh api is not found"

class Polygon(object):
    def __init__(self, s, ll, lcar):
        self.surface = SurfaceID(s)
        self.line_loop = ll
        self.lcar = lcar

class GeomIDBase(int):
   def __repr__(self):
       return self.__class__.__name__+"("+str(int(self))+")"

class VertexID(GeomIDBase):
   def __add__(self, v):
       return VertexID(int(self) + v)

class LineID(GeomIDBase):
   def __add__(self, v):
       return LineID(int(self) + v)
    
class SurfaceID(GeomIDBase):
   def __add__(self, v):
       return SurfaceID(int(self) + v)
   
class VolumeID(GeomIDBase):   
   def __add__(self, v):
       return VolumeID(int(self) + v)
    
class LineLoopID(GeomIDBase):   
   def __add__(self, v):
       return LineLoopID(int(self) + v)
    
class SurfaceLoopID(GeomIDBase):   
   def __add__(self, v):
       return SurfaceLoopID(int(self) + v)
   
class Geometry(object):
    def __init__(self, *args, **kwargs):
        self._point_loc = {}
        gmsh.initialize()

        gmsh.option.setNumber("General.Terminal", 1)
        
        modelname = kwargs.pop("modelname", "model1")
        gmsh.model.add(modelname)

        self.model = gmsh.model
        self.factory = gmsh.model.occ

        self.p = VertexID(0) 
        self.l = LineID(0)   
        self.ll = LineLoopID(0)  
        self.s = SurfaceID(0)   
        self.sl = SurfaceLoopID(0)
        self.v = VolumeID(0)      

        self.dim = -1  # track mesh goometry dimension

        self._point = {}
        self._point_mask = []
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        
    def set_factory(self, factory_type):
        pass

    @staticmethod
    def write(filename):
        gmsh.write(filename)
        
    @staticmethod        
    def finalize():        
        gmsh.finalize()
        
    def add_point(self, p, lcar=0.0, mask=True):
        p = tuple(p)
        if not p in self._point_loc:
            self.p = self.p + 1
            self.factory.addPoint(p[0], p[1], p[2], lcar, self.p)            
            self._point_loc[p] = self.p
            
        if self.dim < 1: self.dim=0
        p_id = self._point_loc[p]
        if mask : self._point_mask.append(p_id)
        return p_id
        
    def add_line(self, p1, p2):
        self.l = self.l + 1      
        self.factory.addLine(p1, p2, self.l)
        
        if self.dim < 1: self.dim=1
        return self.l
      
    def add_circle_arc(self, p2, pc, p3):
        if self.dim < 1: self.dim=1      


    def add_spline(self, pts):
        self.l = self.l + 1             
        self.factory.addSpline(pts, self.l)
        return self.l

    def add_plane_surface(self, tags):
        # tags : 1st element exterier, others makes hole
        tags = list(np.atleast_1d(tags))
        self.s = self.s+1   
        self.factory.addPlaneSurface(tags, self.s)
        if self.dim < 2: self.dim=2
        return self.s       

    def add_line_loop(self, pts):
        self.ll = self.ll+1
        self.factory.addCurveLoop(list(pts), self.ll)
        if self.dim < 1: self.dim=1
        return self.ll
        
    def add_surface_loop(self, *ca):
        self.sl = self.sl+1             
        if self.dim < 2: self.dim=2
        return self.sl        
      
    def add_volume(self, *ca):
        self.dim=3
        pass
      

    def add_polygon(self, pos, lcar = 0.0):
        pts = [self.add_point(p, lcar=lcar) for p in pos]
        lns = [self.add_line(pts[i], pts[i+1]) for i in range(len(pts)-1)]
        lns.append(self.add_line(pts[-1], pts[0]))
        ll = self.add_line_loop(lns)
        sl = self.add_plane_surface((ll,))
        
        ret =  Polygon(sl, ll, lcar)
        if self.dim < 2: self.dim=2
        return ret

    def _boolean_xxx(self, m, input_entity, tool_entity,
                     removeObject=True, removeTool=True, delete=False):
       
        def get_dimtag(entity):
           dimtags = []
           for en in entitiy:
               if isintstance(en, LineID):
                   dimtags.append((1, int(en))
               elif isintstance(en, SurfaceID):
                   dimtags.append((2, int(en))
               elif hasattr(en, 'surface'):
                   dimtags.append((2, int(en))
               elif isintstance(en, VolumeID):
                   dimtags.append((3, int(en))
               else:
                   assert False, "Illegal entity"
           return dimtags

        dimtag1 = self._boolean_get_dimtag(input_entity)
        dimtag2 = self._boolean_get_dimtag(tool_entity)
                               
        if delete:
             removeObject=True
             removeTool=True

        m = getattr(self, factory, name)                      
        dimtag3 = m(dimtag1, dimtag2, removeObject=True, removeTool=True)

        out = []
        for dim, tag in dimtag3:
            if dim == 3:
               out.append(VolumeID(tag))                               
            elif dim == 2:
               out.append(SurfaceID(tag))                               
            elif dim == 1:
               out.append(LineID(tag))
        
    def boolean_intersection(self, input_entity, tool_entity,
                             removeObject=True, removeTool=True, delete=False):
        return self._boolean_xxx('intersect', input_entity, tool_entity,
                          removeObject=True, removeTool=True, delete=False):
                               
    def boolean_union(self, input_entity, tool_entity,
                      removeObject=True, removeTool=True, delete=False):
        return self._boolean_xxx('fuse', input_entity, tool_entity,
                          removeObject=True, removeTool=True, delete=False):
     
    def boolean_difference(self, input_entity, tool_entity,
                           removeObject=True, removeTool=True, delete=False):
        return self._boolean_xxx('cut', input_entity, tool_entity,
                          removeObject=True, removeTool=True, delete=False):

    def boolean_fragments(self, input_entity, tool_entity,
                          removeObject=True, removeTool=True, delete=False):
        return self._boolean_xxx('fragment', input_entity, tool_entity,
                          removeObject=True, removeTool=True, delete=False):

        
    def apply_fragments(self):
        self.factory.synchronize()        
        if self.dim == 0: return

        dimtags =  self.model.getEntities(self.dim)
        if len(dimtags) != 1:
            self.factory.fragment(dimtags[:1], dimtags[1:],
                                  removeObject=True, removeTool=True)
            
            self.factory.synchronize()
            
        if self.dim > 1:
           dimtags =  self.model.getEntities(self.dim-1)
           if len(dimtags) != 1:           
               self.factory.fragment(dimtags[:1], dimtags[1:],
                                  removeObject=True, removeTool=True)
               self.factory.synchronize()               
        self.factory.synchronize()               
        if self.dim > 2:
           dimtags =  self.model.getEntities(self.dim-2)
           if len(dimtags) != 1:           
               self.factory.fragment(dimtags[:1], dimtags[1:],
                                  removeObject=True, removeTool=True)
               self.factory.synchronize()
               
                             
    def call_synchronize(self):
        self.factory.synchronize()        
        self.p = len(self.model.getEntities(0))     
        
            

          

        
    
