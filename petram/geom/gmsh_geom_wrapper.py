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

def id2dimtag(en):
    if isinstance(en, LineID):
        return (1, int(en))
    elif isinstance(en, SurfaceID):
        return (2, int(en))
    elif hasattr(en, 'surface'):
        return (2, int(en))
    elif isinstance(en, VolumeID):
        return (3, int(en))
    else:
        assert False, "Illegal entity"

def dimtag2id(dimtags):        
    out3 = []; out2 = []; out1 = []
    for dim, tag in dimtags:
        if dim == 3 and not tag in out3:
           out3.append(VolumeID(tag))                               
        elif dim == 2 and not tag in out2:
           out2.append(SurfaceID(tag))                               
        elif dim == 1 and not tag in out1:
           out1.append(LineID(tag))
    return out3 + out2 + out1                     
   
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

    @property
    def dim(self):
        if len(self.model.getEntities(3)) > 0: return 3
        if len(self.model.getEntities(2)) > 0: return 2
        if len(self.model.getEntities(1)) > 0: return 1
        return 0
     
    def add_point(self, p, lcar=0.0, mask=True):
        p = tuple(p)
        if not p in self._point_loc:
            self.p = self.p + 1
            self.factory.addPoint(p[0], p[1], p[2], lcar, self.p)            
            self._point_loc[p] = self.p
            
        p_id = self._point_loc[p]
        if mask : self._point_mask.append(p_id)
        return p_id
        
    def add_line(self, p1, p2):
        self.l = self.l + 1      
        self.factory.addLine(p1, p2, self.l)
        
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
        tags = list(np.atleast_1d(pts))       
        self.ll = self.ll+1
        self.factory.addCurveLoop(tags, self.ll)
        return self.ll
        
    def add_surface_loop(self, sl):
        tags = list(np.atleast_1d(sl))              
        self.sl = self.sl+1
        self.factory.addSurfaceLoop(tags, self.sl)        
        return self.sl
      
    def add_volume(self, shells):
        tags = list(np.atleast_1d(shells))              
        self.v = self.v+1
        self.factory.addVolume(tags, self.v)
        return self.sl
      

    def add_polygon(self, pos, lcar = 0.0):
        pts = [self.add_point(p, lcar=lcar) for p in pos]
        lns = [self.add_line(pts[i], pts[i+1]) for i in range(len(pts)-1)]
        lns.append(self.add_line(pts[-1], pts[0]))
        ll = self.add_line_loop(lns)
        sl = self.add_plane_surface((ll,))
        ret =  Polygon(sl, ll, lcar)
        return ret

    def _boolean_xxx(self, m, input_entity, tool_entity,
                     removeObject=False, removeTool=False, delete=False):
       
        def get_dimtag(entity):
           dimtags = []
           for en in entity:
               dimtags.append(id2dimtag(en))
           return dimtags

        dimtag1 = get_dimtag(input_entity)
        dimtag2 = get_dimtag(tool_entity)
                               
        if delete:
             removeObject=True
             removeTool=True

        m = getattr(self.factory, m)
        print("remove", removeObject, removeTool)
        
        dimtag3, dimtagMap = m(dimtag1, dimtag2,
                               removeObject=removeObject,
                               removeTool=removeTool)

        return dimtag2id(dimtag3)                
        
    def boolean_intersection(self, input_entity, tool_entity,
                             removeObject=False, removeTool=False, delete=False):
        return self._boolean_xxx('intersect', input_entity, tool_entity,
                                 removeObject=removeObject, removeTool=removeTool,
                                 delete=delete)
                               
    def boolean_union(self, input_entity, tool_entity,
                      removeObject=False, removeTool=False, delete=False):
        return self._boolean_xxx('fuse', input_entity, tool_entity,
                                 removeObject=removeObject, removeTool=removeTool,
                                 delete=delete)

     
    def boolean_difference(self, input_entity, tool_entity,
                           removeObject=False, removeTool=False, delete=False):
        return self._boolean_xxx('cut', input_entity, tool_entity,
                                 removeObject=removeObject, removeTool=removeTool,
                                 delete=delete)


    def boolean_fragments(self, input_entity, tool_entity,
                          removeObject=False, removeTool=False, delete=False):
        return self._boolean_xxx('fragment', input_entity, tool_entity,
                                 removeObject=removeObject, removeTool=removeTool,
                                 delete=delete)
        
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

    def remove(self, entity, recursive=False):
        dimtags = []
        for en in entity:
            dimtags.append(id2dimtag(en))
        self.factory.remove(dimtags, recursive=recursive)
        return []
     
    def copy(self, entity):
        dimtags = []
        for en in entity:
            dimtags.append(id2dimtag(en))
        dimtags2 = self.factory.copy(dimtags)
        return dimtag2id(dimtags2)                                
        
    def rotate(self, entity, x, y, z, ax, ay, az, angle):
        dimtags = []
        for en in entity:
            dimtags.append(id2dimtag(en))
        dimtags2 = self.factory.rotate(dimtags, x, y, z, ax, ay, az, angle)
        return dimtag2id(dimtags2)                        

    def translate(self, entity, dx, dy, dz):
        dimtags = []
        for en in entity:
            dimtags.append(id2dimtag(en))
        self.factory.translate(dimtags, dx, dy, dz)
        return []

    def dilate(self, entity, x, y, z, a, b, c):
        dimtags = []
        for en in entity:
            dimtags.append(id2dimtag(en))
        self.factory.dilate(dimtags, x, y, z, a, b, c)
        return []
     
    def symmetrize(self, entity, a, b, c, d):
        dimtags = []
        for en in entity:
            dimtags.append(id2dimtag(en))
        self.factory.symmetrize(dimtags, a, b, c, d)
        return []
       

'''
   def addCircleArc(startTag, centerTag, endTag, tag=-1, nx=0., ny=0., nz=0.):        
   def addEllipse(x, y, z, r1, r2, tag=-1, angle1=0., angle2=2*pi):
   def addEllipseArc(startTag, centerTag, majorTag, endTag, tag=-1, nx=0., ny=0., nz=0.)
   def addBezier(pointTags, tag=-1)
   def addCurveLoop(curveTags, tag=-1)
   def addSurfaceFilling(wireTags, tag=-1, sphereCenterTag=-1)
   def twist(dimTags, x, y, z, dx, dy, dz, ax, ay, az, angle, numElements=[], heights=[], recombine=False)
   def addDisk(xc, yc, zc, rx, ry, tag=-1)
   def addSphere(xc, yc, zc, radius, tag=-1, angle1=-pi/2, angle2=pi/2, angle3=2*pi)
   def addBox(x, y, z, dx, dy, dz, tag=-1)
   def addCylinder(x, y, z, dx, dy, dz, r, tag=-1, angle=2*pi)
   def addCone(x, y, z, dx, dy, dz, r1, r2, tag=-1, angle=2*pi)
   def addWedge(x, y, z, dx, dy, dz, tag=-1, ltx=0.)
   def addTorus(x, y, z, r1, r2, tag=-1, angle=2*pi)

   def addThruSections(wireTags, tag=-1, makeSolid=True, makeRuled=False)
   def addThickSolid(volumeTag, excludeSurfaceTags, offset, tag=-1)
