
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
   def __neg__(self):
       return LineID(-int(self))
    
class SurfaceID(GeomIDBase):
   def __add__(self, v):
       return SurfaceID(int(self) + v)
   def __neg__(self):
       return SurfaceID(-int(self))
   
class VolumeID(GeomIDBase):   
   def __add__(self, v):
       return VolumeID(int(self) + v)
   def __neg__(self):
       return VolumeID(-int(self))
    
class LineLoopID(GeomIDBase):   
   def __add__(self, v):
       return LineLoopID(int(self) + v)
   def __neg__(self):
       return LineLoopID(-int(self))
    
class SurfaceLoopID(GeomIDBase):   
   def __add__(self, v):
       return SurfaceLoopID(int(self) + v)
   def __neg__(self):
       return SurfaceLoopID(-int(self))

def id2dimtag(en):
    if isinstance(en, VertexID):
        return (0, int(en))
    elif isinstance(en, LineID):
        return (1, int(en))
    elif isinstance(en, SurfaceID):
        return (2, int(en))
    elif hasattr(en, 'surface'):
        return (2, int(en))
    elif isinstance(en, VolumeID):
        return (3, int(en))
    else:
        assert False, "Illegal entity"+str(en)

def get_dimtag(entity):
    dimtags = []
    for en in entity:
        dimtags.append(id2dimtag(en))
    return dimtags

def dimtag2id(dimtags):        
    out3 = []; out2 = []; out1 = []; out0 = []
    for dim, tag in dimtags:
        if dim == 3 and not tag in out3:
           out3.append(VolumeID(tag))                               
        elif dim == 2 and not tag in out2:
           out2.append(SurfaceID(tag))                               
        elif dim == 1 and not tag in out1:
           out1.append(LineID(tag))
        elif dim == 0 and not tag in out1:
           out0.append(VertexID(tag))
    return out3 + out2 + out1 + out0

class Geometry(object):
    def __init__(self, *args, **kwargs):
        self._point_loc = {}
        gmsh.initialize()

        gmsh.option.setNumber("General.Terminal", 1)
        
        modelname = kwargs.pop("modelname", "model1")
        gmsh.model.add(modelname)

        self.model = gmsh.model
        self.factory = gmsh.model.occ

        #self.p = VertexID(0) 
        #self.l = LineID(0)   
        #self.ll = LineLoopID(0)  
        #self.s = SurfaceID(0)   
        #self.sl = SurfaceLoopID(0)
        #self.v = VolumeID(0)      

        self._point = {}
        self._point_mask = []
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        
    def set_factory(self, factory_type):
        pass

    def clear(self):
        gmsh.clear()
        
    def getBoundingBox(self):
        xmax = -np.inf
        xmin =  np.inf
        ymax = -np.inf
        ymin =  np.inf
        zmax = -np.inf
        zmin =  np.inf
        
        def update_maxmin(dim, tag, xmin, ymin, zmin, xmax, ymax, zmax):
            x1, y1, z1, x2, y2, z2 = self.model.getBoundingBox(dim, tag)           
            xmax = np.max([xmax, x2])
            ymax = np.max([ymax, y2])
            zmax = np.max([zmax, z2])
            xmin = np.min([xmin, x1])
            ymin = np.min([ymin, y1])
            zmin = np.min([zmin, z1])
            return xmin, ymin, zmin, xmax, ymax, zmax
         
        if len(self.model.getEntities(3)) != 0:
           for dim, tag in self.model.getEntities(3):
              xmin, ymin, zmin, xmax, ymax, zmax = update_maxmin(dim, tag,
                                                                 xmin, ymin, zmin,
                                                                 xmax, ymax, zmax)
        elif len(self.model.getEntities(2)) != 0:
           for dim, tag in self.model.getEntities(2):
              xmin, ymin, zmin, xmax, ymax, zmax = update_maxmin(dim, tag,
                                                                 xmin, ymin, zmin,
                                                                 xmax, ymax, zmax)
        elif len(self.model.getEntities(1)) != 0:
           for dim, tag in self.model.getEntities(1):
              xmin, ymin, zmin, xmax, ymax, zmax = update_maxmin(dim, tag,
                                                                 xmin, ymin, zmin,
                                                                 xmax, ymax, zmax)
        else:
           for dim, tag in self.model.getEntities(0):
              xmin, ymin, zmin, xmax, ymax, zmax = update_maxmin(dim, tag,
                                                                 xmin, ymin, zmin,
                                                                 xmax, ymax, zmax)

        return xmin, ymin, zmin, xmax, ymax, zmax
     
    def getObjSizes(self):
        size = []
        for dim, tag in self.model.getEntities():
            x1, y1, z1, x2, y2, z2 = self.model.getBoundingBox(dim, tag)
            s = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
            size.append((dim, tag, s))
        return size
     
    def getVertexCL(self):
        from collections import defaultdict
        
        lcar = defaultdict(lambda: np.inf)
        
        for dim, tag in self.model.getEntities(1):
            x1, y1, z1, x2, y2, z2 = self.model.getBoundingBox(dim, tag)
            s = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
            bdimtags = self.model.getBoundary(((dim, tag,),), oriented=False)
            for bdim, btag in bdimtags:
                lcar[btag] = min((lcar[btag], s))
        return dict(lcar)
       
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
        #if not p in self._point_loc:
        pp = self.factory.addPoint(p[0], p[1], p[2], lcar)
        self._point_loc[p] = VertexID(pp)
        #print("made point ", pp, p)
            
        p_id = self._point_loc[p]
        self._point[p_id]=np.array(p)
        if mask : self._point_mask.append(p_id)
        return p_id
        
    def add_line(self, p1, p2):
        l = self.factory.addLine(p1, p2)
        return LineID(l)
      
    def add_circle_arc(self, p2, pc, p3):
        l = self.factory.addCircleArc(p2, pc, p3)
        return LineID(l)        

    def add_spline(self, pts, remove_control=True):
        l = self.factory.addSpline(pts)
        if remove_control:
           dimtags = [(0, x) for x in pts[1:-1]]
           self.factory.remove(dimtags)
        return LineID(l)

    def add_plane_surface(self, tags):
        # tags : 1st element exterier, others makes hole
        tags = list(np.atleast_1d(tags))
        #self.factory.synchronize()                        
        s = self.factory.addPlaneSurface(tags)
        #self.factory.synchronize()                                
        return SurfaceID(s)

    def add_surface_filling(self, ll):
        s = self.factory.addSurfaceFilling(ll)
        return SurfaceID(s)
       
    def add_line_loop(self, pts, sign=None):
        tags = list(np.atleast_1d(pts))
        if sign is not None:
           for k, v in enumerate(sign):
               if not v: tags[k] = -tags[k]
              
        print("line loop", tags)               

        #self.factory.synchronize()                                
        #en1 = self.model.getEntities(1)
        
        ll = self.factory.addWire(tags, checkClosed=True)

        #self.factory.synchronize()                                
        #en2 = self.model.getEntities(1)
        #if len(en1) != len(en2):
        #  print("removing", tags[-1])
        #  self.factory.remove(((1, abs(tags[-1])),))

        # (note)
        #   somehow, addWire create a duplicated line sometimes
        #   here I delete input lines to enforce re-numbering.
        #
        dimtags = [(1, x) for x in tags]
        self.factory.remove(dimtags)
           
        #self.factory.synchronize()
        #en3 = self.model.getEntities(1)
        #print(en1, en2, en3)
        return LineLoopID(ll)
        
    def add_surface_loop(self, sl):
        tags = list(np.atleast_1d(sl))

        sl = self.factory.addSurfaceLoop(tags)
        return SurfaceLoopID(sl)

    def add_sphere(self, x, y, z, radius):
        v = self.factory.addSphere(x, y, z, radius)
        return VolumeID(v)
     
    def add_cone(self, x, y, z, dx, dy, dz, r1, r2, angle):
        v = self.factory.addCone(x, y, z, dx, dy, dz, r1, r2, angle=angle)
        return VolumeID(v)
   
    def add_wedge(self, x, y, z, dx, dy, dz, ltx):
        v = self.factory.addWedge(x, y, z, dx, dy, dz, -1, ltx)
        return VolumeID(v)

    def add_cylinder(self, x, y, z, dx, dy, dz, r, angle):
        v = self.factory.addCylinder(x, y, z, dx, dy, dz, r, angle=angle)
        return VolumeID(v)
     
    def add_torus(self, x, y, z, r1, r2, angle):
        v = self.factory.addTorus(x, y, z, r1, r2, -1, angle)
        return VolumeID(v)
        
    def add_volume(self, shells):
        tags = list(np.atleast_1d(shells))              
        v = self.factory.addVolume(tags)
        return VolumeID(v)
     
    def add_ellipse_arc(self, startTag, centerTag, endTag):
        a =  self._point[startTag] - self._point[centerTag]
        b =  self._point[endTag] - self._point[centerTag]
        if np.sum(a*a) > np.sum(b*b):
            l = self.factory.addEllipseArc(startTag, centerTag, endTag)
        else:
            l = self.factory.addEllipseArc(endTag, centerTag, startTag)
        return LineID(l)
                      
    def add_polygon(self, pos, lcar = 0.0):
        pts = [self.add_point(p, lcar=lcar) for p in pos]
        lns = [self.add_line(pts[i], pts[i+1]) for i in range(len(pts)-1)]
        lns.append(self.add_line(pts[-1], pts[0]))
        ll = self.add_line_loop(lns)
        sl = self.add_plane_surface((ll,))
        ret =  Polygon(sl, ll, lcar)
        return ret

    def fillet(self, volumes, curves, radii, removeVolume=True):
        volumeTags = list(np.atleast_1d(volumes))
        curveTags = list(np.atleast_1d(curves))
        radii = list(np.atleast_1d(radii))
        print(volumeTags, curveTags, radii)
        outTags = self.factory.fillet(volumeTags,
                                curveTags,
                                radii,
                                removeVolume=removeVolume)
        return [VolumeID(v[1]) for v in outTags]        

    def chamfer(self, volumes, curves, surfaces, distances, removeVolume=True):
        volumeTags = list(np.atleast_1d(volumes))
        curveTags = list(np.atleast_1d(curves))
        surfaceTags = list(np.atleast_1d(surfaces))        
        distances = list(np.atleast_1d(distances))
        outTags = self.factory.chamfer(volumeTags,
                                       curveTags,
                                       surfacesTags,
                                       distances,
                                       removeVolume=removeVolume)
        return [VolumeID(v[1]) for v in outTags]        
     
    def import_shapes(self, fileName, highestDimOnly=True, format=""):
        out_dimtags = self.factory.importShapes(fileName,
                                                highestDimOnly=highestDimOnly,
                                                format="")
        return dimtag2id(out_dimtags)
     
    def _boolean_xxx(self, m, input_entity, tool_entity,
                     removeObject=False, removeTool=False,
                     delete=False):
       
        dimtag1 = get_dimtag(input_entity)
        dimtag2 = get_dimtag(tool_entity)
                               
        if delete:
             removeObject=True
             removeTool=True

        m = getattr(self.factory, m)
        dimtag3, dimtagMap = m(dimtag1, dimtag2,
                               removeObject=removeObject,
                               removeTool=removeTool)
        self.factory.synchronize()
        self.model.getEntities()
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
     
    def boolean_union2d(self, input_entity, tool_entity,
                      removeObject=False, removeTool=False, delete=False):
       
        def get_dimtag(entity):
           dimtags = []
           for en in entity:
               dimtags.append(id2dimtag(en))
           return dimtags
       

        all_entity = input_entity + tool_entity     
        out_entity = self._boolean_xxx('fuse', input_entity, tool_entity,
                                 removeObject=removeObject, removeTool=removeTool,
                                 delete=delete)
        self.factory.synchronize()                

        xmax = -np.inf
        xmin =  np.inf
        ymax = -np.inf
        ymin =  np.inf
        zmax = -np.inf
        zmin =  np.inf
        
        def update_maxmin(dim, tag, xmin, ymin, zmin, xmax, ymax, zmax):
            x1, y1, z1, x2, y2, z2 = self.model.getBoundingBox(dim, tag)           
            xmax = np.max([xmax, x2])
            ymax = np.max([ymax, y2])
            zmax = np.max([zmax, z2])
            xmin = np.min([xmin, x1])
            ymin = np.min([ymin, y1])
            zmin = np.min([zmin, z1])
            return xmin, ymin, zmin, xmax, ymax, zmax
        
        out_dimtag = get_dimtag(out_entity)
        for dim, tag in out_dimtag:
            xmin, ymin, zmin, xmax, ymax, zmax = update_maxmin(dim, tag,
                                                               xmin, ymin, zmin,
                                                               xmax, ymax, zmax)
        dprint1("bounding box", xmin, ymin, zmin, xmax, ymax, zmax)
        
        dx = xmax-xmin
        dy = ymax-ymin
        bbx = self.factory.addRectangle(xmin-dx/10., ymin-dy/10., (zmin+zmax)/2.,
                                        dx*1.2, dy*1.2)
        out_dimtag2, dimtagMap = self.factory.cut(((2, bbx),), out_dimtag)
        #print(out_dimtag2)
        bbx = self.factory.addRectangle(xmin-dx/10., ymin-dy/10., (zmin+zmax)/2.,
                                        dx*1.2, dy*1.2)
        out_dimtag3, dimtagMap = self.factory.cut(((2,bbx),), out_dimtag2)
        self.factory.synchronize()                       
        return dimtag2id(out_dimtag3)                        

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
        self.factory.rotate(dimtags, x, y, z, ax, ay, az, angle)
        return []

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

    def extrude(self, entity, translation_axis=None,
                point_on_axis=None,
                rotation_axis=None,
                angle = 0):
       
        #for en in entity:
        dimtags = [id2dimtag(entity)]
        if translation_axis is not None:
            tax = translation_axis
            dimtags2 = self.factory.extrude(dimtags, tax[0], tax[1], tax[2],)
        else:
            pax = point_on_axis
            rax = rotation_axis
            dimtags2 = self.factory.revolve(dimtags, pax[0], pax[1], pax[2],
                                          rax[0], rax[1], rax[2], angle)
        dprint1("extrude out", dimtags2)          
        return dimtag2id(dimtags2)

    def call_synchronize(self):
        self.factory.synchronize()        
        #self.p = len(self.model.getEntities(0))

    def recursive_getBdry(self, dimtag, bdr=None):
        if bdr is None: bdr = []
        bb = self.model.getBoundary((dimtag,), oriented=False,)
        bdr.extend(bb)

        for x in bb:
           if x[0] > 0:
               bdr = self.recursive_getBdry(x, bdr=bdr)
        return bdr
       
    def get_unique_entity(self, entity):
        entity = [x for x in entity if isinstance(x, GeomIDBase)
                  and not isinstance(x, LineLoopID)
                  and not isinstance(x, SurfaceLoopID)]
        dimtags = get_dimtag(entity)
        dimtags = [x for x in sorted(dimtags)]

        self.factory.synchronize()

        allent = self.model.getEntities()
        print('all', allent)
        dimtags = [x for x in dimtags if x in allent]
        outdimtags = dimtags[:]

        print('input', dimtags)
        for dimtag in reversed(dimtags):
           bdimtags = self.recursive_getBdry(dimtag)
           print("checking", dimtag, bdimtags)
           for x in bdimtags:
               if x in outdimtags:
                   idx = outdimtags.index(x)
                   del outdimtags[idx]
               if not x in allent:
                   idx = outdimtags.index(x)
                   del outdimtags[idx]
        print('output', outdimtags)                                                     
        return dimtag2id(outdimtags)
'''
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
'''
