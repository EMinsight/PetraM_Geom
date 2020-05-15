class GeomIDBase(int):
   def __repr__(self):
       return self.__class__.__name__+"("+str(int(self))+")"
   def to_dimtag(self):
       print ((self.dim, int(self)))
       return (self.dim, int(self))
   
class VertexID(GeomIDBase):
   dim = 0
   def __add__(self, v):
       return VertexID(int(self) + v)

class LineID(GeomIDBase):
   dim = 1    
   def __add__(self, v):
       return LineID(int(self) + v)
   def __neg__(self):
       return LineID(-int(self))
    
class SurfaceID(GeomIDBase):
   dim = 2    
   def __add__(self, v):
       return SurfaceID(int(self) + v)
   def __neg__(self):
       return SurfaceID(-int(self))
   
class VolumeID(GeomIDBase):
   dim = 3    
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
