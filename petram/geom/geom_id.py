class GeomIDBase(int):
    def __init__(self, v=0):
        super(GeomIDBase, self).__init__()
        super(GeomIDBase, self).__eq__(int(v))

    def __repr__(self):
        return self.__class__.__name__+"("+str(int(self))+")"
    def to_dimtag(self):
        print((self.dim, int(self)))
        return (self.dim, int(self))

class VertexID(GeomIDBase):
    dim = 0
    idx = 0
    name = 'pt'
    def __add__(self, v):
        return VertexID(int(self) + v)

class LineID(GeomIDBase):
    dim = 1
    idx = 1
    name = 'ln'
    def __add__(self, v):
        return LineID(int(self) + v)
    def __neg__(self):
        return LineID(-int(self))

class LineLoopID(GeomIDBase):
    idx = 2
    name = 'll'
    def __add__(self, v):
        return LineLoopID(int(self) + v)
    def __neg__(self):
        return LineLoopID(-int(self))

class SurfaceID(GeomIDBase):
    dim = 2
    idx = 3
    name = 'sf'
    def __add__(self, v):
        return SurfaceID(int(self) + v)
    def __neg__(self):
        return SurfaceID(-int(self))

class SurfaceLoopID(GeomIDBase):
    idx = 4
    name = 'sl'
    def __add__(self, v):
        return SurfaceLoopID(int(self) + v)
    def __neg__(self):
        return SurfaceLoopID(-int(self))

class VolumeID(GeomIDBase):
    dim = 3
    idx = 5
    name = 'vol'
    def __add__(self, v):
        return VolumeID(int(self) + v)
    def __neg__(self):
        return VolumeID(-int(self))
