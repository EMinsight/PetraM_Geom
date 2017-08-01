import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('Primitives')

from petram.phys.vtable import VtableElement, Vtable
from petram.geom.geom_model import Geom


cdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('radius', VtableElement('radius', type='float',
                                   guilabel = 'r',
                                   default = 1.0, 
                                   tip = "radius" )),
          ('normal', VtableElement('normal', type='float',
                                   guilabel = 'n',
                                   suffix =('x', 'y', 'z'),
                                   default = [0, 0, 1], 
                                   tip = "normal vector" )),)
class Circle(Geom):
    vt = Vtable(cdata)

rdata =  (('center', VtableElement('center', type='float',
                             guilabel = 'Center',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('elen1', VtableElement('elen1', type='float',
                                   guilabel = 'A',
                                   default = 1.0, 
                                   tip = "length of edge (1)" )),
          ('elen2', VtableElement('elen2', type='float',
                                   guilabel = 'B',
                                   default = 1.0, 
                                   tip = "length of edge (2)" )),
          ('normal', VtableElement('normal', type='float',
                                   guilabel = 'n',
                                   suffix =('x', 'y', 'z'),
                                   default = [0, 0, 1], 
                                   tip = "normal vector" )),)

class Rect(Geom):
    vt = Vtable(rdata)
    
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
                              tip = "Z" )),)

class Polygon(Geom):
    vt = Vtable(pdata)

edata =  (('paxis', VtableElement('paxis', type='float',
                             guilabel = 'Point on Axis',
                             suffix =('x', 'y', 'z'),
                             default = [0,0,0],
                             tip = "Center of Circle" )),
          ('taxis', VtableElement('taxis', type='float',
                                   guilabel = 'Translation Axis',
                                   suffix =('x', 'y', 'z'),
                                   default = [0, 0, 1],
                                   tip = "translation axis" )),
          ('raxis', VtableElement('elen2', type='float',
                                   guilabel = 'Rotation Axis',
                                   suffix =('x', 'y', 'z'),
                                   default = [0, 0, 1],
                                   tip = "rotation axis" )),
          ('angle', VtableElement('angle', type='float',
                                   guilabel = 'Angle',
                                   default = 0.0,
                                   tip = "angle")),)

class Extrude(Geom):    
    vt = Vtable(edata)
    
    
