'''

   gmsh_primitive

   This file define GUI interface. See gmsh_geom_wrapper for actual
   geometry generation routine

'''
from petram.geom.gmsh_config import has_gmsh
from petram.geom.vtable_geom import *
from petram.geom.gmsh_geom_model import use_gmsh_api
from petram.geom.gmsh_geom_model import get_geom_key
from petram.geom.gmsh_geom_model import GmshPrimitiveBase as GeomPB
from petram.phys.vtable import VtableElement, Vtable
import numpy as np
from collections import defaultdict

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshPrimitives')


invalid_pdata = (('NotUsedValue', VtableElement('NotUsedValue', type='array',
                                                guilabel='not_implemented',
                                                default='0.0',
                                                tip="This panel is not implemented")),)
pdata = (('xarr', VtableElement('xarr', type='array',
                                guilabel='X',
                                default='0.0',
                                tip="X")),
         ('yarr', VtableElement('yarr', type='array',
                                guilabel='Y',
                                default='0.0',
                                tip="Y")),
         ('zarr', VtableElement('zarr', type='array',
                                guilabel='Z',
                                default='0.0',
                                tip="Z")),)


def get_numbers(objs, targets):
    return [objs[t] if t in objs else int(t) for t in targets]


class Point(GeomPB):
    vt = Vtable(pdata)


cdata = (('center', VtableElement('center', type='float',
                                  guilabel='Center',
                                  suffix=('x', 'y', 'z'),
                                  default=[0, 0, 0],
                                  tip="Center of Circle")),
         ('ax1', VtableElement('ax1', type='float',
                               guilabel='axis1',
                               suffix=('x', 'y', 'z'),
                               default=[1, 0, 0],
                               tip="axis 1")),
         ('ax2', VtableElement('ax2', type='float',
                               guilabel='axis2',
                               suffix=('x', 'y', 'z'),
                               default=[0, 1, 0],
                               tip="axis 2")),
         ('radius', VtableElement('radius', type='float',
                                  guilabel='r',
                                  default=1.0,
                                  tip="radius")),)


class Circle(GeomPB):
    vt = Vtable(cdata)


rdata = (('corner', VtableElement('corner', type='float',
                                  guilabel='Corner',
                                  suffix=('x', 'y', 'z'),
                                  default=[0, 0, 0],
                                  tip="Center of Circle")),
         ('edge1', VtableElement('edge1', type='float',
                                 guilabel='Edge(1)',
                                 suffix=('x', 'y', 'z'),
                                 default=[1, 0, 0],
                                 tip="Edge of rectangle")),
         ('edge2', VtableElement('edge2', type='float',
                                 guilabel='Edge(2)',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 1, 0],
                                 tip="Edge of rectangle")),)


class Rect(GeomPB):
    vt = Vtable(rdata)


rdata = (('corner', VtableElement('corner', type='float',
                                  guilabel='Corner',
                                  suffix=('x', 'y', 'z'),
                                  default=[0, 0, 0],
                                  tip="Center of Circle")),
         ('edge1', VtableElement('edge1', type='float',
                                 guilabel='Edge(1)',
                                 suffix=('x', 'y', 'z'),
                                 default=[1, 0, 0],
                                 tip="Edge of rectangle")),
         ('edge2', VtableElement('edge2', type='float',
                                 guilabel='Edge(2)',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 1, 0],
                                 tip="Edge of rectangle")),
         ('edge3', VtableElement('edge3', type='float',
                                 guilabel='Edge(3)',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 0, 1],
                                 tip="Edge of rectangle")),)


class Box(GeomPB):
    vt = Vtable(rdata)


vtdata = (('center', VtableElement('center', type='float',
                                   guilabel='Center',
                                   suffix=('x', 'y', 'z'),
                                   default=[0, 0, 0],
                                   tip="Center of Circle")),
          ('x_radius', VtableElement('x_radius', type='float',
                                     guilabel='radius (X)',
                                     default=1.0,
                                     tip="radius in X direction")),
          ('y_radius', VtableElement('y_radius', type='float',
                                     guilabel='radius (Y)',
                                     default=1.0,
                                     tip="radius in Y direction")),
          ('z_rarius', VtableElement('z_radius', type='float',
                                     guilabel='radius (Z)',
                                     default=1.0,
                                     tip="radius in Z direction")),
          ('angle1', VtableElement('angle1', type='float',
                                   guilabel='polar opening(1)',
                                   default=-90,
                                   tip="poloar opening angle (start)")),
          ('angle2', VtableElement('angle2', type='float',
                                   guilabel='polar opening(2)',
                                   default=90.,
                                   tip="polar opening angle (stop)")),
          ('angle3', VtableElement('angle3', type='float',
                                   guilabel='azimuthal opening',
                                   default=360.,
                                   tip="azimuthal opening")),)


class Ball(GeomPB):
    vt = Vtable(vtdata)


vtdata = (('center', VtableElement('center', type='float',
                                   guilabel='Center',
                                   suffix=('x', 'y', 'z'),
                                   default=[0, 0, 0],
                                   tip="Center of Circle")),
          ('axis', VtableElement('axis', type='float',
                                 guilabel='Axis',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 0, 1],
                                 tip="Center of Circle")),
          ('r1', VtableElement('r1', type='float',
                               guilabel='radius1',
                               default=1.0,
                               tip="r1")),
          ('r2', VtableElement('r2', type='float',
                               guilabel='radius2',
                               default=0.0,
                               tip="r2")),
          ('angle', VtableElement('angle', type='float',
                                  guilabel='Angle',
                                  default=360,
                                  tip="angle")),)


class Cone(GeomPB):
    vt = Vtable(vtdata)


vtdata = (('center', VtableElement('center', type='float',
                                   guilabel='Center',
                                   suffix=('x', 'y', 'z'),
                                   default=[0, 0, 0],
                                   tip="Center of Circle")),
          ('axis', VtableElement('axis', type='float',
                                 guilabel='Axis',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 0, 1],
                                 tip="Center of Circle")),
          ('radius', VtableElement('radius', type='float',
                                   guilabel='Radius',
                                   default=1.0,
                                   tip="radius")),
          ('angle', VtableElement('angle', type='float',
                                  guilabel='Angle',
                                  default=360,
                                  tip="angle")),)


class Cylinder(GeomPB):
    vt = Vtable(vtdata)


vtdata = (('corner', VtableElement('corner', type='float',
                                   guilabel='Corner',
                                   suffix=('x', 'y', 'z'),
                                   default=[0, 0, 0],
                                   tip="Center of Circle")),
          ('dxdydz', VtableElement('dxdydz', type='float',
                                   guilabel='Size',
                                   suffix=('x', 'y', 'z'),
                                   default=[1, 1, 0.1],
                                   tip="Size of Wedge")),
          ('ltx', VtableElement('ltx', type='float',
                                guilabel='Top Extend',
                                default=0.0,
                                tip="r1")),)


class Wedge(GeomPB):
    vt = Vtable(vtdata)


vtdata = (('center', VtableElement('center', type='float',
                                   guilabel='Center',
                                   suffix=('x', 'y', 'z'),
                                   default=[0, 0, 0],
                                   tip="Center of Circle")),
          ('r1', VtableElement('r1', type='float',
                               guilabel='radius1',
                               default=1.0,
                               tip="r1")),
          ('r2', VtableElement('r2', type='float',
                               guilabel='radius2',
                               default=0.0,
                               tip="r2")),
          ('angle', VtableElement('angle', type='float',
                                  guilabel='Angle',
                                  default=360,
                                  tip="angle")),
          ('keep_interior', VtableElement('keep_interior', type='bool',
                                          guilabel='keep interior surfaces',
                                          default=True,
                                          tip="Keep Intrior Surfaces ")), )


class Torus(GeomPB):
    vt = Vtable(vtdata)


pdata = (('xarr', VtableElement('xarr', type='array',
                                guilabel='X',
                                default='0.0',
                                tip="X")),
         ('yarr', VtableElement('yarr', type='array',
                                guilabel='Y',
                                default='0.0',
                                tip="Y")),
         ('zarr', VtableElement('zarr', type='array',
                                guilabel='Z',
                                default='0.0',
                                tip="Z")),)


class Line(GeomPB):
    vt = Vtable(pdata)

    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        self.vt.attribute_set(v)
        v["make_spline"] = False
        v["periodic"] = False
        return v

    def panel1_param(self):
        ll = GeomPB.panel1_param(self)
        ll.append(["Spline",
                   self.make_spline, 3, {"text": ""}])
        ll.append(["Close loop(OCC only)",
                   self.periodic, 3, {"text": ""}])
        return ll

    def get_panel1_value(self):
        v = GeomPB.get_panel1_value(self)
        return v + [self.make_spline, self.periodic]

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        GeomPB.import_panel1_value(self, v[:-2])
        self.make_spline = v[-2]
        self.periodic = v[-1]

    def panel1_tip(self):
        tip = GeomPB.panel1_tip(self)
        return tip + ['make spline curve', 'close loop']

    def _make_value_or_expression(self):
        #xarr, yarr, zarr,  lcar = self.vt.make_value_or_expression(self)
        return self.vt.make_value_or_expression(self)

    def add_geom_sequence(self, geom):
        gui_name = self.fullname()
        self.vt.preprocess_params(self)
        gui_param = list(self._make_value_or_expression()) + \
            [self.make_spline, self.periodic]
        geom_name = self.__class__.__name__
        geom.add_sequence(gui_name, gui_param, geom_name)


class Polygon(GeomPB):
    vt = Vtable(pdata)


ldata = (('points', VtableElement('pts', type='string',
                                  guilabel='Points',
                                  default="",
                                  tip="points to be connected")), )


class Spline(GeomPB):
    vt = Vtable(ldata)


class CreateLine(GeomPB):
    vt = Vtable(ldata)


ldata = (('lines', VtableElement('lines', type='string',
                                 guilabel='Lines',
                                 default="",
                                 tip="lines to be connected")), )


class LineLoop(GeomPB):
    vt = Vtable(ldata)


ldata = (('lines', VtableElement('lines', type='string',
                                 guilabel='Lines',
                                 default="",
                                 tip="lines to be connected")),
         ('isplane', VtableElement('isplane_org', type='bool',
                                   guilabel='Use surface filling',
                                   default=False,
                                   tip="Surface filling")), )


class CreateSurface(GeomPB):
    vt = Vtable(ldata)


ldata = (('surfs', VtableElement('surfs', type='string',
                                 guilabel='Surfaces',
                                 default="",
                                 tip="surfacess to be connected")), )


class SurfaceLoop(GeomPB):
    vt = Vtable(ldata)


class CreateVolume(GeomPB):
    vt = Vtable(ldata)


edata = (('ex_target', VtableElement('ex_target', type='string',
                                     guilabel='Target',
                                     default="",
                                     tip="extrusion target")),
         ('taxis', VtableElement_Direction('taxis', type='float',
                                           guilabel='Translation Axis',
                                           suffix=('x', 'y', 'z'),
                                           default=[0, 0, 1],
                                           tip="translation axis")),
         ('ex_len', VtableElement('exlen', type='float',
                                  guilabel='Length',
                                  default=1.0,
                                  tip="Extrusion length")), )


class Extrude(GeomPB):
    vt = Vtable(edata)


edata = (('ex_target', VtableElement('ex_target', type='string',
                                     guilabel='Target',
                                     default="",
                                     tip="extrusion target")),
         ('paxis', VtableElement('paxis', type='float',
                                 guilabel='Point on Axis',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 0, 0],
                                 tip="point on axis")),
         ('raxis', VtableElement('raxis', type='float',
                                 guilabel='Translation Axis',
                                 suffix=('x', 'y', 'z'),
                                 default=[0, 0, 1],
                                 tip="translation axis")),
         ('angle', VtableElement('angel', type='float',
                                 guilabel='angle',
                                 default=180.,
                                 tip="Extrusion length")), )


class Revolve(GeomPB):
    vt = Vtable(edata)


edata = (('ex_target', VtableElement('ex_target', type='string',
                                     guilabel='Target',
                                     default="",
                                     tip="extrusion target")),
         ('path', VtableElement('path', type='string',
                                guilabel='Path (Lines)',
                                default="",
                                tip="sweep path")),)


class Sweep(GeomPB):
    vt = Vtable(edata)


'''
 objection transformations
'''
data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('dx', VtableElement('dx', type='float',
                              guilabel='dx',
                              default=0.0,
                              tip="x-displacement")),
         ('dy', VtableElement('dy', type='float',
                              guilabel='dy',
                              default=0.0,
                              tip="x-displacement")),
         ('dz', VtableElement('dz', type='float',
                              guilabel='dz',
                              default=0.0,
                              tip="z-displacement")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Move(GeomPB):  # tanslate in gmsh
    vt = Vtable(data0)


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('ctr_rot', VtableElement('ctr_rot', type='float',
                                   guilabel='Center',
                                   suffix=('x', 'y', 'z'),
                                   default=[0, 0, 0],
                                   tip="point on revolustion axis")),
         ('ax_rot', VtableElement('ax_rot', type='float',
                                  guilabel='Axis',
                                  suffix=('x', 'y', 'z'),
                                  default=[0., 0., 1.0],
                                  tip="direction of revolustion axis")),
         ('angle', VtableElement('angle', type='float',
                                 guilabel='angle',
                                 default=180.0,
                                 tip="angle of revoluiton")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Rotate(GeomPB):
    vt = Vtable(data0)


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('ctr_scale', VtableElement('ctr_scale', type='float',
                                     guilabel='Center',
                                     suffix=('x', 'y', 'z'),
                                     default=[0, 0, 0],
                                     tip="center of scale")),
         ('size_scale', VtableElement('size_scale', type='float',
                                      guilabel='Scale',
                                      suffix=('x', 'y', 'z'),
                                      default=[1., 1., 1.],
                                      tip="scale size")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Scale(GeomPB):  # Dilate in gmsh
    vt = Vtable(data0)


class Array(GeomPB):
    data0 = (('target_object', VtableElement('target_object', type='string',
                                             guilabel='Object',
                                             default="",
                                             tip="object to move")),
             ('array_count', VtableElement('array_count', type='int',
                                           guilabel='Count', default=1,
                                           tip="Center of Circle")),
             ('displacement', VtableElement('displacement', type='float',
                                            guilabel='displacement',
                                            suffix=('x', 'y', 'z'),
                                            default=[1, 0, 0],
                                            tip="displacemnt")),)
    vt = Vtable(data0)


class ArrayRot(GeomPB):
    data0 = (('target_object', VtableElement('target_object', type='string',
                                             guilabel='Object',
                                             default="",
                                             tip="object to move")),
             ('array_count', VtableElement('array_count', type='int',
                                           guilabel='Count', default=1,
                                           tip="Center of Circle")),
             ('ctr_rot', VtableElement('ctr_rot', type='float',
                                       guilabel='Center',
                                       suffix=('x', 'y', 'z'),
                                       default=[0, 0, 0],
                                       tip="point on revolustion axis")),
             ('ax_rot', VtableElement('ax_rot', type='float',
                                      guilabel='Axis',
                                      suffix=('x', 'y', 'z'),
                                      default=[0., 0., 1.0],
                                      tip="direction of revolustion axis")),
             ('angle', VtableElement('angle', type='float',
                                     guilabel='angle',
                                     default=180.0,
                                     tip="angle of revoluiton")),)
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'ArrayRot'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('flip_ax', VtableElement('flip_ax', type='float',
                                   guilabel='Flip Axis X',
                                   default=0.0,
                                   tip="direction on flip axis")),
         ('flip_ay', VtableElement('flip_ay', type='float',
                                   guilabel='Flip Axis Y',
                                   default=0.0,
                                   tip="direction on flip axis")),
         ('flip_az', VtableElement('flip_az', type='float',
                                   guilabel='Flip Axis Z',
                                   default=1.0,
                                   tip="direction on flip axis")),
         ('flip_d', VtableElement('flip_d', type='float',
                                  guilabel='Offset',
                                  default=0.0,
                                  tip="direction of revolustion axis")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Flip(GeomPB):
    vt = Vtable(data0)


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Volume',
                                         default="",
                                         tip="object to add fillet")),
         ('curves', VtableElement('curves', type='string',
                                  guilabel='Curves',
                                  default="",
                                  tip="curves to add fillet")),
         ('radius', VtableElement('radisu', type='float',
                                  guilabel='Radius',
                                  default=1.0,
                                  tip="radisu")),)


class Fillet(GeomPB):
    vt = Vtable(data0)


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Volume',
                                         default="",
                                         tip="object to add chamfer")),
         ('curves', VtableElement('curves', type='string',
                                  guilabel='Curves',
                                  default="",
                                  tip="curves to add chamfer")),
         ('distance', VtableElement('distance', type='array',
                                    guilabel='Distance',
                                    default="1.0",
                                    tip="distance")),
         ('surfaces', VtableElement('surfaces', type='string',
                                    guilabel='Surfaces',
                                    default="",
                                    tip="distance")),)


class Chamfer(GeomPB):
    vt = Vtable(data0)


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")), )


class Copy(GeomPB):
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Copy object'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('recursive', VtableElement('recursive', type='bool',
                                     guilabel='Recursive',
                                     default=True,
                                     tip="delete recursively")), )


class Remove(GeomPB):
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Delete object'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),)


class Remove2(GeomPB):
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Keep object'


class GeomPB_Bool(GeomPB):
    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        self.vt.attribute_set(v)
        v["delete_input"] = True
        v["delete_tool"] = True
        v["keep_highest"] = False
        return v

    def panel1_param(self):
        ll = GeomPB.panel1_param(self)
        ll.append(["Delete Input",
                   self.delete_input, 3, {"text": ""}])
        ll.append(["Delete Tool",
                   self.delete_tool, 3, {"text": ""}])
        ll.append(["Keep highest dim only",
                   self.keep_highest, 3, {"text": ""}])
        ll.append(
            [None, "Note: Use prefix v(volume), f(face), l(line), p(point)", 2, None])
        return ll

    def get_panel1_value(self):
        v = GeomPB.get_panel1_value(self)
        return v + [self.delete_input,
                    self.delete_tool, self.keep_highest, None]

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        GeomPB.import_panel1_value(self, v[:-4])
        self.delete_input = v[-4]
        self.delete_tool = v[-3]
        self.keep_highest = v[-2]

    def panel1_tip(self):
        tip = GeomPB.panel1_tip(self)
        return tip + ['delete input objects'] + \
            ['delete tool objects'] + ['keep highest dim. objects']

    def add_geom_sequence(self, geom):
        gui_name = self.fullname()
        gui_param = (list(self.vt.make_value_or_expression(self)) +
                     [self.delete_input,
                      self.delete_tool, self.keep_highest])
        geom_name = self.__class__.__name__
        geom.add_sequence(gui_name, gui_param, geom_name)


ddata = (('objplus', VtableElement('objplus', type='string',
                                   guilabel='+',
                                   default="",
                                   tip="added objects")),
         ('objminus', VtableElement('objminus', type='string',
                                    guilabel='-',
                                    default="",
                                    tip="objects to be subtracted")),)


class Difference(GeomPB_Bool):
    vt = Vtable(ddata)


udata = (('objplus', VtableElement('obj1', type='string',
                                   guilabel='Input',
                                   default="",
                                   tip="objects")),
         ('tool_object', VtableElement('tool_object', type='string',
                                       guilabel='Tool Obj.',
                                       default="",
                                       tip="object to move")),)


class Union(GeomPB_Bool):
    vt = Vtable(udata)


class Union2D(GeomPB_Bool):
    vt = Vtable(udata)

    @classmethod
    def fancy_menu_name(self):
        return 'Union'


class Intersection(GeomPB_Bool):
    vt = Vtable(udata)


class Fragments(GeomPB_Bool):
    vt = Vtable(udata)


pdata = (('xarr', VtableElement('xarr', type='array',
                                guilabel='X',
                                default='0.0',
                                tip="X")),
         ('yarr', VtableElement('yarr', type='array',
                                guilabel='Y',
                                default='0.0',
                                tip="Y")),)


class Point2D(GeomPB):
    vt = Vtable(pdata)

    @classmethod
    def fancy_menu_name(self):
        return 'Point'


pdata = (('xarr', VtableElement('xarr', type='array',
                                guilabel='X',
                                default='0.0',
                                tip="X")),
         ('yarr', VtableElement('yarr', type='array',
                                guilabel='Y',
                                default='0.0',
                                tip="Y")),)


class Line2D(Line):
    vt = Vtable(pdata)

    def _make_value_or_expression(self):
        #xarr, yarr, zarr = self.vt.make_value_or_expression(self)
        xarr, yarr = self.vt.make_value_or_expression(self)
        zarr = [0.0 for x in yarr]
        return xarr, yarr, zarr

    @classmethod
    def fancy_menu_name(self):
        return 'Line'


cdata = (('center', VtableElement('center', type='float',
                                  guilabel='Center',
                                  suffix=('x', 'y', ),
                                  default=[0, 0, ],
                                  tip="Center of Circle")),
         ('ax1', VtableElement('ax1', type='float',
                               guilabel='axis1',
                               suffix=('x', 'y',),
                               default=[1, 0, ],
                               tip="axis 1")),
         ('ax2', VtableElement('ax2', type='float',
                               guilabel='axis2',
                               suffix=('x', 'y', ),
                               default=[0, 1, ],
                               tip="axis 2")),
         ('radius', VtableElement('radius', type='float',
                                  guilabel='r',
                                  default=1.0,
                                  tip="radius")),)


class Circle2D(GeomPB):
    vt = Vtable(cdata)

    @classmethod
    def fancy_menu_name(self):
        return 'Circle'


cdata = (('center', VtableElement('center', type='float',
                                  guilabel='Center',
                                  suffix=('x', 'y', ),
                                  default=[0, 0, ],
                                  tip="Center of Circle")),
         ('ax1', VtableElement('ax1', type='float',
                               guilabel='axis1',
                               suffix=('x', 'y',),
                               default=[1, 0, ],
                               tip="axis 1")),
         ('ax2', VtableElement('ax2', type='float',
                               guilabel='axis2',
                               suffix=('x', 'y', ),
                               default=[0, 1, ],
                               tip="axis 2")),
         ('radius', VtableElement('radius', type='float',
                                  guilabel='r',
                                  default=1.0,
                                  tip="radius")),
         ('angle1', VtableElement('angle1', type='float',
                                  guilabel='angle1',
                                  default=0.0,
                                  tip="radius")),
         ('angle2', VtableElement('angle2', type='float',
                                  guilabel='angle2',
                                  default=90.0,
                                  tip="radius")),
         ('fillarc', VtableElement('fillarc', type='bool',
                                   guilabel='fill',
                                   default=True,
                                   tip="fill arc")), )


class Arc2D(GeomPB):
    vt = Vtable(cdata)

    @classmethod
    def fancy_menu_name(self):
        return 'Arc'


rdata = (('corner', VtableElement('corner', type='float',
                                  guilabel='Corner',
                                  suffix=('x', 'y'),
                                  default=[0, 0],
                                  tip="Center of Circle")),
         ('edge1', VtableElement('edge1', type='float',
                                 guilabel='Edge(1)',
                                 suffix=('x', 'y'),
                                 default=[1, 0],
                                 tip="Edge of rectangle")),
         ('edge2', VtableElement('edge2', type='float',
                                 guilabel='Edge(2)',
                                 suffix=('x', 'y'),
                                 default=[0, 1],
                                 tip="Edge of rectangle")),)


class Rect2D(GeomPB):
    vt = Vtable(rdata)

    @classmethod
    def fancy_menu_name(self):
        return 'Rect'


pdata = (('xarr', VtableElement('xarr', type='array',
                                guilabel='X',
                                default='0.0',
                                tip="X")),
         ('yarr', VtableElement('yarr', type='array',
                                guilabel='Y',
                                default='0.0',
                                tip="Y")),)


class Polygon2D(GeomPB):
    vt = Vtable(pdata)

    @classmethod
    def fancy_menu_name(self):
        return 'Polygon'


class Spline2D(Spline):

    @classmethod
    def fancy_menu_name(self):
        return 'Spline'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('dx', VtableElement('dx', type='float',
                              guilabel='dx',
                              default=0.0,
                              tip="x-displacement")),
         ('dy', VtableElement('dy', type='float',
                              guilabel='dy',
                              default=0.0,
                              tip="y-displacement")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Move2D(GeomPB):  # tanslate in gmsh
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Move'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('ctr_rot', VtableElement('ctr_rot', type='float',
                                   guilabel='Center',
                                   suffix=('x', 'y',),
                                   default=[0, 0],
                                   tip="point on revolustion axis")),
         ('angle', VtableElement('angle', type='float',
                                 guilabel='angle',
                                 default=180.0,
                                 tip="angle of revoluiton")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Rotate2D(GeomPB):
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Rotate'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('flip_ax', VtableElement('flip_ax', type='float',
                                   guilabel='Flip Axis X',
                                   default=0.0,
                                   tip="direction on flip axis")),
         ('flip_ay', VtableElement('flip_ay', type='float',
                                   guilabel='Flip Axis Y',
                                   default=1.0,
                                   tip="direction on flip axis")),
         ('flip_d', VtableElement('flip_d', type='float',
                                  guilabel='Offset',
                                  default=0.0,
                                  tip="direction of revolustion axis")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Flip2D(GeomPB):
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Flip'


data0 = (('target_object', VtableElement('target_object', type='string',
                                         guilabel='Object',
                                         default="",
                                         tip="object to move")),
         ('ctr_scale', VtableElement('ctr_scale', type='float',
                                     guilabel='Center',
                                     suffix=('x', 'y'),
                                     default=[0., 0.],
                                     tip="center of scale")),
         ('size_scale', VtableElement('size_scale', type='float',
                                      guilabel='Scale',
                                      suffix=('x', 'y'),
                                      default=[1., 1.],
                                      tip="scale size")),
         ('keep_org', VtableElement('kepp_org', type='bool',
                                    guilabel='Copy',
                                    default=True,
                                    tip="Keep original")), )


class Scale2D(GeomPB):  # Dilate in gmsh
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Scale'


class Array2D(GeomPB):
    data0 = (('target_object', VtableElement('target_object', type='string',
                                             guilabel='Object',
                                             default="",
                                             tip="object to move")),
             ('array_count', VtableElement('array_count', type='int',
                                           guilabel='Count', default=1,
                                           tip="Center of Circle")),
             ('displacement', VtableElement('displacement', type='float',
                                            guilabel='displacement',
                                            suffix=('x', 'y',),
                                            default=[1, 0, ],
                                            tip="displacemnt")),)
    vt = Vtable(data0)

    @classmethod
    def fancy_menu_name(self):
        return 'Array'


edata = (('ex_target', VtableElement('ex_target', type='string',
                                     guilabel='Target',
                                     default="",
                                     tip="target object")),
         ('taxis', VtableElement_Plane('taxis', type='float',
                                       guilabel='cut plane',
                                       suffix=('a', 'b', 'c', 'd'),
                                       default=[1, 0, 0, 0],
                                       tip="Cutting Plane")),)


class SplitByPlane(GeomPB):
    vt = Vtable(edata)


data0 = (('center', VtableElement('center', type='float',
                                  guilabel='Center',
                                  suffix=('x', 'y', 'z'),
                                  default=[0, 0, 0],
                                  tip="Center of WP")),
         ('ax1', VtableElement('ax1', type='float',
                               guilabel='1st Axis',
                               suffix=('x', 'y', 'z'),
                               default=[0, 1, 0],
                               tip="1st Axis")),
         ('ax2', VtableElement('ax2', type='float',
                               guilabel='2nd Axis',
                               suffix=('x', 'y', 'z'),
                               default=[0, 0, 1],
                               tip="2nd Axis")),)


class WorkPlane(GeomPB):
    vt = Vtable(data0)

    def get_possible_child(self):
        return [Point2D, Line2D, Circle2D, Arc2D, Rect2D, Polygon2D, Spline2D,
                Move2D, Rotate2D, Flip2D, Scale2D, Array2D,
                Union2D, Intersection, Difference, Fragments, Copy, Remove,
                CreateLine, CreateSurface]

    def get_possible_child_menu(self):
        return [("", Point2D), ("", Line2D), ("", Circle2D), ("", Arc2D),
                ("", Rect2D), ("", Polygon2D), ("", Spline2D),
                ("", CreateLine), ("", CreateSurface),
                ("", Copy), ("", Remove),
                ("Translate...", Move2D), ("", Rotate2D),
                ("", Flip2D), ("", Scale2D), ("!", Array2D),
                ("Boolean...", Union2D),
                ("", Intersection), ("", Difference), ("!", Fragments),
                ]

    @classmethod
    def fancy_menu_name(self):
        return 'WP By Coordinates'

    @classmethod
    def fancy_tree_name(self):
        return 'WorkPlane'


data0 = (('center', VtableElement('pts1', type='string',
                                  guilabel='center',
                                  default="",
                                  tip="center of WP")),
         ('ax1', VtableElement('ax1', type='string',
                               guilabel='1st axis',
                               default="",
                               tip="point on the 1st axis")),
         ('ax2', VtableElement('ax2', type='string',
                               guilabel='point on surface',
                               default="",
                               tip="point on the surface")),
         ('flip1', VtableElement('flip1', type='bool',
                                 guilabel='flip 1st axis',
                                 default=False,
                                 tip="flip 1st axis")),
         ('flip2', VtableElement('flip2', type='bool',
                                 guilabel='flip 2nd axis',
                                 default=False,
                                 tip="flip 2nd axis")), )


class WorkPlaneByPoints(GeomPB):
    vt = Vtable(data0)

    def get_possible_child(self):
        return [Point2D, Line2D, Circle2D, Arc2D, Rect2D, Polygon2D, Spline2D,
                Move2D, Rotate2D, Flip2D, Scale2D, Array2D,
                Union2D, Intersection, Difference, Fragments, Copy, Remove,
                CreateLine, CreateSurface]

    def get_possible_child_menu(self):
        return [("", Point2D), ("", Line2D), ("", Circle2D), ("", Arc2D),
                ("", Rect2D), ("", Polygon2D), ("", Spline2D),
                ("", CreateLine), ("", CreateSurface),
                ("", Copy), ("", Remove),
                ("Translate...", Move2D), ("", Rotate2D),
                ("", Flip2D), ("", Scale2D), ("!", Array2D),
                ("Boolean...", Union2D),
                ("", Intersection), ("", Difference), ("!", Fragments),
                ]

    @classmethod
    def fancy_menu_name(self):
        return 'WP By Points'

    @classmethod
    def fancy_tree_name(self):
        return 'WorkPlane'


cad_fix_cb = [None, None, 36, {"col": 4,
                               "labels": ["Degenerated",
                                          "SmallEdges",
                                          "SmallFaces",
                                          "sewFaces"]}]
cad_fix_tol = ["Tolerance", 1e-8, 300, None]
cad_fix_elp0 = [cad_fix_cb, cad_fix_tol]

cad_fix_elp = [None, None, 27, [
    {"text": "import fixer"}, {"elp": cad_fix_elp0}]]


class healCAD(GeomPB):
    vt = Vtable(tuple())

    def panel1_param(self):
        from wx import BU_EXACTFIT
        b1 = {"label": "S", "func": self.onBuildBefore,
              "noexpand": True, "style": BU_EXACTFIT}
        b2 = {"label": "R", "func": self.onBuildAfter,
              "noexpand": True, "style": BU_EXACTFIT}
        wc = "ANY|*|STEP|*.stp|IGES|*.igs"

        ll = [[None, None, 241, {'buttons': [b1, b2],
                                 'alignright':True,
                                 'noexpand': True}, ],
              ["Entity", "", 0, {}],
              cad_fix_cb,
              cad_fix_tol]

        return ll

    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        v["fix_entity"] = ""
        v["use_fix"] = False
        v["use_fix_param"] = [True] * 4
        v["use_fix_tol"] = 1e-8
        return v

    def get_panel1_value(self):
        return [None, self.fix_entity, self.use_fix_param, self.use_fix_tol]

    def preprocess_params(self, engine):
        return

    def import_panel1_value(self, v):
        self.fix_entity = str(v[1])
        self.use_fix_param = [x[1] for x in v[2]]
        self.use_fix_tol = float(v[3])

    def panel1_tip(self):
        return [None, None, None, None]

    def add_geom_sequence(self, geom):
        gui_name = self.fullname()
        gui_param = (self.fix_entity, self.use_fix_param, self.use_fix_tol)
        geom_name = self.__class__.__name__
        geom.add_sequence(gui_name, gui_param, geom_name)

    @classmethod
    def fancy_menu_name(self):
        return "healCAD"

# we make BREP import separately so that we can add Brep specific
# interface later....


class CADImport(GeomPB):
    vt = Vtable(tuple())

    def panel1_param(self):
        from wx import BU_EXACTFIT
        b1 = {"label": "S", "func": self.onBuildBefore,
              "noexpand": True, "style": BU_EXACTFIT}
        b2 = {"label": "R", "func": self.onBuildAfter,
              "noexpand": True, "style": BU_EXACTFIT}
        wc = "ANY|*|STEP|*.stp|IGES|*.igs"

        ll = [[None, None, 241, {'buttons': [b1, b2],
                                 'alignright':True,
                                 'noexpand': True}, ],
              ["File(STEP/IGES)", None, 45, {'wildcard': wc}],
              cad_fix_elp,
              [None, True, 3, {'text': "Highest Dim Only"}],
              ["Unit", None, 0, {}],
              [None, "Unit: "", M, MM, INCH, KM, CM, FT, MIL, UM...", 2, None],
              ]

        return ll

    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        v["cad_file"] = ""
        v["use_fix"] = False
        v["use_fix_param"] = [True] * 4
        v["use_fix_tol"] = 1e-8
        v["highestdimonly"] = True
        v["import_unit"] = ""
        return v

    def get_panel1_value(self):
        value = [self.use_fix, [self.use_fix_param, self.use_fix_tol]]
        return [None, self.cad_file, value,
                self.highestdimonly, self.import_unit, None]

    def preprocess_params(self, engine):
        return

    def import_panel1_value(self, v):
        self.cad_file = str(v[1])
        self.use_fix = v[2][0]
        self.use_fix_param = [x[1] for x in v[2][1][0]]
        self.use_fix_tol = float(v[2][1][1])
        self.highestdimonly = bool(v[3])
        self.import_unit = v[4]

    def panel1_tip(self):
        return [None, None, None]

    def add_geom_sequence(self, geom):
        gui_name = self.fullname()
        gui_param = (self.cad_file, self.use_fix, self.use_fix_param, self.use_fix_tol,
                     self.highestdimonly, self.import_unit)
        geom_name = self.__class__.__name__
        geom.add_sequence(gui_name, gui_param, geom_name)

    @classmethod
    def fancy_menu_name(self):
        return "STEP/IGS"

# we make BREP import separately so that we can add Brep specific
# interface later....


class BrepImport(GeomPB):
    vt = Vtable(tuple())

    def panel1_param(self):
        from wx import BU_EXACTFIT
        b1 = {"label": "S", "func": self.onBuildBefore,
              "noexpand": True, "style": BU_EXACTFIT}
        b2 = {"label": "R", "func": self.onBuildAfter,
              "noexpand": True, "style": BU_EXACTFIT}
        wc = "ANY|*|Brep|*.brep"
        ll = [[None, None, 241, {'buttons': [b1, b2],
                                 'alignright':True,
                                 'noexpand': True}, ],
              ["File(.brep)", None, 45, {'wildcard': wc}],
              cad_fix_elp,
              [None, True, 3, {'text': "Highest Dim Only"}],
              ]
        return ll

    def attribute_set(self, v):
        v = super(GeomPB, self).attribute_set(v)
        v["cad_file"] = ""
        v["use_fix"] = False
        v["use_fix_param"] = [True] * 4
        v["use_fix_tol"] = 1e-8
        v["highestdimonly"] = True
        return v

    def get_panel1_value(self):
        value = [self.use_fix, [self.use_fix_param, self.use_fix_tol]]
        return [None, self.cad_file, value, self.highestdimonly]

    def preprocess_params(self, engine):
        return

    def import_panel1_value(self, v):
        self.cad_file = str(v[1])
        self.use_fix = v[2][0]
        self.use_fix_param = [x[1] for x in v[2][1][0]]
        self.use_fix_tol = float(v[2][1][1])
        self.highestdimonly = bool(v[3])

    def panel1_tip(self):
        return [None, None, None]

    def add_geom_sequence(self, geom):
        gui_name = self.fullname()
        gui_param = (
            self.cad_file,
            self.use_fix,
            self.use_fix_param,
            self.use_fix_tol,
            self.highestdimonly,
        )
        geom_name = self.__class__.__name__
        geom.add_sequence(gui_name, gui_param, geom_name)

    @classmethod
    def fancy_menu_name(self):
        return 'Brep'
