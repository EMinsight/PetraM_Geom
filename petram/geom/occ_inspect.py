import numpy as np

from petram.geom.occ_cbook import *

def xyz2txt(c):
    return ", ".join([str(x) for x in (c.X(), c.Y(), c.Z())])

def shape_property_txt(bt, shape):
    if isinstance(shape, TopoDS_Vertex):
        pnt = bt.Pnt(shape)
        c_txt = xyz2txt(pnt)
        txt = ['Point:',
               '  Coords:\t' + c_txt]

    if isinstance(shape, TopoDS_Edge):
        loc = TopLoc_Location()
        curve, first, last = bt.Curve(shape, loc)
        is_closed = curve.IsClosed()
        is_periodic = curve.IsPeriodic()
        curve, kind = downcast_curve(curve)
        length = measure_edge_length(shape)
        txt = ['Curve:',
               '  Kind:\t' + kind,
               '  Length:\t' + str(length), 
               '  Parameter:\t' + str([first, last]),
               '  Closed:\t' + str(is_closed),
               '  Periodic:\t' + str(is_periodic)]

        if curve.IsKind('Geom_Conic'):
            txt_c = xyz2txt(curve.Location())
            txt_n = xyz2txt(curve.Axis().Direction())
            txt.extend(['  Center:\t' + txt_c,
                        '  Normal:\t' + txt_n])
        if curve.IsKind('Geom_Circle'):
            r = curve.Radius()
            txt.extend(['  Radius:\t' + str(r)])
        if curve.IsKind('Geom_Ellipse'):
            r1 = curve.MajorRadius()
            r2 = curve.MinorRadius()
            txt_f1 = xyz2txt(curve.Focus1())
            txt_f2 = xyz2txt(curve.Focus2())
            txt.extend(['  Radius1:\t' + str(r1),
                        '  Radius2:\t' + str(r2),
                        '  Focus1:\t' + txt_f1,
                        '  Focus2:\t' + txt_f2])
        if curve.IsKind('Geom_Hypabola'):
            r1 = curve.MajorRadius()
            r2 = curve.MinorRadius()
            txt_f1 = xyz2txt(curve.Focus1())
            txt_f2 = xyz2txt(curve.Focus2())
            txt.extend(['  Radius1\t:' + str(r1),
                        '  Radius2:\t' + str(r2),
                        '  Focus1:\t' + txt_f1,
                        '  Focus2:\t' + txt_f2])
        if curve.IsKind('Geom_Parabola'):
            txt_f = xyz2txt(curve.Focus())
            txt.extend(['  Focus:\t' + txt_f])
        if curve.IsKind('Geom_BSplineCurve'):
            txt.extend(['  Start:\t' + xyz2txt(curve.StartPoint()),
                        '  End:\t' + xyz2txt(curve.EndPoint()),
                        '  #Knots:\t' + str(curve.NbKnots()),
                        '  #Poles:\t' + str(curve.NbPoles())])
        if curve.IsKind('Geom_BezierCurve'):
            txt.extend(['  Start:\t' + xyz2txt(curve.StartPoint()),
                        '  End:\t' + xyz2txt(curve.EndPoint()),
                        '  #Poles:\t' + str(curve.NbPoles())])
        if curve.IsKind('Geom_TrimmedCurve'):
            txt.extend(['  Start:\t' + xyz2txt(curve.StartPoint()),
                        '  End:\t' + xyz2txt(curve.EndPoint()),])


    if isinstance(shape, TopoDS_Face):
        surf = bt.Surface(shape)
        u1, u2, v1, v2 = surf.Bounds()
        is_uperiodic = surf.IsUPeriodic()
        is_vperiodic = surf.IsVPeriodic()

        system = GProp_GProps()
        brepgprop_SurfaceProperties(shape, system)
        surfacecount = system.Mass()

        surf, kind = downcast_surface(surf)
        print(surf, kind)
        txt = ['Surface:',
               ' Kind:\t' + kind,
               ' Area:\t' + str(surfacecount),
               ' U-Parameter:\t' + str([u1, u2]),
               ' V-Parameter:\t' + str([v1, v2]),
               ' Periodic (U,V):\t' + str([is_uperiodic, is_vperiodic]),]

        if surf.IsKind('Geom_Plane'):
            a, b, c, d = surf.Coefficient()
            txt = ', '.join([str(x) for x in (a, b, c, d)])
            txt.extend(['  Coefficient:\t' + txt])

    if isinstance(shape, TopoDS_Solid):
        txt = ['',]

    return '\n'.join(txt)

def shape_inspector(shape, inspect_type, shapes):

    print("inspection ", shape, inspect_type, shapes)
    bt = BRep_Tool()

    ret = ''
    data = None

    if inspect_type == 'property':
        prop = [shape_property_txt(bt, s) for s in shapes]
        return ' \n\n'.join(prop), ''

    if inspect_type == 'smallface':
        args, topolist = shapes
        thr = args[0]
        nsmall, smax, faces, areas = check_shape_area(shape, thr,
                                                      return_area=True)
        gids = [topolist.find_gid(f) for f in faces]
        txt = ',\n'.join([str(gid) + " (area = "+str(a) + ")"
                          for gid, a in zip(gids, areas)])

        txt = txt + '\n smax = ' + str(smax)

        if nsmall != 0:
            data = gids
        return txt, data

    if inspect_type == 'shortedge':
        args, topolist = shapes
        thr = args[0]
        nsmall, lmax, edges, ll = check_shape_length(shape, thr,
                                                     return_area=True)
        gids = [topolist.find_gid(e) for e in edges]
        txt = ',\n'.join([str(gid) + " (L = "+str(l) + ")"
                          for gid, l in zip(gids, ll)])

        txt = txt + '\n smax = ' + str(lmax)

        if nsmall != 0:
            data = gids
        return txt, data

    if inspect_type == 'dist_p_s':
        # distance between point and surface

        pnt = bt.Pnt(shapes[0])
        p1 = np.array((pnt.X(), pnt.Y(), pnt.Z(),))
        surf = bt.Surface(shapes[1])

        pj = GeomAPI_ProjectPointOnSurf(pnt, surf)
        print("number of solution ", pj.NbPoints())

        pnt = pj.NearestPoint()
        p2 = np.array((pnt.X(), pnt.Y(), pnt.Z(),))

        print("number of solution ", pj.NbPoints())

        dist = np.sqrt(np.sum((p1 - p2)**2))
        ret = dist
    return ret, data

