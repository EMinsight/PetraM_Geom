'''
 OCC heal Shape

   the healing process here is genearated based on OCC_Internals::_healShape
   in GMSH, but re-written in Python
 
'''
import time

from petram.geom.occ_cbook import *

from OCC.Core.ShapeBuild import ShapeBuild_ReShape
from OCC.Core.ShapeExtend import (ShapeExtend_OK,
                                  ShapeExtend_DONE1,
                                  ShapeExtend_DONE2,
                                  ShapeExtend_DONE3,
                                  ShapeExtend_DONE4,
                                  ShapeExtend_DONE5,
                                  ShapeExtend_DONE6,
                                  ShapeExtend_DONE7,
                                  ShapeExtend_DONE8,
                                  ShapeExtend_FAIL1,
                                  ShapeExtend_FAIL2,
                                  ShapeExtend_FAIL3)

from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.BRepLib import breplib_OrientClosedSolid


def _fix_Degenerated(shape, verbose=False):
    if verbose:
        print(" - Fixing degenerated edges and faces")

    rebuild = ShapeBuild_ReShape()
    for edge in iter_shape(shape, 'edge'):
        if BRep_Tool.Degenerated(edge):
            rebuild.Remove(edge)
    shape = rebuild.Apply(shape)

    mapper = TopTools_IndexedMapOfShape()
    topexp_MapShapes(shape, TopAbs_FACE, mapper)

    rebuild = ShapeBuild_ReShape()
    sff = ShapeFix_Face()
    for face in iter_shape(shape, 'face'):
        sff = ShapeFix_Face(face)
        sff.SetFixAddNaturalBoundMode(True)
        sff.SetFixSmallAreaWireMode(True)
        sff.Perform()

        if (sff.Status(ShapeExtend_DONE1) or sff.Status(ShapeExtend_DONE2) or
            sff.Status(ShapeExtend_DONE3) or sff.Status(ShapeExtend_DONE4) or
            sff.Status(ShapeExtend_DONE5)):
            pass

        if verbose:
            print(" . Repaired face ", mapper.FindIndex(face))
            if sff.Status(ShapeExtend_DONE1):
                print(" . Some wires are fixed")
            elif sff.Status(ShapeExtend_DONE2):
                print(" . Orientation of wires fixed")
            elif sff.Status(ShapeExtend_DONE3):
                print(" . Missing seam added")
            elif sff.Status(ShapeExtend_DONE4):
                print(" . Small area wire removed")
            elif sff.Status(ShapeExtend_DONE5):
                print(" . Natural bounds added")

        newface = sff.Face()
        rebuild.Replace(face, newface)
    shape = rebuild.Apply(shape)

    rebuild = ShapeBuild_ReShape()
    for edge in iter_shape(shape, 'edge'):
        if BRep_Tool.Degenerated(edge):
            rebuild.Remove(edge)
    shape = rebuild.Apply(shape)

    return shape

def _fix_SmallEdges(shape, verbose=False, tolerance=1e-6):
    if verbose:
        print(" - Fixing small edges")

    mapper = TopTools_IndexedMapOfShape()
    topexp_MapShapes(shape, TopAbs_WIRE, mapper)

    sfw = ShapeFix_Wire()
    rebuild = ShapeBuild_ReShape()
    for oldwire, face in iterdouble_shape(shape, 'wire'):
        sfw = ShapeFix_Wire(oldwire, face, tolerance)

        sfw.SetModifyTopologyMode(True)
        sfw.SetClosedWireMode(True)

        replace = False
        replace = (sfw.FixReorder() or replace)
        replace = (sfw.FixConnected() or replace)

        num_fixsmall = sfw.FixSmall(False, tolerance)
        if num_fixsmall > 0:
            if not (sfw.StatusSmall(ShapeExtend_FAIL1) or
                    sfw.StatusSmall(ShapeExtend_FAIL2) or
                    sfw.StatusSmall(ShapeExtend_FAIL3)):
                if verbose:
                    print(" . Fixed small edge in wire %d",
                          mapper.FindIndex(oldwire))
            replace = True
        elif sfw.StatusSmall(ShapeExtend_FAIL1):
            print('\n'.join(["Failed to fix small edge in wire %d, edge cannot be checked ",
                             "(no 3d curve and no pcurve)"]),
                  mapper.FindIndex(oldwire))

        elif sfw.StatusSmall(ShapeExtend_FAIL2):
            print('\n'.join(["Failed to fix small edge in wire %d, ",
                             "edge is null-length and has different vertives at begin and ",
                             "end, and lockvtx is True or ModifiyTopologyMode is False"]),
                  mapper.FindIndex(oldwire))

        elif sfw.StatusSmall(ShapeExtend_FAIL3):
            print("Failed to fix small edge in wire, CheckConnected has failed",
                  mapper.FindIndex(oldwire))

        replace = (sfw.FixEdgeCurves() or replace)
        replace = (sfw.FixDegenerated() or replace)
        replace = (sfw.FixSelfIntersection() or replace)
        replace = (sfw.FixLacking(True) or replace)

        if replace:
            newwire = sfw.Wire()
            rebuild.Replace(oldwire, newwire)

    shape = rebuild.Apply(shape)

    rebuild = ShapeBuild_ReShape()
    for edge in iter_shape(shape, 'edge'):
        if BRep_Tool.Degenerated(edge):
            rebuild.Remove(edge)
    shape = rebuild.Apply(shape)

    sfwf = ShapeFix_Wireframe()
    sfwf.SetPrecision(tolerance)
    sfwf.Load(shape)
    sfwf.SetModeDropSmallEdges(True)

    num_fix = sfwf.FixSmallEdges()
    if num_fix > 0:
        if verbose:
            print(" - Fixing wire frames")
            if sfwf.StatusSmallEdges(ShapeExtend_OK):
                print(" . No small edges found")
            if sfwf.StatusSmallEdges(ShapeExtend_DONE1):
                print(" . Some small edges fixed")
            if sfwf.StatusSmallEdges(ShapeExtend_FAIL1):
                print(" . Failed to fix some small edges")
        shape = sfwf.Shape()
    return shape


def _fix_SmallFaces(shape, verbose=False, tolerance=1e-6):
    if verbose:
        print(" - Fixing small faces")
    sffsm = ShapeFix_FixSmallFace()
    sffsm.Init(shape)
    sffsm.SetPrecision(tolerance)
    sffsm.Perform()

    shape = sffsm.FixShape()
    return shape


def _sew_Faces(shape, verbose=False, tolerance=1e-6):
    if verbose:
        print(" - Sew faces")

    sewedObj = BRepBuilderAPI_Sewing(tolerance)
    for face in iter_shape(shape, 'face'):
        sewedObj.Add(face)

    sewedObj.Perform()

    if not sewedObj.SewedShape().IsNull():
        shape = sewedObj.SewedShape()
    else:
        if verbose:
            print(" . Could not sew")

    return shape

def _make_Solids(shape, verbose=False, tolerance=1e-6):
    ms = BRepBuilderAPI_MakeSolid()

    count = 0
    for shell in iter_shape(shape, 'shell'):
        ms.Add(shell)
        count = count + 1
        
    shape = ms.Solid()
    if count == 0:
        print(" . Could not make solid (no shell)")
    else:
        ba = BRepCheck_Analyzer(shape)
        if ba.IsValid():
            sfs = ShapeFix_Shape()
            sfs.Init(ms)
            sfs.SetPrecision(tolerance)
            sfs.SetMaxTolerance(tolerance)
            sfs.Perform()

            shape = sfs.Shape()

            for solid in iter_shape(shape, solid):
                newsolid = solid
                breplib_OrientClosedSolid(newsolid)
                rebuild = ShapeBuild_ReShape()
                rebuild.Replace(solid, newsolid)
                shape = rebuild.Apply(shape, TopAbs_COMPSOLID)
        else:
            if verbose:
                print(" . Could not make solid")

    return shape

def make_new_compound(shape, builder=None):

    if builder is None:
        builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)

    ex1 = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex1.More():
        builder.Add(comp, ex1.Current())
        ex1.Next()
    ex1 = TopExp_Explorer(shape, TopAbs_FACE, TopAbs_SHELL)
    while ex1.More():
        builder.Add(comp, ex1.Current())
        ex1.Next()
    ex1 = TopExp_Explorer(shape, TopAbs_EDGE, TopAbs_WIRE)
    while ex1.More():
        builder.Add(comp, ex1.Current())
        ex1.Next()
    ex1 = TopExp_Explorer(shape, TopAbs_VERTEX, TopAbs_EDGE)
    while ex1.More():
        builder.Add(comp, ex1.Current())
        ex1.Next()
    return comp


def heal_shape(shape, scaling=1.0, fixDegenerated=False,
               fixSmallEdges=False, fixSmallFaces=False,
               sewFaces=False, makeSolids=False,
               verbose=False, tolerance=1e-8):
    
    names = ('compound', 'compsolid', 'solid',
             'shell', 'face', 'wire', 'edge', 'vertex')
    
    old_maps = prep_maps(shape, return_all=True, return_compound=True)
    print([old_maps[n].Size() for n in names])
        
    if scaling != 1.0:
        if verbose:
            print("Scaling geometry (factor = " + str(scaling) + ")")
        t = gp_Trsf()

        t.SetScaleFactor(scaling)
        trsf = BRepBuilderAPI_Transform(t)
        trsf.Perform(shape, False)

        assert trsf.IsDone(), "Can not perform scaling"

        shape = trsf.Shape()
        shape = make_new_compound(shape)        

    if verbose:
        t1 = time.perf_counter()
        print("Healing shapes (tolerance = ", str(tolerance) + ")")

    surfacecount = 0
    for face in iter_shape(shape, 'face'):
        system = GProp_GProps()
        brepgprop_SurfaceProperties(face, system)
        surfacecount += system.Mass()

    do_final = False
    if fixDegenerated:
        shape = _fix_Degenerated(shape, verbose=verbose)
        do_final = True

    if fixSmallEdges:
        shape = _fix_SmallEdges(shape, verbose=verbose, tolerance=tolerance)
        do_final = True        

    if fixSmallFaces:
        shape = _fix_SmallFaces(shape, verbose=False, tolerance=tolerance)
        do_final = True
        
    if sewFaces:
        shape = _sew_Faces(shape, verbose=False, tolerance=tolerance)
        do_final = True

    if do_final:
        rebuild = ShapeBuild_ReShape()
        for edge in iter_shape(shape, 'edge'):
            if BRep_Tool.Degenerated(edge):
                rebuild.Remove(edge)
        shape = rebuild.Apply(shape)

    if makeSolids:
        shape = _make_Solids(shape, verbose=False, tolerance=tolerance)

    surfacecount2 = 0
    for face in iter_shape(shape, 'face'):
        system = GProp_GProps()
        brepgprop_SurfaceProperties(face, system)
        surfacecount2 += system.Mass()

    new_maps = prep_maps(shape, return_all=True, return_compound=True)

    if verbose:
        t2 = time.perf_counter()
        print("Done healing shapes (Wall : ", str(t2-t1), " s")
        for name in names:
            n1 = old_maps[name].Size()
            n2 = new_maps[name].Size()
            print("    - " + name.capitalize() + ",\t" +
                  str(n1) + " --> " + str(n2))
        print("    - Total surface area,  \t" +
              str(surfacecount) +  " --> " + str(surfacecount2))

    return shape
