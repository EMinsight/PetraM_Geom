import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/box2.brep') 
'''
mw.add('transfinite_edge', '7', nseg = 6)
mw.add('transfinite_edge', '3', nseg = 7, bump=0.2)
mw.add('transfinite_edge', '11', nseg = 7, progression=1.03)
mw.add('freeface', '4', resolution=8)
mw.add('extrude_face', '1', '4', '3', nlayers=10)
#mw.add('freevolume', '1', resolution=7)
mw.add('freevolume', '2', resolution=7)
'''
mw.add('freevolume', '2', resolution=7)
mw.add('extrude_face', '1', '5', '2', nlayers=10)

mw.generate(dim = 3)# finalize=True)
mw.save_mesh_all('tmp')


