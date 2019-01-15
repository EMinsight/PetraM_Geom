import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/torus2.brep') 
#mw.add('transfinite_edge', '81', nseg = 7)
#mw.add('freeface', '17', resolution=3)
#mw.add('copyface', '17', '45')

mw.add('freeface', '8, 17', resolution=3)

mw.add('extrude_face', '2, 3', '8, 17', '10, 18', nlayers=10)
mw.add('extrude_face', '4, 6', '10, 18', '25, 40', nlayers=4)
mw.add('extrude_face', '5, 7', '25, 40', '33, 45', nlayers=10)
'''
#mw.add('freevolume', 'remaining', resolution=7)
#mw.add('freeface', '6', resolution=3)
'''
mw.add('freevolume', 'remaining', resolution=4)
mw.generate(dim = 3, finalize=True)
mw.save_mesh_all('tmp')


