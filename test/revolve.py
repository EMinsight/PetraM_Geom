import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/torus.brep') 
mw.add('freeface', '2, 10', resolution=3)
#mw.add('copyface', '2, 10', '4, 11')
mw.add('revolve_face', '1, 2', '2, 10', '4, 11', nlayers=10)
#mw.add('revolve_face', '1', '2', '4', nlayers=10)
'''
mw.add('transfinite_edge', '7', nseg = 6)
mw.add('transfinite_edge', '3', nseg = 7, bump=0.2)
mw.add('transfinite_edge', '11', nseg = 7, progression=1.03)
mw.add('freeface', '4', resolution=8)
mw.add('extrude_face', '1', '4', '3', nlayers=10)
#mw.add('freevolume', '1', resolution=7)
mw.add('freevolume', '2', resolution=7)
'''
#mw.add('freevolume', 'remaining', resolution=7)
mw.generate(dim = 2)# finalize=True)
mw.save_mesh_all('tmp')


