import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/box2.brep') 
#mw.add('cl', [2], 0.02)
mw.add('transfinite_edge', '7', nseg = 5)
mw.add('transfinite_edge', '3', nseg = 15, bump=0.2)
mw.add('transfinite_edge', '11', nseg = 15, progression=1.03)
mw.add('freeface', '4', resolution=5)
#mw.add('freeface', '2', resolution=5)
mw.add('copyface','4','3')
#mw.add('copyface', '2', '10')
#mw.add('freeface', '1,3,4,6', resolution=15)
#mw.add('freeface', 'remaining', resolution=10)
#mw.add('freevolume', '1,2', resolution=10)
mw.generate(dim = 3, finalize=True)
mw.save_mesh_all('tmp')


