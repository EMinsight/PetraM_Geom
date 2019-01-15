import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/box2.brep') 
#mw.add('cl', [2], 0.02)
#mw.add('freeedge', [7], resolution=15)
mw.add('freeface', '2', resolution=5)
mw.add('copyface', '2', '5')
mw.add('copyface', '2', '10')
mw.add('freeface', '11,3,8,6', resolution=20)
mw.add('freevolume', '1,2', resolution=10)
mw.generate(dim = 3)# finalize=True)
mw.save_mesh_all('tmp')

