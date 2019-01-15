import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/box_circle.brep') 
mw.add('cl', '5', 0.03)
#mw.add('freeedge', [7], resolution=15)
mw.add('freeedge', '2, 3', resolution=10)
mw.add('freeface', '1', resolution=5, maxsize=0.15)
mw.generate(finalize=True)
mw.save_mesh_all('tmp')


