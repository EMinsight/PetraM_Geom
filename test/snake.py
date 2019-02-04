import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/snake.brep') 
mw.add('transfinite_edge', '10', nseg = 20)
#print(gmsh.model.getEntities())
#gmsh.model.mesh.setTransfiniteCurve(10, 20)
#gmsh.model.mesh.generate(1)
#mw.add('freeface', '2, 10', resolution=3)
#mw.add('copyface', '2, 10', '4, 11')
#mw.add('revolve_face', '1, 2', '2, 10', '4, 11', nlayers=10)
#mw.add('revolve_face', '1', '2', '4', nlayers=10)
#mw.add('freevolume', 'remaining', resolution=7)
mw.generate(dim = 1)# finalize=True)
mw.save_mesh_all('tmp')


