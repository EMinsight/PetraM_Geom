import gmsh
import petram.mesh.gmsh_mesh_wrapper 
import petram.geom.geom_utils
reload(petram.geom.geom_utils)
reload(petram.mesh.gmsh_mesh_wrapper)
from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper

mw = GMSHMeshWrapper()
model1 = mw.load_brep('brep/revolve2.brep') 
#mw.add('transfinite_edge', '81', nseg = 7)
#mw.add('freeface', '17', resolution=3)
#mw.add('copyface', '17', '45')
mw.add('transfinite_edge', '10', nseg = 12)
mw.add('freeface', '6', resolution=3)
#mw.add('extrude_face', '1', '6', '5', nlayers=5)
#mw.add('copyface', '6', '5', axan = ([0, 1., 1], -60))
#mw.add('freevolume', 'remaining', resolution=4)
mw.add('revolve_face', '1', '6', '5', axan = ([0, 1., 1], -60))
'''
#mw.add('freevolume', 'remaining', resolution=7)
#mw.add('freeface', '6', resolution=3)
'''

mw.generate(dim = 3)#, finalize=True)
mw.save_mesh_all('tmp')


