gmsh_init = False

try:
   import gmsh
   has_gmsh = True   
except ImportError:
   has_gmsh = False
   assert False, "gmsh api is not found"

if not gmsh_init:   
   gmsh.initialize()
gmsh_init = True

