## Petra-M geometry

Geometry editor interface

### Build/Install
   These steps are common to all PetraM module installation

   1) decide install destination and set enviromental variables, PetraM
   
      example)
         - export PetraM=/usr/local/PetraM  (public place)
         - export PetraM=/Users/shiraiwa/sandbox_root/PetraM (user directory)

   2) add PYTHONPATH
      - if you dont have, make python package directory 
           mkdir -p $PetraM/lib/python2.7/site-packages
      - set PYTHONPATH
           export PYTHONPATH=$PetraM/lib/python2.7/site-packages:$PYTHONPATH

   4) prepare Makefile
       use a file in Makefile_templates.
       Makefile_default may work well.

   5) build
      - make 

   6) install
      - make install
