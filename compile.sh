CC=$(which mpicc) CXX=$(which mpicxx) cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/library -DBUILD_DOC=OFF  -DTEST=ON -DBOOST_MP=OFF -DMAKE_PYAMI=ON -DPYTHON_LIBRARY_DIR=/path/to/install/python/module/ ..
