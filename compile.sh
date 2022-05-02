CC=$(which mpicc) CXX=$(which mpicxx) cmake -DCMAKE_INSTALL_PREFIX=/where/to/install/libami/install -DBUILD_DOC=ON  -DTEST=ON -DBOOST_MP=OFF ..
