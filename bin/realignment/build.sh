cd realign
gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h

#bioconda
#meta.yaml
#  build:
#    - {{ compiler('c') }}
#    - {{ compiler('cxx') }}
#  host:
#    - boost-cpp

#build.sh
#$CXX -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
#$CXX -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp -I $PREFIX/include -L $PREFIX/lib

#ONT
#cd realign
#$CC -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h
