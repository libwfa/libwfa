# libwfa
Wave-function analysis tool library for quantum chemical applications

Code and ideas by F. Plasser, M. Wormit, S. A. Bäppler, B. Thomitzni, and A. Dreuw

Please contact one of the authors if you are interested in interfacing libwfa to your quantum chemistry code. See LICENSE for licensing information.

To compile the standalone version, type:
~~~~
./configure
cd build
make
~~~~

Testing requires the libtest library, which is currently not included. Compile libtest, then type:
~~~~
cd libwfa/tests
../build/tests/libwfa_tests
~~~~

