# libwfa
Wave-function analysis tool library for quantum chemical applications

Code and ideas by F. Plasser, M. Wormit, S. A. Bäppler, B. Thomitzni, and A. Dreuw

Further contributions by F. Chen, P. Pokhilko, and A. I. Krylov

Please contact one of the authors if you are interested in interfacing libwfa to your quantum chemistry code. See LICENSE for licensing information.

To compile the standalone version, type:
~~~~
./configure
cd build
make
~~~~

Testing of the standalone library is only possible via the libtest library, which is part of Q-Chem. Compile libtest, then type:
~~~~
cd libwfa/tests
../build/tests/libwfa_tests
~~~~

If you are using libwfa as part of OpenMolcas, you can use
~~~~
pymolcas verify 835
~~~~
for verification of the OpenMolcas/WFA interface.
