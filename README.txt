gfortran correccion_scalar.f90 $HOME/looptools/lib64/libooptools.a -o correccion_scalar

gfortran -O3 polylog_module.f90 main.f90 -o main

gfortran -O3 polylog_module.f90 correccion_scalar.f90 $HOME/looptools/lib64/libooptools.a -o correccion_scalar



si no funciona, usar:

cd ~/LoopTools-2.15
make clean
./configure
make