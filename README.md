
Los programas se deben inicializar en el sistema operativo Linux.
Para instalar la librería LoopTools se debe colocar el comando:

cd ~/LoopTools-2.15
make clean
./configure
make

Para correr el programa que se desee, se debe usar el comando:

gfortran -O3 polylog_module.f90 correccion_scalar.f90 $HOME/looptools/lib64/libooptools.a -o correccion_scalar

cambiando el nombre del programa al que se desee.
