
Los programas se deben inicializar en el sistema operativo Linux.

Ademas, para poder correr los programas de fortran se debe tener instalada la libreria de LoopTools.

Para compilar el programa que se desee, se debe usar el comando:

gfortran -O3 polylog_module.f90 correccion_scalar.f90 $HOME/looptools/lib64/libooptools.a -o correccion_scalar

cambiando el nombre del programa y la ruta escrita.