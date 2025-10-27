#!bin/bash
cd ./src
rm *.so *.mod *.o
gfortran -c Constants.f90
f2py -m Constants -c Constants.f90
#wait
cd ./Dynamics
f2py -m Dynamics_reverse -c Dynamics_reverse.f90 -I ../ --f90flags='-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
f2py -m Dynamics_forward -c Dynamics_forward.f90 -I ../ --f90flags='-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
#wait
cd ../Electron
f2py -m FS_electron_weno5 -c FS_electron_weno5.f90 -I ../ --f90flags='-Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
f2py -m FS_electron_fullhide -lgomp -c FS_electron_fullhide.f90 -I ../ --f90flags='-fopenmp -Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
#wait
cd ../Interpolation
rm *.so
f2py -m SED_interpolation -lgomp -c SED_interpolation.f90  -I ../ --f90flags='-fopenmp -Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
f2py -m SED_interpolation_structured -lgomp -c SED_interpolation_structured.f90  -I ../ --f90flags='-fopenmp -Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
#wait
cd ../Radiation
rm *.so
f2py -m Annihilation -lgomp -c Annihilation.f90 -I ../ --f90flags='-fopenmp -Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
f2py -m Seed_reverse -lgomp -c Seed_reverse.f90 -I ../ --f90flags='-fopenmp -Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
f2py -m SSC_spec -lgomp -c SSC_spec.f90 -I ../ --f90flags='-fopenmp -Ofast -march=native -funroll-loops -ffast-math -fno-signed-zeros -fno-trapping-math'
cd ..
cd ..
echo "Compile complete!"
