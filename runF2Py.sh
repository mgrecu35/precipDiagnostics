gfortran -c -fPIC hydromet/MCRPparameters.f90
gfortran -c -fPIC hydromet/gammaf.f90
gfortran -c -fPIC hydromet/heymsplatt.f90
gfortran -c -fPIC hydromet/watoptic.f90
gfortran -c -fPIC hydromet/mie_sphere.f90
gfortran -c -fPIC hydromet/absorb.f
gfortran -c -fPIC radtran/*.f90
gfortran -c -fPIC plug-ins/pwvsat.f90
gfortran -c -fPIC radtran/band_dble.f90
gfortran -fPIC -shared -o libsdsu_dep.so MCRPparameters.o gammaf.o heymsplatt.o watoptic.o mie_sphere.o absorb.o band_dble.o  
python -m numpy.f2py -c -I/Users/mgrecu/precipDiagnostics/ -L/Users/mgrecu/precipDiagnostics -m sdsu_tables hydromet/microp_set.f90  hydromet/mie_cldw.f90 radtran/emissivity-sp.f pwvsat.o plug-ins/qv2rh.f90 plug-ins/q2Wkg.f90 radtran/radtran_tau_dble.f radtran/rosen.f fortran_rad_tran.f90 -lsdsu_dep 
