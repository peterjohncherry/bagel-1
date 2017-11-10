../configure CXXFLAGS=" -Wextra -Wall -Wno-placement-new -Wno-strict-aliasing  -Wno-deprecated-declarations -Wno-return-type -Wno-sign-compare -Wno-unused -Werror -mavx" \
--enable-mkl \
--with-mpi="intel" \
--with-include="-I/home/peter/UTILS/BOOST/IMPI_2018-boost_1_63_0 -I/opt/intel/mkl/include"\
 LDFLAGS="-L/home/peter/UTILS/BOOST/IMPI_2018-boost_1_63_0/stage/lib -L/opt/intel/mkl/lib/intel64 -L/opt/intel/mkl/bin/ " \
 CC="/opt/intel/impi/2018.0.128/intel64/bin/mpicc" CC="/opt/intel/impi/2018.0.128/intel64/bin/mpicxx"\
