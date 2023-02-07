#ifdef FFTW
module FFTW3
  use, intrinsic :: iso_c_binding
#ifdef LIBPBSA
  include 'fftw3.f03'
#else
#ifdef MPI
  include 'fftw3-mpi.f03'
#else
  include 'fftw3.f03'
#endif
#endif
end module FFTW3
#endif
