module fft_drv
  use quant_util, only: discr, dpr, dpc
  use singleton, only: fftn
#ifdef FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
  type (C_PTR), save :: planf, planb
  real(kind=dpr), save :: fftw3_fac
#endif

contains
  !-----------------------------------------------------------------------------
  subroutine fft_ini(d,rwx,nrx)
    implicit none
    type (discr), intent (in) :: d
    integer, dimension(:), intent (inout) :: nrx 
    complex(kind=dpc), dimension(:), intent (inout) :: rwx 
    integer :: nr1
#ifdef FFTW3
    nr1 = d%nrt(d%ndim)
    nrx(1:d%ndim) = d%nr(d%ndim:1:-1) + 1
    planf = fftw_plan_dft(d%ndim,nrx,rwx,rwx,FFTW_FORWARD,FFTW_MEASURE)
    planb = fftw_plan_dft(d%ndim,nrx,rwx,rwx,FFTW_BACKWARD,FFTW_MEASURE)
    fftw3_fac = 1.0_dpr/sqrt(real(nr1,kind=dpr))
#else
    nrx(1:d%ndim) = d%nr(1:d%ndim) + 1
#endif
  end subroutine fft_ini
  !-----------------------------------------------------------------------------
  subroutine fft_end
    implicit none
#ifdef FFTW3
    call fftw_destroy_plan(planf)
    call fftw_destroy_plan(planb)
#endif
  end subroutine fft_end
  !-----------------------------------------------------------------------------
  subroutine fft_fwrd(rwx,nrx)
    implicit none
    integer, dimension(:),           intent (in)    :: nrx 
    complex(kind=dpc), dimension(:), intent (inout) :: rwx 
#ifdef FFTW3
    call fftw_execute_dft(planf,rwx,rwx)
    rwx = rwx*fftw3_fac
#else
    call fftn(rwx,nrx)
#endif
  end subroutine fft_fwrd
  !-----------------------------------------------------------------------------
  subroutine fft_bwrd(rwx,nrx)
    implicit none
    integer, dimension(:),           intent (in)    :: nrx 
    complex(kind=dpc), dimension(:), intent (inout) :: rwx 
#ifdef FFTW3
    call fftw_execute_dft(planb,rwx,rwx)
    rwx = rwx*fftw3_fac
#else
    call fftn(rwx,nrx,inv=.true.)
#endif
  end subroutine fft_bwrd
end module fft_drv
