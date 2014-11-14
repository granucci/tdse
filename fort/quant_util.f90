module quant_util
  implicit none
  !
  integer, parameter :: dpr=kind(1.d0)
  integer, parameter :: dpc=kind((1.0_dpr,1.0_dpr))
  !
  type laser
     integer :: lasercw
     integer :: lasertyp
     real(kind=dpr) :: omega0
     real(kind=dpr), dimension(3) :: E0
     real(kind=dpr), dimension(20) :: laserparm
     integer :: filed
  end type laser
  !
  type discr
     integer           :: ndim
     integer           :: nstati,tcycles,nprt,istati,filem,nprt_wp
     integer           :: filewp,iwrt,adiabatize
     logical           :: hcomplex,radiation
     real(kind=dpr)    :: time,ene_add,thrwp
     real(kind=dpr)    :: t0,ttot
     character(len=10) :: wpack
     integer, allocatable, dimension (:)        :: nr,nrt
     real(kind=dpr), allocatable, dimension (:) :: dr,rmin,rmax
     real(kind=dpr), allocatable, dimension (:) :: massa,omega,r0,p0
     integer, allocatable, dimension (:,:)          :: inddim
     complex (kind=dpc), allocatable, dimension (:) :: fftkin
     real (kind=dpc), allocatable, dimension (:)    :: fftkin2
     type(laser) :: lsr
  end type discr
  !
  type hamil
     real(kind=dpr), allocatable    :: h(:,:)
     real(kind=dpr), allocatable    :: dia(:,:,:)
     complex(kind=dpr), allocatable :: h_c(:,:)
     complex(kind=dpr), allocatable :: dia_c(:,:,:)
     real(kind=dpr), allocatable    :: e(:,:)
     complex(kind=dpc), allocatable :: uteu(:,:,:)
     real(kind=dpr), allocatable    :: edip(:,:)
  end type hamil
  !
  real (kind=dpr), parameter :: fstoau = 41.34137221718_dpr
  real (kind=dpr), parameter :: autocm = 219474.6313705_dpr
  !
  real (kind=dpr), parameter :: c_au = 137.035999074_dpr
  !
  integer, parameter :: mdim=50
  !
contains
  recursive subroutine indset(d)
    implicit none
    type(discr), intent (inout) :: d 
    integer, dimension (mdim), save :: indx
    integer, save                   :: k,nlev=0
    integer, dimension (d%ndim) :: nr1
    integer                     :: i,m,ndim

    ndim=d%ndim
    nr1(1:ndim)=d%nr(1:ndim)+1

    nlev = nlev + 1
    indx(nlev)=0
    do
      indx(nlev) = indx(nlev) + 1
      if (nlev < ndim) then
         call indset(d)
      else
         k = indx(1)
         do m=1,nlev-1
           k = k + (indx(m+1)-1)*iprod(1,m,nr1)
         end do
         d%inddim(k,1:ndim)=indx(1:ndim)
         !d%r_dim(k,1:ndim) = d%rmin(1:ndim) + (indx(1:ndim)-1)*d%dr(1:ndim)
      endif
      if (indx(nlev) >= nr1(nlev)) then
         nlev=nlev-1
         exit
      endif
    end do
  end subroutine indset
  !------------------------------------------------------
  subroutine ind_fft_kin(d)
    implicit none
    type(discr), intent (inout) :: d 
    complex (kind=dpc), dimension (d%ndim) :: dw
    complex (kind=dpc)                     :: fac
    real (kind=dpc), dimension (d%ndim)    :: dw2
    real (kind=dpc)                        :: fac2
    real (kind=dpr)                        :: pi,denom,due
    integer, dimension (d%ndim)            :: nr1,nr2,ix,indkin
    integer                                :: i,j,ndim,nrtot
    !
    ndim=d%ndim
    nrtot=d%nrt(ndim)
    nr1(1:ndim)=d%nr(1:ndim)+1
    nr2(1:ndim)=nr1(1:ndim)/2 
    !
    pi=4.0_dpr*atan(1.0_dpr)
    due=2.0_dpr
    fac=cmplx(0,1,kind=dpc)*d%time*due*pi**2
    fac2=due*pi**2
    do i=1,ndim
      denom=d%massa(i)*(nr1(i)*d%dr(i))**2
      dw(i)=fac/denom
      dw2(i)=fac2/denom
    end do


    do i=1,nrtot
       ix=d%inddim(i,:)
       do j=1,ndim
         if (ix(j) > nr2(j)) then
            indkin(j)=(ix(j)-nr1(j)-1)**2
         else
            indkin(j)=(ix(j)-1)**2
         endif
       end do
       d%fftkin(i)=exp(-sum(dw(1:ndim)*indkin(1:ndim)))
       d%fftkin2(i)=sum(dw2(1:ndim)*indkin(1:ndim))
    end do
  end subroutine ind_fft_kin
  !------------------------------------------------------
  integer function iprod(i,j,a)
    implicit none
    integer, intent (in)                :: i,j
    integer, dimension (*), intent (in) :: a
    integer :: m

    iprod=1
    do m=i,j
      iprod=iprod*a(m)
    end do
  end function iprod
  !------------------------------------------------------
  real (kind=dpr) function aprod(i,j,a)
    implicit none
    integer, intent (in)                       :: i,j
    real(kind=dpr), dimension (*), intent (in) :: a
    integer :: m

    aprod=1.0_dpr
    do m=i,j
      aprod=aprod*a(m)
    end do
  end function aprod
end module quant_util
