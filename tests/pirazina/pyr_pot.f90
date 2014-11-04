program pyr_pot
  implicit none
  integer, parameter :: ndim=3
  real(kind=kind(1.d0)), dimension (ndim) :: k1, k2, w
  real(kind=kind(1.d0)) :: lambda, e1, e2
  real(kind=kind(1.d0)), parameter :: autoev=27.2113957d0
  character (len=1), dimension (0:9), parameter :: digit =&
      &(/'0','1','2','3','4','5','6','7','8','9'/)
  real(kind=kind(1.d0)), dimension(ndim) :: xmin,xmax
  integer, dimension(ndim)               :: nx
  namelist /DAT/ xmin,xmax,nx
  real(kind=kind(1.d0)), dimension (ndim) :: r,h,dx
  integer :: i,j,k

  xmin=-10.d0
  xmax=10.d0
  nx=2
  read(5,DAT)

!-- parametri originali in eV
  w(1) = 0.126d0
  w(2) = 0.074d0
  w(3) = 0.118d0
  k1(1) =  0.037d0
  k1(2) = -0.105d0
  k1(3) =  0.d0
  k2(1) = -0.254d0
  k2(2) =  0.149d0
  k2(3) =  0.d0
  lambda = 0.262d0
  e1 = 3.94d0
  e2 = 4.84d0

!-- passaggio a a.u.
  w=w/autoev
  k1=k1/autoev
  k2=k2/autoev
  lambda=lambda/autoev
  e1=e1/autoev
  e2=e2/autoev

  !write(6,*) 'm1=',1.d0/w(1)
  !write(6,*) 'm2=',1.d0/w(2)
  !write(6,*) 'm3=',1.d0/w(3)
  
  do i=1,ndim
     dx(i) = (xmax(i)-xmin(i))/nx(i)
  end do

  do i=1,nx(3)+1
     r(3)=xmin(3)+(i-1)*dx(3)
     do j=1,nx(2)+1
        r(2)=xmin(2)+(j-1)*dx(2)
        do k=1,nx(1)+1
           r(1)=xmin(1)+(k-1)*dx(1)
           call hdia(r,h)
           write(6,'(6f15.9)') r,h
        end do
     end do
  end do

contains
  !----------------------------------------------------------------
  subroutine hdia(x,val)
    implicit none
    real(kind=kind(1.d0)), dimension(ndim), intent (in)    :: x
    real(kind=kind(1.d0)), dimension(ndim), intent (inout) :: val
    real(kind=kind(1.d0)) :: wxx

    wxx=sum(w(:)*x(:)**2)*0.5d0

    val(1) = e1 + wxx + sum(k1*x)
    val(2) = lambda*x(3)
    val(3) = e2 + wxx + sum(k2*x)

  end subroutine hdia
  
end program pyr_pot
