program con_int
  implicit none
  !
  ! See JCP 104, 5517 (1996).
  !
  integer, parameter :: ndim=2
  real(kind=kind(1.d0)), dimension(ndim) :: xmin,xmax,dx
  integer, dimension(ndim)               :: nx
  real(kind=kind(1.d0)) :: alpha,beta,gamma,kx,ky,delta,x1,x2,x3,mx,my
  namelist /DAT/ alpha,beta,gamma,kx,ky,delta,x1,x2,x3,mx,my,xmin,xmax,nx
  character(len=2)                 :: aa
  real(kind=kind(1.d0)), parameter :: fstoau=41.34137221718d0
  real(kind=kind(1.d0)), parameter :: amutoau=1822.888515d0
  real(kind=kind(1.d0)), parameter :: zero=0.d0, due=2.d0, quattro=4.d0
  real(kind=kind(1.d0)) :: x,y,h11,h12,h22,eexp
  integer :: i,j,k
  !
  read(5,DAT)
  !
  do i=1,ndim
     dx(i) = (xmax(i)-xmin(i))/nx(i)
  end do

     do j=1,nx(2)+1
        y=xmin(2)+(j-1)*dx(2)
        do k=1,nx(1)+1
           x=xmin(1)+(k-1)*dx(1)
           h11  = (kx*(x-x1)**2 + ky*y*y)/due
           h22  = (kx*(x-x2)**2 + ky*y*y)/due + delta
           eexp = exp(-alpha*(x-x3)**2 - beta*y*y)
           h12  = gamma*y*eexp
           write(6,'(6f15.10)') x,y,h11,h12,h22
        end do
     end do
  !
  !
end program con_int
