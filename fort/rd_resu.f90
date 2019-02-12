program rd_resu
  use quant_util, only: dpr,dpc,discr,ind_fft_kin,aprod
  implicit none
  complex (kind=dpc), parameter :: zzero  = (0.0_dpr,0.0_dpr)
  character (len=120) :: filewp,card
  integer             :: narg,k,icod,ndim,nrtot,nstati,a,i,j
  real (kind=dpr)     :: ti,tini,tend,rmin,drtot,rdiss
  type (discr)        :: d
  complex (kind=dpc), allocatable, dimension (:,:) :: rw
  real (kind=dpr), allocatable, dimension (:)      :: r
  real (kind=dpr), allocatable, dimension (:)      :: ek,ep,pop
  real (kind=dpr), allocatable, dimension (:,:)    :: e_pes,ravg,rsd
  !
  narg = command_argument_count()
  if (narg <= 1) then
     write(6,'(a)') ' Uso: rd_resu <file_wp> time_ini [time_end, rdiss]'
     stop 12
  endif
  !
  call get_command_argument(1, filewp)
  open(1,file=filewp,form='unformatted',status='old')
  rewind 1
  read(1) d%ndim,d%nstati,d%istati,d%hcomplex,d%adiabatize
  ndim   = d%ndim
  nstati = d%nstati
  allocate (d%r0(ndim),d%p0(ndim),d%massa(ndim),d%omega(ndim))
  allocate (d%dr(ndim),d%rmin(ndim),d%nrt(ndim),d%nr(ndim))
  read(1) d%rmin,d%dr,d%nrt,d%massa,d%omega,d%r0,d%p0
  nrtot  = d%nrt(ndim)
  allocate (d%inddim(nrtot,ndim),d%fftkin(nrtot),d%fftkin2(nrtot))
  allocate (e_pes(nrtot,nstati))
  read(1) d%inddim,e_pes
  !
  do i=1,ndim
    if (i == 1) then
       d%nr(1) = d%nrt(1) - 1
    else
       d%nr(i) = d%nrt(i)/d%nrt(i-1) - 1
    endif
  end do
  call ind_fft_kin(d)
  drtot=aprod(1,ndim,d%dr)
  !
  allocate (rw(nrtot,nstati),r(ndim))
  allocate (ek(nstati),ep(nstati),pop(nstati))
  allocate (rsd(ndim,nstati),ravg(ndim,nstati))
  !
  open(2,file='WPDATA',form='formatted',status='unknown')
  rewind 2
  !
  rdiss=-1.0_dpr
  tini=0.0_dpr
  if (narg >= 2) then
     call get_command_argument(2, card)
     read(card,*)tini
     tend=tini+10000.0_dpr
     if (narg >= 3) then
        call get_command_argument(3, card)
        read(card,*)tend
        if (narg >= 4) then
           call get_command_argument(4, card)
           read(card,*)rdiss
        endif
     endif
  endif
  do 
    read(1,iostat=icod)k,ti,rw
    if (icod /= 0) exit
    if (ti >= tini .and. ti < tend) then
       call ekin(rw,ek,ep,pop)
       call rcalc(rw,pop,ravg,rsd)
       write(2,'(a,i8,a,f15.6,a,e18.9)')'   step=',k,'  time (fs)=',ti,' ekin=',sum(ek)
       write(2,'(a)') ' State      pop       ekin           epot                  <r>            r_sd'
       do i=1,nstati
          write(2,'(i6,f12.6,50E18.9)') i,pop(i),ek(i),ep(i),(ravg(j,i),rsd(j,i),j=1,ndim)
       end do
       write(6,'(a,i8,a,f15.6,a,20f12.6)')'#   step=',k,'  time (fs)=',ti,' pop=',pop
       do i=1,nrtot
         r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
         write(6,'(41f20.10)')r,(rw(i,a),a=1,nstati)
       end do
       if (rdiss > 0.0_dpr .and. ndim == 1) then
          call diss(rw,rdiss,ek,ep,pop,ravg,rsd)
          write(2,'(a,f12.6,a)') ' At dissociation (r >',rdiss,')'
          write(2,'(a)') ' State      pop       ekin           epot                  <r>            r_sd'
          write(2,'(i6,f12.6,4e18.9)') (i,pop(i),ek(i),ep(i),ravg(1,i),rsd(1,i),i=1,nstati)
       endif
       if (narg == 2) then
          exit
       endif
       write(6,*)
    endif
  end do
contains
  !****************************************************************************
  subroutine ekin(rw,ek,ep,pop)
    use singleton
    implicit none
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    real (kind=dpr), dimension (:), intent (inout)      :: ek,ep,pop
    complex (kind=dpc), dimension(size(rw,1)) :: rwfft
    integer, dimension (size(d%nr,1))         :: nrx
    integer :: alpha

    nrx=d%nr+1
    do alpha=1,d%nstati
       rwfft=rw(:,alpha)
       call fftn(rwfft,nrx)
       ek(alpha)=sum(d%fftkin2(:)*abs(rwfft(:))**2)*drtot
       ep(alpha)=sum(e_pes(:,alpha)*abs(rw(:,alpha))**2)*drtot
       pop(alpha)=sum(abs(rw(:,alpha))**2)*drtot
    end do

  end subroutine ekin
  !****************************************************************************
  subroutine rcalc(rw,pop,ravg,rsd)
    implicit none
    complex (kind=dpc), dimension (:,:), intent (in) :: rw
    real (kind=dpr), dimension (:), intent (in)      :: pop
    real (kind=dpr), dimension (:,:), intent (inout) :: ravg,rsd
    real (kind=dpr), dimension (nrtot)    :: r
    integer :: alpha,jdim

    do jdim=1,ndim
       r(1:nrtot) = d%rmin(jdim) + (d%inddim(1:nrtot,jdim)-1)*d%dr(jdim)
       do alpha=1,d%nstati
          ravg(jdim,alpha)=sum(r(:)*abs(rw(:,alpha))**2)*drtot
          rsd(jdim,alpha)=sum(r(:)**2*abs(rw(:,alpha))**2)*drtot
       end do
    end do

    do alpha=1,d%nstati
       if (pop(alpha) < 1.e-10_dpr) then
          ravg(1:ndim,alpha) = 0.0_dpr
          rsd(1:ndim,alpha) = 0.0_dpr
       else
          ravg(1:ndim,alpha) = ravg(1:ndim,alpha)/pop(alpha)
          rsd(1:ndim,alpha) = rsd(1:ndim,alpha)/pop(alpha)
       endif
    end do
    rsd=sqrt(rsd-ravg**2)
  end subroutine rcalc
  !****************************************************************************
  subroutine diss(rw,rdiss,ek,ep,pop,ravg,rsd)
    implicit none
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    real (kind=dpr), dimension (:), intent (inout)      :: ek,ep,pop
    real (kind=dpr), dimension (:,:), intent (inout)    :: ravg,rsd
    real (kind=dpr), intent (inout)                     :: rdiss
    real (kind=dpr), dimension (nrtot)    :: r
    integer :: alpha

    r(1:nrtot) = d%rmin(1) + (d%inddim(1:nrtot,1)-1)*d%dr(1)
    do alpha=1,d%nstati
       where (r < rdiss)
         rw(:,alpha) = zzero
       end where
    end do

    call ekin(rw,ek,ep,pop)
    call rcalc(rw,pop,ravg,rsd)

  end subroutine diss
  !****************************************************************************
end program rd_resu

