program rd_resu
  use quant_util, only: dpr,dpc,discr
  implicit none
  character (len=120) :: filewp,card
  integer             :: narg,k,icod,ndim,nrtot,nstati,a,i
  real (kind=dpr)     :: ti,tini,tend,rmin
  type (discr)        :: d
  complex (kind=dpc), allocatable, dimension (:,:) :: rw
  real (kind=dpr), allocatable, dimension (:)      :: r
  !
  narg = command_argument_count()
  if (narg <= 1) then
     write(6,'(a)') ' Uso: rd_resu <file_wp> time_ini [time_end]'
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
  allocate (d%dr(ndim),d%rmin(ndim),d%nrt(ndim))
  read(1) d%rmin,d%dr,d%nrt,d%massa,d%omega,d%r0,d%p0
  nrtot  = d%nrt(ndim)
  allocate (d%inddim(nrtot,ndim))
  read(1) d%inddim
  !
  allocate (rw(nrtot,nstati),r(ndim))
  !
  if (narg >= 2) then
     call get_command_argument(2, card)
     read(card,*)tini
     tend=tini+10000.0_dpr
     if (narg >= 3) then
        call get_command_argument(3, card)
        read(card,*)tend
     endif
     do 
       read(1,iostat=icod)k,ti,rw
       if (icod /= 0) exit
       if (ti >= tini .and. ti < tend) then
          write(6,'(a,i8,a,f15.6)')'#   step=',k,'  time (fs)=',ti
          do i=1,nrtot
            r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
            write(6,'(41f20.10)')r,(rw(i,a),a=1,nstati)
          end do
          if (narg == 2) then
             exit
          endif
          write(6,*)
       endif
     end do
  endif
  
end program rd_resu
