module quant_efield
  use quant_util, only: dpr,dpc,discr,hamil,fstoau,autocm,c_au
  implicit none
  public :: efield, pulse_ini
  private
  real (kind=dpr), save :: pi,halfpi
  !
contains
  !-------------------------------------------------------------
  subroutine efield(d,hx)
    implicit none
    type (discr), intent (in)    :: d
    type (hamil), intent (inout) :: hx
    real (kind=dpr), dimension (d%ndim)               :: r
    real (kind=dpr), dimension (d%nstati)             :: e
    real (kind=dpr), dimension (3*d%nstati)           :: rw
    real (kind=dpr), dimension (d%nstati,d%nstati)    :: ah,u,ut
    complex (kind=dpc), dimension (d%nstati,d%nstati) :: ahc,uc,utc
    complex (kind=dpc), dimension (d%nstati,d%nstati) :: ute
    complex (kind=dpc), dimension (d%nstati*d%nstati) :: w
    complex (kind=dpc) :: vx
    real (kind=dpr)    :: dh,dscr,e1,e2,num,den,ang,c,s,h11,h12,h22,pval
    integer            :: i,j,k,kj,nstati,ndim,nrtot,alpha,lw,info
    !
    nstati = d%nstati
    !
    call pulse(d,pval)
    !
    ndim  = d%ndim
    nrtot = d%nrt(ndim)
    vx=cmplx(0,1,kind=dpc)*d%time*0.5_dpr
    lw=nstati**2
    !
    if (d%hcomplex) then
       do i=1,nrtot
          kj=0
          do j=1,nstati
             do k=1,j
                kj=kj+1
                ahc(j,k)=hx%h_c(i,kj) - hx%edip(i,kj)*pval
                ahc(k,j)=conjg(hx%h_c(i,kj))
                !write(6,*) 'i,j,k,ahc(j,k)',i,j,k,ahc(j,k)
             end do
          end do
          call zheev('V','L',nstati,ahc,nstati,e,w,lw,rw,info)
          uc = conjg(transpose(ahc))
          !!write(6,"(a/,20(4f12.6/))") ' controllo: real(ahc)  =',((real(ahc(k,j)),j=1,nstati),k=1,nstati)
          !!write(6,"(a/,20(4f12.6/))") ' controllo: aimag(ahc) =',((aimag(ahc(k,j)),j=1,nstati),k=1,nstati)
          !
          utc = ahc
          do alpha=1,d%nstati
             ute(:,alpha)=utc(:,alpha)*exp(-vx*e(alpha))
          end do
          hx%uteu(i,:,:)=matmul(ute,uc)
          !
       end do
    else
       do i=1,nrtot
          kj=0
          do j=1,nstati
             do k=1,j
                kj=kj+1
                ah(j,k)=hx%h(i,kj) - hx%edip(i,kj)*pval
                ah(k,j)=ah(j,k)
             end do
          end do
          call eispd(nstati,nstati,ah,e,ah,u)
          u = transpose(ah)
          !
          ut = ah
          do alpha=1,d%nstati
             ute(:,alpha)=ut(:,alpha)*exp(-vx*e(alpha))
          end do
          hx%uteu(i,:,:)=matmul(ute,u)
          !
       end do
    endif
  end subroutine efield
  !-------------------------------------------------------------
  subroutine pulse(d,pval)
    implicit none
    type (discr), intent (in)      :: d
    real(kind=dpr), intent (inout) :: pval
    real(kind=dpr) :: t,cw,g,t0,tau,invl

    t = d%ttot - d%time*0.5_dpr

    if (d%lsr%lasercw == 0) then
       cw=cos(d%lsr%omega0*t)
    elseif (d%lsr%lasercw == 1) then
       cw=sin(d%lsr%omega0*t)
    elseif (d%lsr%lasercw > 1) then
       write(6,'(a,i8)') " unknown carrier wave.   lasercw=",d%lsr%lasercw
       stop 12
    endif

    invl=0.0_dpr
    if (d%lsr%lasertyp == 1) then
       invl=1.0_dpr
    elseif (d%lsr%lasertyp == 2) then
       g=(t-d%lsr%laserparm(1))/d%lsr%laserparm(2)
       g=exp(-g**2)
       invl=g
    elseif (d%lsr%lasertyp == 3) then
       t0=d%lsr%laserparm(1)
       tau=d%lsr%laserparm(2)
       if ((t >= t0 - tau) .and. (t <= t0 + tau)) then
          invl=cos(halfpi*(t-t0)/tau)
       endif
    elseif (d%lsr%lasertyp > 3) then
       write(6,'(a,i8)') " unknown lasertyp.   lasertyp=",d%lsr%lasertyp
       stop 12
    endif
    pval = cw*invl
  end subroutine pulse
  !-------------------------------------------------------------
  subroutine pulse_ini(d)
    implicit none
    type (discr), intent (inout)      :: d
    real(kind=dpr) :: fwhm,powau,powmw,width,ampl,ampl2,t0,automw,tau

    pi=4.0_dpr*atan(1.0_dpr)
    halfpi=pi*0.5_dpr
    automw=6.4364088024e9_dpr

    ampl2=sum(d%lsr%e0*d%lsr%e0)
    ampl=sqrt(ampl2)
    powau=c_au*ampl2/(8.0_dpr*pi)
    powmw=powau*automw
    if (d%lsr%lasertyp == 1) then
       write(6,"(a)") ' Radiation: continuous wave'
       write(6,"(a,2(d18.9,a))") ' Frequency =',d%lsr%omega0,' a.u.  ', &
         & d%lsr%omega0*autocm,' cm^(-1)'
       write(6,"(a,d18.9,a)") ' Amplitude =',ampl,' a.u.'
       write(6,"(a,2(d18.9,a))") ' Intensity =',powau,' a.u.  ',powmw,' MW/cm^2'
    elseif (d%lsr%lasertyp == 2) then
       t0=d%lsr%laserparm(1)
       fwhm=d%lsr%laserparm(2)
       d%lsr%laserparm(2) = fwhm/sqrt(2.0_dpr*log(2.0_dpr))
       width=sqrt(8.0_dpr*log(2.0_dpr))/d%lsr%laserparm(2)
       write(6,"(a)") ' Radiation: continuous wave with gaussian envelope'
       write(6,"(a,2(d18.9,a))") ' Carrier frequency =',d%lsr%omega0,' a.u.  ',&
         & d%lsr%omega0*autocm,' cm^(-1)'
       write(6,"(a,d18.9,a)") ' Maximum amplitude =',ampl,' a.u.'
       write(6,"(a,2(d18.9,a))") ' Mean time =',t0,' a.u.  ',t0/fstoau,' fs'
       write(6,"(a,2(d18.9,a))") ' FWHM      =',fwhm,' a.u.  ',fwhm/fstoau,' fs'
       write(6,"(a,2(d18.9,a))") ' Peak Intensity =',powau,' a.u.  ',powmw,' MW/cm^2'
       write(6,"(a,2(d18.9,a))") ' Frequency FWHM =',width,' a.u.  ', &
         & width*autocm,' cm^(-1)'
    elseif (d%lsr%lasertyp == 3) then
       t0=d%lsr%laserparm(1)
       tau=d%lsr%laserparm(2)
       write(6,"(a)") ' Radiation: continuous wave with cossq envelope'
       write(6,"(a,2(d18.9,a))") ' Carrier frequency =',d%lsr%omega0,' a.u.  ',&
         & d%lsr%omega0*autocm,' cm^(-1)'
       write(6,"(a,d18.9,a)") ' Maximum amplitude =',ampl,' a.u.'
       write(6,"(a,2(d18.9,a))") ' Mean time =',t0,' a.u.  ',t0/fstoau,' fs'
       write(6,"(a,2(d18.9,a))") ' FWHM      =',tau,' a.u.  ',tau/fstoau,' fs'
       write(6,"(a,2(d18.9,a))") ' Peak Intensity =',powau,' a.u.  ',powmw,' MW/cm^2'
    elseif (d%lsr%lasertyp > 3) then
       write(6,'(a,i8)') " unknown lasertyp.   lasertyp=",d%lsr%lasertyp
       stop 12
    endif

  end subroutine pulse_ini
  !
end module quant_efield
