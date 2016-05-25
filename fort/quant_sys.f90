module quant_sys
  use quant_util, only: dpc,dpr,discr,hamil,fstoau,aprod
  implicit none

  public :: ini_wp, hmol
  private

  real (kind=dpr), parameter :: eps=1.e-8_dpr

contains
  !-------------------------------------------------------------
  subroutine ini_wp(d,rw)
    implicit none
    type (discr), intent (inout)                        :: d
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    complex (kind=dpc), parameter :: zzero = (0.0_dpr,0.0_dpr)
    complex (kind=dpc)  :: val
    character (len=120) :: fnam
    real (kind=dpr)     :: dum,fac2t,pi,rp,wpx
    real (kind=dpr), dimension (d%ndim) :: omega,fac1,fac2,r,rr
    integer             :: i,alpha,icod,ndim,nrtot
    !
    pi = atan(1.0_dpr)*4.0_dpr
    ndim=d%ndim
    !
    do i=1,ndim
      omega(i) = d%omega(i)
      fac1(i)  = d%massa(i)*omega(i)/2.0_dpr
      dum      = sqrt(d%massa(i)*omega(i)/pi)
      fac2(i)  = sqrt(dum)
    end do
    fac2t=aprod(1,ndim,fac2)
    !
    nrtot=d%nrt(ndim)
    if (d%wpack == 'gauss') then
       do alpha=1,d%nstati
          do i=1,nrtot
             val = zzero
             if (alpha == d%istati) then
                r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
                r(:) = r(:) - d%r0(:)
                !r(:) = d%r_dim(i,:) - d%r0(:)
                rr=r**2
                rp=dot_product(d%p0,r)
                dum=fac2t*exp(-dot_product(fac1,rr))
                val = dum*exp(cmplx(0,rp,kind=dpc))
             endif
             rw(i,alpha)=val
          end do
       end do
    elseif (d%wpack == 'old') then
       read(d%filewp)
       read(d%filewp)
       read(d%filewp)
       do
         read(d%filewp,iostat=icod)i,d%t0,rw
         if (icod /= 0) then
            exit
         endif
       end do
       inquire(unit=d%filewp,name=fnam)
       write(6,'(2a)') ' Letti dal file ',trim(fnam)
       write(6,'(a,i8,a,f15.6/)') ' step =',i,' time (fs)=',d%t0
       d%t0=d%t0*fstoau
    elseif (d%wpack == 'read') then
       do i=1,nrtot
         r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
         read(d%fileiniwp,*) rr,wpx
         rw(i,d%istati)=cmplx(wpx,0.0_dpr,kind=dpc)
       end do
       if (sum(abs(r-rr)) > eps) then
          write(6,'(a)') ' ERRORE IN INI_WP: rx DIVERSO DA r.'
          write(6,*)     ' r =',r
          write(6,*)     ' rx=',rr
          stop 12
       endif
    else
       write(6,'(3a)') ' wpack=',d%wpack,' non permesso'
       write(6,'(a)') ' ..... stop'
       stop 12
    endif
  end subroutine ini_wp
  !-------------------------------------------------------------
  subroutine hmol(d,hx)
    implicit none
    type (discr), intent (inout) :: d
    type (hamil), intent (inout) :: hx
    integer         :: i,j,k,kk,nstati,nstri
    integer         :: ndim,nrtot
    real (kind=dpr), dimension (d%ndim) :: r,rx
    real (kind=dpr), dimension ((d%nstati**2+d%nstati)/2,3) :: wd
    !
    nstati = d%nstati
    nstri  = (nstati**2+nstati)/2
    !
    ndim  = d%ndim
    nrtot = d%nrt(ndim)
    !
    if (d%filem < 1) then
       write(6,'(a)') ' Errore in lettura PES'
       stop 12
    else
       rewind d%filem
    endif
    !
    do i=1,nrtot
       r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
       !r=d%r_dim(i,:)
       if (d%hcomplex) then
          !read(d%filem,*) rx,(hx%h_c(i,k),k=1,nstri)
          read(d%filem,*) rx,(wd(k,1),wd(k,2),k=1,nstri)
          hx%h_c(i,1:nstri)=cmplx(wd(1:nstri,1),wd(1:nstri,2),kind=dpc)
       else
          read(d%filem,*) rx,(hx%h(i,k),k=1,nstri)
       endif
       if (sum(abs(r-rx)) > eps) then
          write(6,'(a)') ' ERRORE IN HMOL:'
          write(6,'(a)') ' rx DIVERSO DA r.'
          write(6,*)     ' r =',r
          write(6,*)     ' rx=',rx
          stop 12
       endif
    end do
    !
    if (d%radiation) then
       hx%edip = 0.0_dpr
       do i=1,nrtot
          r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
          read(d%lsr%filed,*) rx,(wd(k,1),k=1,nstri)
          read(d%lsr%filed,*) rx,(wd(k,2),k=1,nstri)
          read(d%lsr%filed,*) rx,(wd(k,3),k=1,nstri)
          do j=1,3
             hx%edip(i,:) = hx%edip(i,:) + d%lsr%e0(j)*wd(:,j)
          end do
       end do
    endif
    !
    if (d%hcomplex) then
       do k=1,nstati
          kk=(k*k+k)/2
          hx%h_c(:,kk) = hx%h_c(:,kk) + d%ene_add
       end do
    else
       do k=1,nstati
          kk=(k*k+k)/2
          hx%h(:,kk) = hx%h(:,kk) + d%ene_add
       end do
    endif
    !
    !
    call adiabt(d,hx)
    !
  end subroutine hmol
  !-------------------------------------------------------------
  subroutine adiabt(d,hx)
    implicit none
    type (discr), intent (in)    :: d
    type (hamil), intent (inout) :: hx
    real (kind=dpr), dimension (d%ndim)               :: r
    real (kind=dpr), dimension (d%nstati)             :: e
    real (kind=dpr), dimension (d%nstati)             :: abspot
    real (kind=dpr), dimension (3*d%nstati)           :: rw
    real (kind=dpr), dimension (d%nstati,d%nstati)    :: ah,u,ut
    complex (kind=dpc), dimension (d%nstati,d%nstati) :: ahc,uc,utc
    complex (kind=dpc), dimension (d%nstati,d%nstati) :: ute
    complex (kind=dpc), dimension (d%nstati*d%nstati) :: w
    complex (kind=dpc) :: vx
    real (kind=dpr)    :: vxr
    real (kind=dpr)    :: dh,dscr,e1,e2,num,den,ang,c,s,h11,h12,h22
    integer            :: i,j,k,kj,nstati,ndim,nrtot,alpha,lw,info
    !
    nstati = d%nstati
    !
    ndim  = d%ndim
    nrtot = d%nrt(ndim)
    vxr=d%time*0.5_dpr
    vx=cmplx(0,1,kind=dpc)*vxr
    lw=nstati**2
    !
    if (d%iwrt > 2) then
       write(6,*) '  R         Adiabatic energies'
    endif
    if (d%hcomplex) then
       do i=1,nrtot
          r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
          kj=0
          do j=1,nstati
             do k=1,j
                kj=kj+1
                ahc(j,k)=hx%h_c(i,kj)
                ahc(k,j)=conjg(hx%h_c(i,kj))
                !write(6,*) 'i,j,k,ahc(j,k)',i,j,k,ahc(j,k)
             end do
             abspot(j)=-aimag(ahc(j,j))
             ahc(j,j)=cmplx(real(ahc(j,j)),0,kind=dpc)
          end do
          call zheev('V','L',nstati,ahc,nstati,e,w,lw,rw,info)
          uc = conjg(transpose(ahc))
          hx%e(i,:) = e
          hx%dia_c(i,:,:) = uc
          !!write(6,"(a/,20(4f12.6/))") ' controllo: real(ahc)  =',((real(ahc(k,j)),j=1,nstati),k=1,nstati)
          !!write(6,"(a/,20(4f12.6/))") ' controllo: aimag(ahc) =',((aimag(ahc(k,j)),j=1,nstati),k=1,nstati)
          !
          utc = ahc
          do alpha=1,d%nstati
             ute(:,alpha)=utc(:,alpha)*exp(-vx*hx%e(i,alpha))
          end do
          hx%uteu(i,:,:)=matmul(ute,uc)
          !
          ! absorbing potential (diagonal imaginary term in the diabatic hamiltonian)
          if (maxval(abs(abspot)) > d%thr_abspot) then
             abspot=exp(vxr*abspot)
             do j=1,nstati
                hx%uteu(i,j,:)=hx%uteu(i,j,:)*abspot(:)
             end do
          endif
          !
          if (d%iwrt > 2) then
             write(6,"(20f15.8)") r,hx%e(i,:)
          endif
       end do
    else
       do i=1,nrtot
          r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
          kj=0
          do j=1,nstati
             do k=1,j
                kj=kj+1
                ah(j,k)=hx%h(i,kj)
                ah(k,j)=hx%h(i,kj)
             end do
          end do
          call eispd(nstati,nstati,ah,e,ah,u)
          u = transpose(ah)
          hx%e(i,:) = e
          hx%dia(i,:,:) = u
          !
          ut = ah
          do alpha=1,d%nstati
             ute(:,alpha)=ut(:,alpha)*exp(-vx*hx%e(i,alpha))
          end do
          hx%uteu(i,:,:)=matmul(ute,u)
          !
          if (d%iwrt > 2) then
             write(6,"(20f15.8)") r,hx%e(i,:)
          endif
       end do
    endif
  end subroutine adiabt
end module quant_sys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Struttura del file file_mol. Esempio: 3 stati      !
!                                                     !
!  r, h(1,1),h(2,1),h(2,2),h(3,1),h(3,2),h(3,3)       !
!                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
