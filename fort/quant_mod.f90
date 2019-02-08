module quant_mod
  use quant_util, only: dpr,dpc,discr,hamil,fstoau,indset, &
                        &  aprod,ind_fft_kin
  use quant_sys, only: hmol,ini_wp
  use quant_efield, only: efield
  use fft_drv, only: fft_ini, fft_end, fft_fwrd, fft_bwrd
  implicit none
  private
  public :: quant_integ
  complex (kind=dpc), allocatable, dimension (:)   :: rwx
  integer, allocatable, dimension (:)              :: nrx
  !
contains
  !-------------------------------------------------------------
  subroutine quant_integ(d,iw)
    implicit none
    type (discr), intent (inout) :: d
    integer, intent (in)         :: iw
    !
    type (hamil) :: hel
    complex (kind=dpc), allocatable, dimension (:,:) :: rw
    integer, dimension (d%ndim) :: nrinv
    integer :: i,npt,nstati,nr1,nstri,nstru
    !
    nstati = d%nstati
    nstri  = (nstati*nstati+nstati)/2
    nstru  = (nstati*nstati-nstati)/2
    nr1    = d%nrt(d%ndim)
    npt    = nr1*nstri
    write(iw,*) 'Numero totale di punti nella griglia =',nr1
    write(iw,*) 'dr =',d%dr
    allocate (rw(nr1,nstati))
    allocate (hel%uteu(nr1,nstati,nstati))
    allocate (hel%e(nr1,nstati))
    if (d%hcomplex) then
       allocate (hel%h_c(nr1,nstri))
       allocate (hel%dia_c(nr1,nstati,nstati))
    else
       allocate (hel%h(nr1,nstri))
       allocate (hel%dia(nr1,nstati,nstati))
    endif
    if (d%radiation) then
       allocate (hel%edip(nr1,nstri))
    endif
    !
    call indset(d)
    call ind_fft_kin(d)
    call wrtgrd(d,iw)
    !
    allocate (rwx(nr1),nrx(d%ndim))
    !
    call fft_ini(d,rwx,nrx)
    !
    d%t0 = 0.0_dpr
    call ini_wp(d,rw)
    call hmol(d,hel)
    if (d%adiabatize > 1) then
       call rw_dia(d,rw,hel,"dia")
    endif
    !
    d%ttot = d%t0
    call quant_resu(d,rw,hel,0,iw)
    do i=1,d%tcycles
       d%ttot = d%ttot + d%time
       call split_op(d,hel,rw)
       if (mod(i,d%nprt)==0) then
          call quant_resu(d,rw,hel,i,iw)
       endif
    end do
    !
    call fft_end
    !
    deallocate (rwx,nrx)
    deallocate (rw)
    deallocate (hel%uteu)
    deallocate (hel%e)
    if (d%hcomplex) then
       deallocate (hel%h_c, hel%dia_c)
    else
       deallocate (hel%h, hel%dia)
    endif
  end subroutine quant_integ
  !-------------------------------------------------------------
  ! trasf == "dia"  => passaggio da adia a dia:  wp = wp * hel%dia
  ! trasf /= "dia"  => passaggio da dia a adia:  wp = wp * (hel%dia)+
  subroutine rw_dia(d,rw,hel,trasf)
    implicit none
    type (discr), intent (in)                           :: d
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    type (hamil), intent (in)                           :: hel
    character (len=*), intent (in)                      :: trasf
    complex (kind=dpc), dimension(d%nstati) :: rwx
    integer :: i,alpha
    !
    if (d%hcomplex) then
       if (trasf == "dia") then
          do i=1,d%nrt(d%ndim)
             rwx=rw(i,:)
             do alpha=1,d%nstati
                rw(i,alpha) = sum(rwx*hel%dia_c(i,:,alpha))
             end do
          end do
       else
          do i=1,d%nrt(d%ndim)
             rwx=rw(i,:)
             do alpha=1,d%nstati
                !rw(i,alpha) = sum(rwx*hel%dia_c(i,alpha,:))
                rw(i,alpha) = sum(rwx*conjg(hel%dia_c(i,alpha,:)))
             end do
          end do
       endif
    else
       if (trasf == "dia") then
          do i=1,d%nrt(d%ndim)
             rwx=rw(i,:)
             do alpha=1,d%nstati
                rw(i,alpha) = sum(rwx*cmplx(hel%dia(i,:,alpha),0,dpc))
             end do
          end do
       else
          do i=1,d%nrt(d%ndim)
             rwx=rw(i,:)
             do alpha=1,d%nstati
                rw(i,alpha) = sum(rwx*cmplx(hel%dia(i,alpha,:),0,dpc))
             end do
          end do
       endif
    endif
    !
  end subroutine rw_dia
  !-------------------------------------------------------------
  subroutine split_op(d,hel,rw)
    implicit none
    type (discr), intent (in)                           :: d
    type (hamil), intent (inout)                        :: hel
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    integer :: i,nrtot
    !
    nrtot=d%nrt(d%ndim)
    !
    if (d%radiation) call efield(d,hel)
    !
    do i=1,nrtot
      rw(i,:)=matmul(hel%uteu(i,:,:),rw(i,:))
    end do
    !
    call split_kin(d,rw)
    !
    do i=1,nrtot
      rw(i,:)=matmul(hel%uteu(i,:,:),rw(i,:))
    end do
    ! 
  end subroutine split_op
  !-------------------------------------------------------------
  subroutine split_kin(d,rw)
    implicit none
    type (discr), intent (in)                           :: d
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    integer :: alpha
    integer, save :: ifirst=0, ndim, nrtot
    !
    if (ifirst == 0) then
       ifirst = 1
       ndim=d%ndim
       nrtot=d%nrt(ndim)
    endif
    !
    do alpha=1,d%nstati
      rwx=rw(:,alpha)
      call fft_fwrd(rwx,nrx)
      rwx(1:nrtot)=rwx(1:nrtot)*d%fftkin(1:nrtot)
      !
      call fft_bwrd(rwx,nrx)
      rw(:,alpha)=rwx
    end do
  end subroutine split_kin
  !-------------------------------------------------------------
  subroutine wrtgrd(d,iw)
    implicit none
    type (discr), intent (in)  :: d
    integer, intent (in)       :: iw
    real (kind=dpr), dimension (d%ndim) :: r
    integer                             :: i,ndim,nrtot

    ndim=d%ndim
    nrtot=d%nrt(ndim)

    if (d%iwrt < 10) return

    write(iw,'(a)') 'n.     indexes        r values'
    do i=1,nrtot
      r(:) = d%rmin(1:ndim) + (d%inddim(i,1:ndim)-1)*d%dr(1:ndim)
      write(iw,'(i10,6i5)',advance='no')i,d%inddim(i,:)
      write(iw,'(6f12.6)')r
    end do

    if (d%iwrt >= 10) then
       stop 10
    endif
  end subroutine wrtgrd
  !-------------------------------------------------------------
  subroutine quant_resu(d,rw,hel,istep,iw)
    implicit none
    type (discr), intent (in)                           :: d
    type (hamil), intent (in)                           :: hel
    complex (kind=dpc), dimension (:,:), intent (inout) :: rw
    integer, intent (in)                                :: istep,iw
    complex (kind=dpc), parameter :: zzero  = (0.0_dpr,0.0_dpr)
    real (kind=dpr), parameter    :: zero   = 0.0_dpr
    integer                :: alpha,i,i1,ab,ip,im,beta,abg,iini,iend
    integer, save          :: ifirst=0, iwpf, ndim, nrtot
    real (kind=dpr), save  :: drtot
    complex (kind=dpc), dimension(size(rw,1)) :: bra,ket
    real (kind=dpr), dimension(d%nstati)      :: pop, ek, pop_dia
    real (kind=dpr), dimension((d%nstati*d%nstati+d%nstati)/2)  :: epx
    real (kind=dpr), allocatable, dimension(:,:)  :: e_pes
    real (kind=dpr)        :: ttotfs,epot,ekin,etot,pp,rx
    logical, save          :: adiab
    logical                :: zwrt_e
    !
    if (ifirst == 0) then
       ifirst = 1
       ndim=d%ndim
       nrtot=d%nrt(ndim)
       drtot=aprod(1,ndim,d%dr)
       adiab=(d%adiabatize == 1).or.(d%adiabatize == 2)
       !
       if (d%wpack /= 'old') then
          allocate(e_pes(nrtot,d%nstati))
          if (adiab) then
            e_pes=hel%e
          else
            if (d%hcomplex) then
              do alpha=1,d%nstati
                ab=(alpha*alpha+alpha)/2
                e_pes(1:nrtot,alpha)=hel%h_c(1:nrtot,ab)
              end do
            else
              do alpha=1,d%nstati
                ab=(alpha*alpha+alpha)/2
                e_pes(1:nrtot,alpha)=hel%h(1:nrtot,ab)
              end do
            endif
          endif
          write(d%filewp)d%ndim,d%nstati,d%istati,d%hcomplex,d%adiabatize
          write(d%filewp)d%rmin,d%dr,d%nrt,d%massa,d%omega,d%r0,d%p0
          write(d%filewp)d%inddim,e_pes
          deallocate (e_pes)
       endif
       !       
       if (d%iwrt > 0 .and. d%iwrt <= 3) then
          write(iw,'(a,2(15x,a),20x,2(15x,a))') &
               & '  Step', 'Time (fs)','pop','Etot','Ekin (diabatic)'
       elseif (d%iwrt > 3) then
          write(iw,'(a,2(15x,a),20x,4(15x,a))') &
               & '  Step', 'Time (fs)','pop','Etot','Ekin (diabatic)', &
               & '  pop (diabatic)        Ekin on diabats', &
               & '  Epot on diabats, triangular'
       else
          write(iw,'(a,2(15x,a))') &
               & '  Step', 'Time (fs)','pop'
       endif
       iwpf=10
       open(iwpf,file='WP',form='formatted',status='unknown')
       rewind iwpf
    endif
    ttotfs=d%ttot/fstoau
    !
    zwrt_e = (d%iwrt > 0) .and. (istep > 0)
    !
    ! --- calcolo energia totale e cinetica
    if (zwrt_e) then
       ek=zero
       epot=zero
       ab=0
       do alpha=1,d%nstati
          rwx=rw(:,alpha)
          call fft_fwrd(rwx,nrx)
          ek(alpha)=sum(d%fftkin2(:)*abs(rwx(:))**2)*drtot
          bra=conjg(rw(:,alpha))
          do beta=1,alpha
             ab=ab+1
             ket=rw(:,beta)
             if (d%hcomplex) then
                epx(ab)=real(sum(bra(:)*hel%h_c(:,ab)*ket(:)),kind=dpr)*drtot
             else
                epx(ab)=real(sum(bra(:)*hel%h(:,ab)*ket(:)),kind=dpr)*drtot
             endif
             if (beta /= alpha) then
               epx(ab)=epx(ab)*2.0_dpr
             endif
             epot=epot+epx(ab)
          end do
       end do
       ekin=sum(ek)
       etot=ekin+epot
       if (d%iwrt > 3) then
          do alpha=1,d%nstati
             pop_dia(alpha)=sum(abs(rw(:,alpha))**2)*drtot
          end do
       endif
    endif
    !
    ! --- se necessario si passa alla rappr. adiabatica
    if (adiab) then
       call rw_dia(d,rw,hel,"adia")
    endif
    !
    ! --- calcolo popolazione degli stati
    do alpha=1,d%nstati
       !pop(alpha)=dot_product(bra,ket)*drtot
       pop(alpha)=sum(abs(rw(:,alpha))**2)*drtot
    end do
    !
    ! --- scrittura risultati
    if (.not.zwrt_e) then
       write(iw,'(i10,50f20.10)') istep,ttotfs,pop
    else
       if (d%iwrt > 3) then
          write(iw,'(i10,50f20.10)') &
            & istep,ttotfs,pop,etot,ekin,pop_dia,ek,epx
       else
          write(iw,'(i10,50f20.10)') istep,ttotfs,pop,etot,ekin
       endif
    endif
    !
    ! --- eventuale scrittura pacchetto d'onda sul file WP
    if (mod(istep,d%nprt_wp)==0) then
       if (d%iwrt > 1 .and. d%ndim == 1) then
          do iini=1,d%nrt(1)
            pp=sum(abs(rw(iini,:))**2)
            if (pp > d%thrwp) then
              exit
            endif
          end do
          do iend=d%nrt(1),1,-1
            pp=sum(abs(rw(iend,:))**2)
            if (pp > d%thrwp) then
              exit
            endif
          end do
      
          do i=iini,iend
             rx  = d%rmin(1) + (i-1)*d%dr(1)
             write(iwpf,'(41f20.10)')rx,(rw(i,alpha),alpha=1,d%nstati)
          end do
          write(iwpf,*)
       endif
       !
       if (d%iwrt > 1) then
          write(d%filewp)istep,ttotfs,rw
       endif
    endif
    !
    ! --- se necessario si ritorna alla base diabatica
    if (adiab) then
       call rw_dia(d,rw,hel,"dia")
    endif
    !
  end subroutine quant_resu
end module quant_mod
