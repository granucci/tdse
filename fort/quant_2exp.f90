module quant_2exp
  use quant_util, only: dpc,dpr,discr,hamil
  implicit none

  public :: f_2exp
  private

contains
  !-------------------------------------------------------------
  subroutine f_2exp(x,ix,hx,d)
    implicit none
    integer, intent (in)            :: ix
    real (kind=dpr), intent (in)    :: x
    type (hamil), intent (inout)    :: hx
    type (discr), intent (in)       :: d
    real (kind=dpr) :: re,rc,a2,alpha1,alpha2,beta,de,b,h12
    namelist /DAT2EXP/ re,rc,a2,alpha1,alpha2,beta,de,b,h12
    real (kind=dpr) :: x2,xc,xxc,hd11,hd12,hd22,dhd11,dhd12,dhd22
    real (kind=dpr) :: dx,e1,e2,tgt,ttgt,g12
    real (kind=dpr),save :: a1
    integer,save         :: ifirst = 0
    !
    if (ifirst == 0) then
       ifirst = 1
       read(5,nml=DAT2EXP)
       a1 = a2*exp((alpha1-alpha2)*rc) - de*exp(alpha1*rc)
    endif
    !
    x2=x+x
    xc=x-rc
    xxc=xc**2
    !
    hd11=a1*exp(-alpha1*x)+dE
    hd22=a2*exp(-alpha2*x)
    hd12=b*exp(-beta*xxc) + h12*sin(x)**2
    !
    dhd11=-a1*alpha1*exp(-alpha1*x)
    dhd22=-a2*alpha2*exp(-alpha2*x)
    dhd12=-2.0_dpr*b*beta*xc*exp(-beta*xxc) + h12*sin(x2)
    !
    dx=sqrt((hd22-hd11)**2 + 4.0_dpr*hd12**2)
    e1=(hd11+hd22 - dx)/2.0_dpr
    e2=(hd11+hd22 + dx)/2.0_dpr
    tgt=(e1-hd11)/hd12
    ttgt=tgt**2
    g12=(tgt*(dhd22-dhd11) +(1.0_dpr-ttgt)*dhd12)/(dx*(1.0_dpr+ttgt))
    !
    if (d%nonad_coupl <= 0) then
       hx%h(ix,1) = hd11
       hx%h(ix,2) = hd12
       hx%h(ix,3) = hd11
    else
       hx%h(ix,1) = e1
       hx%h(ix,2) = 0.0_dpr
       hx%h(ix,3) = e2
       hx%g(ix,1) = g12
       ! hx%t viene eventualmente calcolato in hmol
    endif
  end subroutine f_2exp
end module quant_2exp
