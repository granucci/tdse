program quant_ev
!-------------------------------------------------------------------------
! Risoluzione della TDSE su una griglia.
! Una o piu' dimensioni (r), uno o piu` stati.
!
! L'hamiltoniano viene dato in input in base diabatica e l'evoluzione
! e` effettuata col metodo Split-Operator.
! 
! Tutto in unita` atomiche (se non altrimenti specificato).
!
! Namelist &DAT 
! NDIM          = numero di coordinate
! RMIN          = min. valore di r (vettore di dimensione NDIM)
! RMAX          = max. valore di r (vettore di dimensione NDIM)
! NR            = numero di punti della griglia. La distanza dr fra
!                 un punto e l'altro della griglia e` dr=(RMAX-RMIN)/NR.
!                 Vettore di dimensione NDIM
! MASSA         = masse del sistema (vettore di dimensione NDIM)
! NSTATI        = numero di stati
! ISTATI        = stato di partenza (sul quale si trova il pacchetto 
!                 iniziale). Tale stato di partenza e` in
!                 rappresentazione diabatica se ADIABATIZE < 2,
!                 altrimenti e` in rappr. adiabatica.
! TIME          = time step
! TCYCLES       = numero totale di time steps
! FILE_WP       = nome del file (di output) non formattato contenente
!                 i risultati dell'evoluzione temporale.
!                 Default FILE_WP='WPTOT'.
! WPACK         = 'gauss' pacchetto iniziale gaussiano caratterizzato
!                         da: 
!                                     <r> = R0 
!                                     <p> = P0
!                           <r^2> - <r>^2 = 1/(2*MASSA*OMEGA)
!                  'old' pacchetto iniziale letto dalla fine del file
!                        FILE_WP: si continua un'evoluzione temporale 
!                        effettuata in precedenza. I nuovi risultati
!                        vengono appesi alla fine del file FILE_WP.
!                 'read' pacchetto iniziale letto dal file formattato 
!                        specificato in FILE_INIWP.
! OMEGA         = vedere WPACK.
! R0            = vedere WPACK.
! P0            = vedere WPACK.
! FILE_MOL      = nome del file (di input) dal quale vengono lette 
!                 energie e accoppiamenti, in rappresentazione diabatica. 
!                 Si veda il modulo "quant_sys" per la struttura di tale
!                 file.
! FILE_INIWP    = nome del file formattato (di input) dal quale viene letto 
!                 il pacchetto iniziale se WPACK='read'. FILE_INIWP e` 
!                 costituito da tanti record quanti sono i punti della 
!                 griglia; ogni record e` del tipo:
!                 r  wpx
!                 dove r(1:ndim) contiene il valore delle coordinate nel dato
!                 punto della griglia e "wpx" e` il valore (reale) del pacchetto
!                 d'onda nel punto considerato per lo stato ISTATI.
!                 Solo se WPACK='read'.
! HCOMPLEX      = F, l'hamiltoniano diabatico e` reale (default);
!                 T, l'hamiltoniano diabatico e` complesso.
! ADIABATIZE    = 0 cond. iniziali e risultati in rappr. diabatica (Default).
!                 1 condizioni iniziali in rappresentazione diabatica,
!                   risultati in rappresentazione adiabatica.
!                 2 condizioni iniziali e risultati sono dati in 
!                   rappresentazione adiabatica.
!                 3 condizioni iniziali in rappresentazione adiabatica,
!                   risultati in rappresentazione diabatica.
! NPRT          = numero di timestep fra una scrittura e l'altra.
!                 Default NPRT=1.
! NPRT_WP       = come NPRT, ma per la funzione d'onda.
!                 Default NPRT_WP=NPRT.
! ENE_ADD       = fattore da aggiungere alle energie. Default
!                 ENE_ADD=0. 
! IWRT          = 0 scrittura normale (default).
!               > 0 si scrive l'energia totale e cinetica.
!               > 1 si scrivono i risultati dell'evoluzione (ovvero
!                   il pacchetto d'onda) sul file FILE_WP. Se 
!                   NDIM=1 i pacchetti vengono anche scritti
!                   sul file formattato 'WP'.
!               > 2 si scrivono le energie adiabatiche (griglia).
!               > 3 si scrivono le pop. e le energie cinetiche degli
!                   stati diabatici (qualunque sia la rappr. nella quale
!                   vengono dati i risultati).
!               > 9 si scrive la griglia e si interrompe
!                   il calcolo (utile per NDIM > 1).
! THRWP         = soglia per la scrittura dei pacchetti sul file
!                 formattato 'WP'. Si scrivono solo i punti della
!                 griglia caratterizzati:
!                 a) da una popolazione maggiore di THRWP
!                 b) compresi fra due punti con una popolazione 
!                    maggiore di THRWP.
! RADIATION     = F, calcolo normale (default);
!                 T, interazione con la radiazione.
!---- keywords rilevanti solo per RADIATION=T -------------------------
! FILE_DIP      = nome del file (di input) dal quale vengono letti
!                 i dipoli. Ha la stessa struttura di FILE_MOL, eccetto 
!                 che per ogni punto della griglia vengono date tre
!                 schede per le tre componenti del dipolo.
! E0            = vettore tridimensionale che fornisce l'ampiezza 
!                 massima del campo elettrico.
! OMEGA0        = frequenza dell'onda portante. Per il momento la sola
!                 onda portante implementata e` del tipo cos(OMEGA0*t).
! LASERCW       = 0, onda portante tipo cos(OMEGA0*t), default;
!                 1, onda portante tipo sin(OMEGA0*t).
! LASERTYP      = intero che specifica la forma analitica dell'inviluppo.
!                 Il campo e` dato da  (inviluppo) * (onda portante).
!                 1, inviluppo = 1;
!                 2, inviluppo gaussiano = exp(-((t-t0)/tau)**2).
!                 3, inviluppo cossq     = cos(pi*(t-t0)/(2*tau))
!                    nell'intervallo t0-tau <= t <= t0+tau. Al difuori di tale
!                    intervallo l'inviluppo vale 0.
! LASERPARM     = parametri dell'inviluppo.
!                 - caso LASERTYP = 2: inviluppo gaussiano.
!                   t0  = LASERPARM(1)
!                   tau = fwhm/sqrt(2*log(2))
!                   fwhm = LASERPARM(2)
!                 - caso LASERTYP = 3: inviluppo cossq.
!                   t0  = LASERPARM(1)
!                   tau = LASERPARM(2)
!----------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------
  use quant_util, only: dpr,discr,mdim
  use quant_mod, only: quant_integ
  use quant_efield, only: pulse_ini
  implicit none
  integer               :: ndim
  integer               :: nstati,tcycles,nprt,istati,nprt_wp
  integer               :: iwrt,adiabatize,lasertyp,lasercw
  logical               :: hcomplex,radiation
  real(kind=dpr)        :: time,ene_add,thrwp,omega0
  character (len=120)   :: file_mol,file_wp,file_dip,file_iniwp
  character (len=10)    :: wpack
  real(kind=dpr), dimension(3)     :: e0
  real(kind=dpr), dimension(20)    :: laserparm
  integer, dimension (mdim)        :: nr
  real(kind=dpr), dimension (mdim) :: rmin,rmax,massa,omega,r0,p0
  namelist /DAT/ rmin,rmax,nr,massa,nstati,time,tcycles,file_wp,wpack, &
   &  omega,r0,p0,istati,file_mol,nprt,ene_add,thrwp,nprt_wp, &
   &  iwrt,adiabatize,ndim,hcomplex, &
   &  radiation,omega0,e0,lasertyp,lasercw,laserparm,file_dip,file_iniwp
  type (discr)          :: d
  integer               :: iw,i,nrtot
  !
  !
  ndim          = 1
  nprt          = 1
  nprt_wp       = -1
  wpack         = 'gauss'
  ene_add       = 0.0_dpr
  thrwp         = 0.0_dpr
  adiabatize    = 0
  file_wp       = 'WPTOT'
  iwrt          = 0
  hcomplex      = .false.
  radiation     = .false.
  lasercw       = 0
  file_dip      = ''
  read(5,DAT)
  if (nprt_wp == -1) nprt_wp = nprt
  write(6,DAT)
  if (mdim < ndim) then
     write(6,*) 'Stop: increase mdim'
     stop 12
  endif
  allocate (d%nr(ndim),d%r0(ndim),d%p0(ndim),d%massa(ndim),d%omega(ndim))
  allocate (d%dr(ndim),d%rmin(ndim),d%rmax(ndim),d%nrt(ndim))
  d%ndim        = ndim
  d%rmin        = rmin(1:ndim)
  d%rmax        = rmax(1:ndim)
  d%nr          = nr(1:ndim)
  d%nstati      = nstati
  d%istati      = istati
  d%massa       = massa(1:ndim)
  d%time        = time
  d%tcycles     = tcycles
  d%nprt        = nprt
  d%nprt_wp     = nprt_wp
  d%omega       = omega(1:ndim)
  d%r0          = r0(1:ndim)
  d%p0          = p0(1:ndim)
  d%wpack       = wpack
  d%ene_add     = ene_add
  d%thrwp       = thrwp
  d%adiabatize  = adiabatize
  d%iwrt        = iwrt
  d%hcomplex    = hcomplex
  d%radiation   = radiation
  if (d%radiation) then
     d%lsr%omega0 = omega0
     d%lsr%e0 = e0
     d%lsr%lasercw = lasercw
     d%lsr%lasertyp = lasertyp
     d%lsr%laserparm = laserparm
  endif
  !
  iw            = 6
  !
  !
  do i=1,ndim
    d%dr(i) = (rmax(i)-rmin(i))/nr(i)
    if (i == 1) then
       d%nrt(1)      = d%nr(1)+1
    else
       d%nrt(i) = d%nrt(i-1)*(d%nr(i)+1)
    endif
  end do
  !
  d%filem = 2
  open(d%filem,file=file_mol,form='formatted',status='unknown')
  !
  d%filewp = 3
  if (d%wpack == 'old') then
     open(d%filewp,file=file_wp,form='unformatted',status='old')
  else
     open(d%filewp,file=file_wp,form='unformatted',status='unknown')
  endif
  rewind d%filewp
  !
  if (d%wpack == 'read') then
     d%fileiniwp=8
     open(d%fileiniwp,file=file_iniwp,form='formatted',status='old')
     rewind d%fileiniwp
  endif
  !
  if (d%radiation) then
     d%lsr%filed = 4
     open(d%lsr%filed,file=file_dip,form='formatted',status='old')
     call pulse_ini(d)
  endif
  !
  nrtot=d%nrt(ndim)
  allocate (d%inddim(nrtot,ndim),d%fftkin(nrtot),d%fftkin2(nrtot))
  !
  call quant_integ(d,iw)
  !
  deallocate (d%nr,d%r0,d%p0,d%massa,d%omega)
  deallocate (d%dr,d%rmin,d%rmax,d%nrt)
  deallocate (d%inddim,d%fftkin,d%fftkin2)
end program quant_ev

