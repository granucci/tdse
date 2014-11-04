#!/bin/csh
#
setenv rmin   0.2020
#setenv rmax   79.9925
#setenv nr     8399
setenv rmax   210.997500
setenv nr     22189
setenv r0     3.74297601369
setenv p0     0.0
setenv omega  0.00119375984
#setenv omega  0.00359943
#
#setenv mH  1837.2893342685
#setenv mLi 12650.846294100

setenv mH  29166.2162400000
setenv mLi 29166.2162400000

#
#
setenv botoan 0.529177249
setenv autokc 627.50955315
setenv autocm 219474.625
setenv fstoau 41.341372
setenv zero   0.000000
#
#
set mtot=`echo "$mLi+$mH"|bc -l`
set massa=`echo "$mLi*$mH/($mLi+$mH)"|bc -l`
#
set wcm=`echo "$omega*$autocm"|bc -l`
set ecin=`echo "($p0^2 + $massa*$omega/2)/(2*$massa)"|bc -l`
echo "omega (in cm-1) = $wcm        starting kinetic energy (a.u.) = $ecin"
/home/gio/tdse/fort/quant_ev <<EOF >out.quant_radiation
 &DAT
   rmin=$rmin,  rmax=$rmax,  nr=$nr,  massa=$massa,
   omega=$omega, r0=$r0, p0=$p0, 
   nstati=2,  istati=1,    file_mol='pot_dia',
   wpack='gauss', time=0.50, tcycles=58000, 
   nprt=1, nprt_wp=80, thrwp=1.d-6, adiabatize=2, iwrt=2,
   radiation=.true., E0=0.0005,0.0,0.0,
   lasertyp=2, omega0=0.04824,
   file_dip='dipd.dat',
   laserparm=4000.0, 1000.0,
 &END
EOF
