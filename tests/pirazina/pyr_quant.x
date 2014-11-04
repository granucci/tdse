#!/bin/csh
#
setenv rmin   -20.0
setenv rmax   20.0
setenv nr     100
setenv r0     0.0
setenv omega1 0.126
setenv omega2 0.074
setenv omega3 0.118
#
#
setenv botoan 0.529177249
setenv autokc 627.50955315
setenv autocm 219474.625
setenv autoev 27.2113957
setenv fstoau 41.341372
setenv zero   0.000000
#
#
set omega1=`echo "$omega1/$autoev"|bc -l`
set omega2=`echo "$omega2/$autoev"|bc -l`
set omega3=`echo "$omega3/$autoev"|bc -l`
set massa1=`echo "1.0/$omega1"|bc -l`
set massa2=`echo "1.0/$omega2"|bc -l`
set massa3=`echo "1.0/$omega3"|bc -l`
#
./pyr_pot <<EOF >out.pyr_pot
 &DAT
   xmin=$rmin,$rmin,$rmin,
   xmax=$rmax,$rmax,$rmax,
   nx=$nr,$nr,$nr,
 &END
EOF
#
#
/home/gio/tdse/fort/quant_ev <<EOF >out.quant
 &DAT
   ndim=3,
   rmin=$rmin,$rmin,$rmin,
   rmax=$rmax,$rmax,$rmax,
   nr=$nr,$nr,$nr,
   massa=$massa1,$massa2,$massa3,
   omega=$omega1,$omega2,$omega3,
   r0=0.0,0.0,0.0,
   p0=0.0,0.0,0.0,
   iwrt=1, nprt=20,
   nstati=2,  istati=2,  file_mol='out.pyr_pot',
   wpack='gauss', time=1.00, tcycles=100000, 
   adiabatize=1,
 &END
EOF
