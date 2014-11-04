#!/bin/csh
#
setenv xmin   0.20
setenv xmax   10.0
setenv ymin   -3.0
setenv ymax   3.0
setenv nrx    1000
setenv nry    200
setenv mx     20000
setenv my     6666.6666667
setenv kx     0.02 
setenv ky     0.10 
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
./con_int <<EOF >out.con_pot
 &DAT
   xmin=$xmin,$ymin,
   xmax=$xmax,$ymax,
   nx=$nrx,$nry,
   mx=$mx,  my=$my,
   x1=4.0,  x2=3.0,  x3=3.0,  
   kx=$kx,
   ky=$ky,
   delta=0.01,
   alpha=3.0,
   beta=1.5,
   gamma=0.08,
 &END
EOF
#
set omegax=`echo "sqrt($ky/$mx)"|bc -l`
set omegay=`echo "sqrt($ky/$my)"|bc -l`
#
set Dx=`echo "1/sqrt($mx*$omegax)"|bc -l`
set Dy=`echo "1/sqrt($my*$omegay)"|bc -l`
echo "Dx=" $Dx
echo "Dy=" $Dy
#
/home/gloria/tdse/fort/quant_ev <<EOF >out.quant
 &DAT
   ndim=2,
   rmin=$xmin,$ymin,
   rmax=$xmax,$ymax,
   nr=$nrx,$nry,
   massa=$mx,$my,
   omega=$omegax,$omegay,
   r0=2.0,0.0,
   p0=0.0,0.0,
   iwrt=0,
   nstati=2,  istati=1,  file_mol='out.con_pot',
   wpack='gauss', time=0.50, tcycles=10000, 
   adiabatize=1,
 &END
EOF
#
rm out.con_pot
