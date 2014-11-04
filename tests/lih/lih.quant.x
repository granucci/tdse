#!/bin/csh
#
setenv rmin   0.2020
setenv rmax   210.997500
setenv nr     22189
setenv r0     12.0
setenv omega  0.008
setenv vli    0.0
setenv vh     0.0
#
setenv mH  1837.2893342685
setenv mLi 12650.846294100
#
#
setenv botoan 0.529177249
setenv autokc 627.50955315
setenv autocm 219474.625
setenv fstoau 41.341372
setenv zero   0.000000
setenv epsg   0.05     
#
#
set mtot=`echo "$mLi+$mH"|bc -l`
set massa=`echo "$mLi*$mH/($mLi+$mH)"|bc -l`
set p0=`echo "$massa*($vh - $vli)"|bc -l`
goto quant
#
set cvh=0.001
set cvli=`echo "-$cvh*$mH/$mLi"|bc -l`
set dcv=`echo "$cvh-($cvli)"|bc -l`
#
set dr=`echo "($rmax-$rmin)/$nr"|bc -l`
set nr1=`echo "$nr + 1"|bc -l`
set k=0
rm mol_dat
la:
  set k=`echo "$k + 1"|bc -l`
  if ($k > $nr1) goto fine
  set r=`echo "$rmin+($k-1)*$dr"|bc -l`
  set rang=`echo "$r*$botoan"|bc -l`
  sed -e "s/DDDDDD/$rang/" mop_template.dat |\
  sed -e "s/EEEEEE/$cvh/"|\
  sed -e "s/FFFFFF/$cvli/" >mop_inp.dat
  job.x mop_inp
  set e1=`rdyn mop_inp.dyn 17| awk '{if ($4 == 0 && $9 == 1) print $10}'`
  set e2=`rdyn mop_inp.dyn 18| awk '{if ($4 == 0 && $9 == 1) print $10}'`
  set g=`rdyn mop_inp.dyn 32| awk '{if ($4 == 1 && $9 == 1) print $10}'`
  set dt=`rdyn mop_inp.dyn 1| awk '{if ($4 == 1 && $9 == 1) print $10}'`
  set g=`echo "$g/($dcv*$dt*$fstoau)"|bc -l`
  if ($k > 1) then
    set dg=`echo "scale=0; sqrt(($g-($g0))^2)/$epsg"|bc -l`
    if ($dg > 1) set g=`echo "-1*($g)"|bc -l`
  endif
  echo "$r    $e1    $zero    $e2    $g" >>mol_dat
  set g0=$g
  goto la
fine:
#
#
#
quant:
set wcm=`echo "$omega*$autocm"|bc -l`
set ecin=`echo "($p0^2 + $massa*$omega/2)/(2*$massa)"|bc -l`
echo "omega (in cm-1) = $wcm        energia cinetica iniziale (a.u.) = $ecin"
/home/gloria/tdse/fort/quant_ev <<EOF >out.quant
 &DAT
   ndim=1,
   rmin=$rmin,  rmax=$rmax,  nr=$nr,  massa=$massa,
   omega=$omega, r0=$r0, p0=$p0,  ene_add=0.2131,
   nstati=2,  istati=2,   file_mol='pot_dia',
   wpack='gauss', time=0.50, tcycles=800, 
   nprt=80, thrwp=1.d-6, adiabatize=2,
 &END
EOF
