#!/bin/bash

disps='001'

for no in $disps
do

   echo $no
   mkdir disp-$no
   cp dftb_in.hsd disp-$no/
   cp geo.genS-$no disp-$no/geo.gen

   cd disp-$no
   dftb+ > out 2>&1
   cd ..  
done
