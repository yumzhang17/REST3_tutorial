#!/bin/bash

scale_intra=$1 # \lambda^pp
scale_km=$2    # \kappa_m
input=$3       # processed topology file e.g., p53.top
output=$4      # generated new topology

scale_inter=`echo $scale_intra $scale_km | awk '{print $2*sqrt($1)}'`
echo 'scale_km = ' $scale_km
echo 'scale_intra = ' $scale_intra
echo 'scale_inter = ' $scale_inter

bash partial_tempering.sh $scale_intra $scale_inter < $input  > topol-$scale_intra-$scale_inter.top

line1=`grep -n 'nonbond_params' topol-$scale_intra-$scale_inter.top | cut -d":" -f 1 | awk 'NR==1{print}'`
line2=`grep -n 'nonbond_params' topol-$scale_intra-$scale_inter.top | cut -d":" -f 1 | awk 'NR==2{print}'`
line1a=$(($line1+1))
line1a2=$(($line1+2))
line2b=$(($line2-1))
sed -n ''$line1a'h;'$line1a2','$line2b'H;$G;1,'$line1'p;'$line2',$p' topol-$scale_intra-$scale_inter.top > test.top
mv test.top $output
echo "done"
echo ""
