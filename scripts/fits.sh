#!/bin/bash

#data=data/wwnocuts_v2.lhe
data=data/jul17.lhe

#########################################################################################

# Two parameter fits
if false; then
for f in tilt ediff
do
  # all
  ./bin/fit_polariz -i $data -f $f -n 20 -o fit_"$f"_all.pdf
  for e in lab cm
  do
    # binned
    ./bin/fit_polariz -i $data -f $f -n 12 -o fit_"$f"_bins_energy_"$e".pdf --bins-energy-$e 100:1000:25
  done
done
fi

# One parameter fits
if false; then
for f in ediff_abs
do
  # all
  ./bin/fit_polariz -i $data -f $f -n 20 --draw-chisq -o fit_"$f"_all.pdf
  for e in lab cm
  do
    # binned
    ./bin/fit_polariz -i $data -f $f -n 12 --draw-chisq -o fit_"$f"_bins_energy_"$e".pdf --bins-energy-$e 100:1000:25
  done
done
fi

#########################################################################################

f=op

# One parameter fits
if true; then
for e in lab # cm
do
  for n in 10 15 50; do
    ./bin/fit_polariz -i $data -f $f -n $n --draw-chisq -o fit_"$f"_bins_energy_"$e"_n$n.pdf --bins-energy-$e 100:1000:25 -d weighted/opang/energy_$e/nocut2
  done
done
fi

#########################################################################################

# Fits withing |cos(production angle)| bins
if false; then
./bin/fit_polariz -i $data -f tilt  --bins-cos-prod 0:1:0.1 -o fit_tilt_bins_cos_prod.pdf
./bin/fit_polariz -i $data -f ediff --bins-cos-prod 0:1:0.1 -o fit_ediff_bins_cos_prod.pdf -n 50
./bin/fit_polariz -i $data -f ediff_abs --bins-cos-prod 0:1:0.1 -o fit_ediff_abs_bins_cos_prod.pdf -n 50
fi

# Opening angle fits
if false; then
./bin/fit_polariz -i $data -f op --bins-energy-lab 100:1000:25 -o fit_op_bins_energy_lab.pdf -d weighted/opang/energy_lab/nocut2 --draw-chisq
./bin/fit_polariz -i $data -f op --bins-energy-cm  100:1000:25 -o fit_op_bins_energy_cm.pdf  -d weighted/opang/energy_cm/nocut

./bin/fit_polariz -i $data -f op_norm --bins-energy 100:1000:25 -o fit_op_norm_bins_energy.pdf -d weighted/opang_norm/energy/nocut -n150
fi


