#!/bin/sh

for VAR in {0..24}
do
#echo '/^"$VAR"/p'
sed -n '/^'"${VAR}\s"'/p' energies.dat >> outfile_$VAR
awk '{print $NF}' outfile_$VAR >> crystalEnergy_$VAR
#rm outfile_$VAR
done
