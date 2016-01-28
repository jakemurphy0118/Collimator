#! /bin/bash
myPath='/nas02/home/c/b/cbartram/CALIOPE-AWESOME/Analysis/Energy_Per_Crystal_2photon/crystalEnergyFiles/crystalEnergyData/crystalEnergy_'
myPath2='/nas02/home/c/b/cbartram/CALIOPE-AWESOME/Analysis/Energy_Per_Crystal_2photon/crystalEnergyFiles/crystalEnergyHists/crystalEnergyHist_'
for VAR in {0..23}
do
INVAR=$VAR
INVAR=$myPath$INVAR
OUTVAR=$myPath2$VAR
root -b -q './energyHist.C("'$INVAR'","'$OUTVAR'",'$VAR')'
done
