#CALIOPE

##Geant Simulation Code for CALIOPE
The point of this simulation is to look at the decay of orthopositronium in the APEX array located at TUNL at Duke University. We seek to determine whether or not our experimental configuration will be achieve sensitivites necessary to be competitive with existing searches for CP violation in the charged lepton sector.

###Command Line Arguments
THIS IS NOT WORKING!!
The user can set the spin of the orthopositronium.
The spin will determine the angular distribution of the normal to the decay plane with respect to spin quantization axis (z-axis).
There are three choices: m=-1, m=0, and m=+1.

###Checks and Directory Structure
You can test the output of several pieces of the code. These are set up to run easily and quickly for testing purposes and/or sanity checks. They are as follows:
1) Theta Distribution: This is contained in the directory ThetaDist. You can run the program 'theta'. This generates the ROOT output file 'fillrandom.root', which can be examined by typing 'root phiHist.C'.
2) Phi Distribution: This is contained in the directory PhiDist. You can run the program 'theta'. This generates the ROOT output file 'fillrandom.root', which can be examined by typing 'root thetaHist.C'.
3) One can also check the theta distribution of a run using the directory 'ThetaDistCheck'. This will grab the initial momenta of all particles and check to make sure that geant outputs the expected theta distribution. If everything is working fine, the histogram should look identical to the input distribution in the 'ThetaDist' directory.

###Design

How it works:

The main CALIOPE code starts up the user interface. 
It initializes detector geometry.
It initializes physics list.

Implementation of theta distribution:
Coordinate axes drawn to see be able to better visualize theta distribution
User decides the spin of the decaying positronium
If you run it with the macro 'gps.mac', you can look at two back-to-back gamma rays.

RootOutputData.hh contains all variables in the TTree which can be used to plot relevant data.
The macro writes a ROOT file 511keVgamma.root.

The source code includes the following key components:
DetectorConstruction.cc contains contains the code to build the APEX array.

####Processing

A typical run has, on average, more than a million events. The code is set up to be run in parallel.

The CALIOPE.cc code directs all output to ROOT files. This location is in the main() function of CALIOPE.cc. Right now, it goes to the scratch space on the Killdevil and Kure clusters. This should be adjusted to specify locations for other computing clusters.

PRINTING ENERGY SPECTRA OF INDIVIDUAL CRYSTALS

The code is designed to be run on the UNC university's clusters (though it can be run on any other cluster). The script 'cleanup.sh' is first run to clean out the appropriate folder ('scratch space') on Killdevil and/or Kure. Then, 'batch.sh' is run, producing a multitude of files in the scratch space. 'Batch.sh' uses the linux 'bsub' command from LSF. The script 'merge.sh' merges the files in the scratch space into one giant file. This is the 'finale' root file which we will process.

The script 'processRoot.sh' grabs the hit energies for each crystal 'copy number' and prints them to a file. It sums all the hit energies for each crystal in a single event. This file should then be copied to one of the sub-directories of the Analysis directory for further processing.

To make the energy spectra, first go to the relevant directory under the Analysis directory. Run the getCrystalEnergies.sh script. This should output several files of the name crystalEnergy_* and outputfile_*. Copy these files to the crystalEnergyFiles/crystalEnergyData directory. Then, back in the main directory, run the 'loopFiles.sh' script. This will generate the corresponding histograms for the energy spectra in the folder crystalEnergyFiles/crystalEnergyHists. 

PRINTING THE TOTAL ENERGY SPECTRA

The total energy spectra can likewise be printed using BLAH BLAH BLAH.

LOOKING FOR COINCIDENT EVENTS
