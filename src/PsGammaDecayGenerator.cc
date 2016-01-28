// 
// The following class PsDecayGenerator generates momenta vectors for decays of oPs.
// It handles the following types of decays:
//
// 1. o-Ps Angular distribution for CP-violating decays.
// 2. o-Ps decays using Ore-Powell Energy distribution and Bernreuther angular distributions, as predicted by QED
//
// All units are in terms of electron rest masses., ie. energy=1 corresponds to one electron rest mass.
// The spin/magnetic field axis is in the +z direction in the coordinate system of the returned vectors.
// 
// Author: R. Henning
// Date: 3/5/2015
// 
#include "Randomize.hh"
#include <iostream>
#include <TF2.h>
#include <TF1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TROOT.h>
#include "G4PhysicalConstants.hh"
#include "globals.hh"
#include "PsGammaDecayGenerator.hh"
using namespace CLHEP;
// ---------------------------------------------------------------------------------- //
PsGammaDecayGenerator::PsGammaDecayGenerator()
{
}

PsGammaDecayGenerator::~PsGammaDecayGenerator()
{
}
Double_t PsGammaDecayGenerator::oPsAngularDistribution_m_1(Double_t *theta, Double_t *p)
{
// Distribution of angle between normal of decay plane and spin axis for m=+-1. 

	Double_t c = TMath::Cos(*theta); 
	return 0.5*(3-c*c);
}

// ---------------------------------------------------------------------------------- //
 		
Double_t PsGammaDecayGenerator::oPsAngularDistribution_m_0(Double_t *theta, Double_t *p)
{
// Distribution of angle between normal of decay plane and spin axis for m=0. 

	Double_t c = TMath::Cos(*theta); 
	return 1.0+c*c;
}
		
// ---------------------------------------------------------------------------------- //
 
Double_t PsGammaDecayGenerator::oPsEnergyDistribution_OP(Double_t *k, Double_t *p)
{
// Energy distribution of single gamma-rays emitted during o-Ps decay assuming SM physics
// k is units of electron masses.
// From Ore, Powell, Phys. Rev. 75, 1696 (1949)

	Double_t a = 1 - *k;
	Double_t b = 2 - *k;
	return 2*(*k*a/(b*b) - 2*a*a/(b*b)*TMath::Log(a) + b/(*k) + 2*a/((*k)*(*k))*TMath::Log(a)); 
}

// ---------------------------------------------------------------------------------- //

Double_t PsGammaDecayGenerator::oPsEnergyDistribution_PS(Double_t *k, Double_t *p)
{
// Energy distribution of single gamma-rays emitted during o-Ps decay from phase space only
// k is units of electron masses.
// From Ore, Powell, Phys. Rev. 75, 1696 (1949)

	Double_t kk = (*k)*(*k);
	return kk*(6-6*(*k)+kk); 
}

// ---------------------------------------------------------------------------------- //
 
void PsGammaDecayGenerator::Initialize()
{

	foPsAngularDistribution_m_1 = new TF1("oPsAngularDistribution_m_1", this, &PsGammaDecayGenerator::oPsAngularDistribution_m_1, 0, TMath::Pi(), 0);
	foPsAngularDistribution_m_0 = new TF1("oPsAngularDistribution_m_0", this, &PsGammaDecayGenerator::oPsAngularDistribution_m_0, 0, TMath::Pi(), 0);
	foPsAngularDistribution_CP_theta_phi = new TF2("oPsAngularDistribution_CP_theta_phi", "1+sin(2*y)*cos(x)", -TMath::Pi(), TMath::Pi(), 0, TMath::Pi());
	foPsEnergyDistribution_OP = new TF1("oPsEnergyDistribution_OP", this, &PsGammaDecayGenerator::oPsEnergyDistribution_OP, 0.001, 0.999, 0); //Doesn't like 0,1
	foPsEnergyDistribution_PS = new TF1("oPsEnergyDistribution_PS", this, &PsGammaDecayGenerator::oPsEnergyDistribution_PS, 0., 1., 0);

	foPsAngularDistribution_CP_theta_phi->SetNpx(200);
	foPsAngularDistribution_CP_theta_phi->SetNpy(100);
	
	delete gRandom;
	// static interface
	const G4long* table_entry;
	G4long theSeed;

	//theSeed = CLHEP::HepJamesRandom::getSeed();
	theSeed = CLHEP::HepRandom::getTheSeed();
	//CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	const long* seeds;
	//	CLHEP::HepRandom::getTheTableSeeds(seeds,theSeed);
	seeds = CLHEP::HepRandom::getTheSeeds();
	//	G4cout << "SEED: " << theSeed << " " <<  seeds[theSeed] << G4endl;
	gRandom = new TRandom3(theSeed); // Use UUID as seed. 
	/*	if(table_entry!=NULL)
	  {
	    gRandom = new TRandom3(*table_entry); // Use UUID as seed. 
	  }
	else
	  {
	    gRandom = new TRandom3(0); // Use UUID as seed. 
	    }*/

}

// ---------------------------------------------------------------------------------- //

Int_t PsGammaDecayGenerator::GenerateoPs(TVector3 *k1, TVector3 *k2, TVector3 *k3, Int_t m)
{
    
	DecayPlaneKinematics(k1, k2, k3, foPsEnergyDistribution_OP->GetRandom(),
                         foPsEnergyDistribution_OP->GetRandom());
	
	Double_t angle = (m == 0) ? foPsAngularDistribution_m_0->GetRandom() : foPsAngularDistribution_m_1->GetRandom();
	//Double_t angle=gRandom->Rndm()*TMath::Pi();
	k1->RotateX(angle);
	k2->RotateX(angle);
	k3->RotateX(angle);
	
	angle = gRandom->Rndm()*TMath::Pi()*2;
	k1->RotateZ(angle);
	k2->RotateZ(angle);
	k3->RotateZ(angle);
	
	return 1;
}


// ---------------------------------------------------------------------------------- //
// ---------------------------------------------------------------------------------- //

Int_t PsGammaDecayGenerator::GenerateoPs_CP(TVector3 *k1, TVector3 *k2, TVector3 *k3)
{
    
	// The energy distribution of the photons is model-dependent, hence we don't know what
	// it will be. For now we just pick them using phase
	// space calculations from Ore-Powell (ie. the matrix element is a constant).
	//
	
	DecayPlaneKinematics(k1, k2, k3, foPsEnergyDistribution_PS->GetRandom(),
                         foPsEnergyDistribution_PS->GetRandom());

	TVector3 zAxis(0,0,1);
	
	// Find theta and phi. Psi is determined by kinematics. 
	Double_t theta, phi;
	foPsAngularDistribution_CP_theta_phi->GetRandom2(phi, theta);
	
	// Find rotation axis to rotate decay plane normal vector by. 
	TVector3 rotationAxis = k1->Cross(zAxis);
	
	// Rotate vectors on decay plane to proper phi value.
	k1->RotateZ(phi);
	k2->RotateZ(phi);
	k3->RotateZ(phi);
	
	// Rotate entire decay plane. Now the angle between the projection of the z-axis on the decay plane
	// and k1 is phi. 
	// k1 x k2 at this point may be pointing in the -z direction. Correct for this. 
	
	if(k1->Cross(*k2) * zAxis < 0) 
		theta = TMath::Pi() - theta;
	
	// Check that the rotation axis is not zero
	if((zAxis.X()!=0)&&(zAxis.Y()!=0)&&(zAxis.Z()!=0)){
	k1->Rotate(theta, rotationAxis);
	k2->Rotate(theta, rotationAxis);
	k3->Rotate(theta, rotationAxis);
	}
	
	return 1;
}


// ---------------------------------------------------------------------------------- //

void PsGammaDecayGenerator::DecayPlaneKinematics(TVector3 *k1, TVector3 *k2, TVector3 *k3, Double_t x1, Double_t x2)
{
// Generates 3 random momentum vectors on the z-plane. The vectors are kinetically constrained so that they sum to zero.
// Two of the vectors will have magnitude x1 and x2, with the third being 2 - x1 - x2
// k1 and k2 will have the largest and second largest momenta, respectively. 
	
	Double_t k1_mag, k2_mag, k3_mag;
	Double_t x3 = 2.0 - x1 - x2;
	
// Sort so that k1_mag > k2_mag > k3_mag
	
	if(x1>x2)
		if(x3>x1) {
			k1_mag = x3;
			k2_mag = x1;
			k3_mag = x2;
		} else if(x2 > x3) {
			k1_mag = x1;
			k2_mag = x2;
			k3_mag = x3;
		} else {
			k1_mag = x1;
			k2_mag = x3;
			k3_mag = x2;
			}
	else if(x3 > x2) {
			k1_mag = x3;
			k2_mag = x2;
			k3_mag = x1;
		} else if(x1 > x3) {
			k1_mag = x2;
			k2_mag = x1;
			k3_mag = x3;
		} else {
			k1_mag = x2;
			k2_mag = x3;
			k3_mag = x1;
			}

// Create three vectors on z-plane using kinematic constraints. k3 points along the -x axis for simplification.
// Since we have ranked the vectors, the equations below will generate the directions with a certain handedness, ie. k1, k2 and k3
// will always be in the same rotational order when viewed from, say the +z axis. We break this bias by 
// randomizing the sign of k1_y.
    
	Double_t k1_x = -0.5*((k2_mag*k2_mag - k1_mag*k1_mag)/k3_mag - k3_mag);
	Double_t k1_y = TMath::Sqrt(k1_mag*k1_mag - k1_x*k1_x) * ((gRandom->Rndm() > 0.5) ? -1 : 1);
    
	k1->SetXYZ(k1_x, k1_y, 0.);
	k2->SetXYZ(k3_mag - k1_x, -k1_y, 0.);
	k3->SetXYZ(-k3_mag, 0., 0.);

	
// Rotate by random angle around z-axis. k3 no longer points in -x direction. 
	Double_t angle = gRandom->Rndm()*TMath::Pi()*2;
	k1->RotateZ(angle);
	k2->RotateZ(angle);
	k3->RotateZ(angle);
}

// ---------------------------------------------------------------------------------- //
// ---------------------------------------------------------------------------------- //
		
		
