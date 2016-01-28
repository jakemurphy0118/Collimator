#include <iostream>
#include "TimeDistribution.hh"
#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <cstring>
#include <math.h>
TimeDistribution::TimeDistribution(){
  //  hist_1 = new TH1F("theta_hist","Angular Distribution of Normal to Decay Plane",num_bins,x_low,x_high);
}

TimeDistribution::~TimeDistribution(){
}

//Distribution for m=0
Double_t TimeDistribution::pPSDecay(Double_t *x, Double_t *par)
{
  Double_t xx=x[0];
  Double_t halfLife = par[0];//0.0866;
  f=(1/halfLife)*exp(-(xx/halfLife));
  return f;  
}

Double_t TimeDistribution::oPSDecay(Double_t *x, Double_t *par)
{
  Double_t xx=x[0];
  Double_t halfLife = par[0];//98.4;
  f=(1/halfLife)*exp(-(xx/halfLife));
  return f;  
}

Double_t TimeDistribution::GetVal(){
  val=hist_1->GetRandom();
  return val;
}

Double_t TimeDistribution::GetBranchingRatio(std::string Distr, double H){
  // Number of 2 gamma rays as compared to 3
  
  double mu = 9.27e-21;
  double g = 2;
  double deltaE=1.28e-15;
  double y = (mu*g*H)/(deltaE/2);
  double gamma_t=(1/(142e-9));
  double gamma_s=(1/(.125e-9));
  double B=1;
  if(Distr.compare("pPSDecay")==0){
    B=(((2+y*y)-2*(sqrt(1+y*y)))/(y*y))*(gamma_s/gamma_t);
  }

  if(Distr.compare("oPSDecay")==0){
    B=(((2+y*y)-2*(sqrt(1+y*y)))/(y*y))*(gamma_s/gamma_t);
  }

  return B;
}

TH1F* TimeDistribution::buildHist(std::string Distr,double H)
{
  //Calculate the half-life of the perturbed triplet state and perturbed singlet state given a magnetic field in CGS
  double mu = 9.27e-21;
  double g = 2;
  double deltaE=1.28e-15;
  double y = (mu*g*H)/(deltaE/2);
  double gamma_t=(1/(142e-9));
  double gamma_s=(1/(.125e-9));
  double gamma=((gamma_t+gamma_s)/2)+(1/sqrt(1+y*y))*((gamma_t-gamma_s)/2);
  double halfLife_10=(1/gamma)*(1e9);
  gamma=((gamma_t-gamma_s)/2)+(1/sqrt(1+y*y))*((gamma_t+gamma_s)/2);
  double halfLife_00=(1/gamma)*(1e9);

  TRandom3 r(0);
  Double_t *x;
  Double_t *par;
  TimeDistribution * fptr = new TimeDistribution(); //create the user function calls
  
  TF1 *f1 = new TF1("f1",fptr,&TimeDistribution::pPSDecay,0,1,1,"TimeDistribution","myfunction");
  f1->SetParameters(0.125,0);

  TF1 *f2 = new TF1("f2",fptr,&TimeDistribution::oPSDecay,0,1,1,"TimeDistribution","myfunction");
  f2->SetParameters(142.,0);

  TF1 *f3 = new TF1("f3",fptr,&TimeDistribution::oPSDecay,0,1,1,"TimeDistribution","myfunction");
  f3->SetParameters(halfLife_10,0);

  TF1 *f4 = new TF1("f4",fptr,&TimeDistribution::oPSDecay,0,1,1,"TimeDistribution","myfunction");
  f4->SetParameters(halfLife_00,0);

  int num_bins=500;
  double x_low=0;
  double x_high=1000;
  hist_1 = new TH1F("theta_hist","Angular Distribution of Normal to Decay Plane",num_bins,x_low,x_high);
  sampling_points=100000;

  if(Distr.compare("pPSDecay")==0)
    {    
      hist_1->FillRandom("f1",sampling_points);
    }
  if(Distr.compare("oPSDecay")==0)
    {      
      hist_1->FillRandom("f2",sampling_points);
    }
  // Histograms for nonzero magnetic fields
  if(Distr.compare("pPSDecay_2Gamma")==0)
    {
      hist_1->FillRandom("f3",sampling_points);
    }
  if(Distr.compare("oPSDecay_2Gamma")==0)
    {
      hist_1->FillRandom("f4",sampling_points);
    }
  if(Distr.compare("pPSDecay_3Gamma")==0)
    {
      hist_1->FillRandom("f3",sampling_points);
    }
  if(Distr.compare("oPSDecay_3Gamma")==0)
    {
      hist_1->FillRandom("f4",sampling_points);
    }
  
  hist_1->GetXaxis()->SetTitle("theta");

  //Get a value from the distribution
  val=hist_1->GetRandom();


  TFile * myfile = new TFile("fillrandom.root","RECREATE");
  hist_1->Write();
  //Get a value from the distribution
  val=hist_1->GetRandom();


  myfile->Close();

  
  int z;


  delete f1;
  delete f2;
  delete fptr;
  return hist_1;

}

