#ifndef TIMEDISTRIBUTION_H
#define TIMEDISTRIBUTION_H
#include <TROOT.h>
#include <cstring>
#include <TH1F.h>
//Any other important ROOT header files

class TimeDistribution{

public:
  TimeDistribution();
  ~TimeDistribution();

  Double_t pPSDecay(Double_t* x, Double_t* par);
  Double_t oPSDecay(Double_t* x, Double_t* par);
  Double_t GetVal();
  Double_t GetBranchingRatio(std::string,double);
  //Build histogram from theta disribution and return a random angle
  TH1F* buildHist(std::string,double);

private:
  TH1F* hist_1;
  Double_t f;
  int x_low;
  int x_high;
  int num_bins;
  int sampling_points;
  Double_t val;
  Double_t angle;
};
#endif
