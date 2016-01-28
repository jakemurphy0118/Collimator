void energyHist(string a,string b,int c) {
  //Macro to make a histogram of the theta distribution
#include <fstream>
#include <sstream>

  ifstream inFile;
  stringstream ss;
  stringstream ss2;
  ss2 << b << ".eps";
  ss << "crystal_" << c;
  std::cout << ss2.str().c_str() << std::endl;
  inFile.open(a.c_str(),ios::in);
  int x_low=0;
  int x_high=1.4;
  int num_bins=30;
  TH1F *energyHist = new TH1F(ss.str().c_str(), "Energy Distribution",num_bins, x_low, x_high);
  double angle;
  double smearedenerg;
  double totenerg;
  while(!inFile.eof()){
    inFile >> totenerg;
    energyHist->Fill(totenerg);
  }

  gStyle->SetOptLogy();
  energyHist->Draw();
  energyHist->GetXaxis()->SetTitle("Energy (MeV)");
  //  gStyle->SetOptLogy();
  gStyle->SetPaperSize(18,10);
  c1->SaveAs(ss2.str().c_str());
  c1->SaveAs("Energy.C");
}
