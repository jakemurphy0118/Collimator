{
#include <sstream>
#include <fstream>
#include <string>
  gROOT->SetBatch(kTRUE);
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TFile f("1.root");
  f->cd(0);
  TIter nextkey(f.GetListOfKeys());
  TKey *key;
  //  std::string namestring;
  std::string classstring;
  ostringstream ss("");
  int nHist=0;
  char namestring[16];
  while ((key=(TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
    printf("at key:%s, object
class:%s\n",key->GetName(),obj->ClassName());
    //    namestring=obj->GetName();
    classstring=obj->ClassName();
    if(classstring.compare("TH2D")==0){
      obj->Draw("colz");
      std::cout << "HUZZAH" << std::endl;
    }
    else{
      obj->Draw("hist");
    }
    sprintf ( namestring, "%d", nHist );
    std::cout << namestring << std::endl;
    strcat(namestring,".png");
    std::cout << namestring << std::endl;
    c1->SaveAs(namestring);
    nHist++;
    ss.clear();
  delete obj;
 }

  asymmetryParameter = -(hAP1->Integral(1,100) - hAP1->Integral(101,200))/hAP1->Integral();
  asymmetryParameter_wS = -(hAP1_wS->Integral(1,100) - hAP1_wS->Integral(101,200))/hAP1_wS->Integral();
  asymmetryParameter_APEX = -(hAP1_APEX->Integral(1,100) - hAP1_APEX->Integral(101,200))/hAP1_APEX->Integral();
  std::cout << "Asymmetry Parameter at Vertex: " << asymmetryParameter << " +- " << 1./TMath::Sqrt(hAP1->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter at Outside Source Holder: " << asymmetryParameter_wS << " +- " << 1./TMath::Sqrt(hAP1_wS->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter with APEX: " << asymmetryParameter_APEX << " +- " << 1./TMath::Sqrt(hAP1_APEX->Integral()) << std::endl;
}
