{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Aug 26 12:27:43 2014) by ROOT version5.34/18
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   c1->Range(-0.125,0.3558439,1.125,3.787105);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *crystal_23 = new TH1F("crystal_23","Energy Distribution",30,0,1);
   crystal_23->SetBinContent(1,1467);
   crystal_23->SetBinContent(2,111);
   crystal_23->SetBinContent(3,30);
   crystal_23->SetBinContent(4,108);
   crystal_23->SetBinContent(5,60);
   crystal_23->SetBinContent(6,138);
   crystal_23->SetBinContent(7,70);
   crystal_23->SetBinContent(8,60);
   crystal_23->SetBinContent(9,131);
   crystal_23->SetBinContent(10,92);
   crystal_23->SetBinContent(11,63);
   crystal_23->SetBinContent(12,20);
   crystal_23->SetBinContent(13,40);
   crystal_23->SetBinContent(14,60);
   crystal_23->SetBinContent(15,50);
   crystal_23->SetBinContent(16,886);
   crystal_23->SetBinContent(17,10);
   crystal_23->SetBinContent(18,10);
   crystal_23->SetBinContent(19,10);
   crystal_23->SetBinContent(20,20);
   crystal_23->SetBinContent(21,10);
   crystal_23->SetBinContent(22,20);
   crystal_23->SetBinContent(23,10);
   crystal_23->SetBinContent(24,40);
   crystal_23->SetBinContent(26,10);
   crystal_23->SetBinContent(27,49);
   crystal_23->SetBinContent(28,20);
   crystal_23->SetBinContent(30,40);
   crystal_23->SetBinContent(31,408);
   crystal_23->SetEntries(4043);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("crystal_23");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 4043   ");
   text = ptstats->AddText("Mean  = 0.2485");
   text = ptstats->AddText("RMS   = 0.2614");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   crystal_23->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(crystal_23);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   crystal_23->SetLineColor(ci);
   crystal_23->GetXaxis()->SetTitle("Energy (MeV)");
   crystal_23->GetXaxis()->SetLabelFont(42);
   crystal_23->GetXaxis()->SetLabelSize(0.035);
   crystal_23->GetXaxis()->SetTitleSize(0.035);
   crystal_23->GetXaxis()->SetTitleFont(42);
   crystal_23->GetYaxis()->SetLabelFont(42);
   crystal_23->GetYaxis()->SetLabelSize(0.035);
   crystal_23->GetYaxis()->SetTitleSize(0.035);
   crystal_23->GetYaxis()->SetTitleFont(42);
   crystal_23->GetZaxis()->SetLabelFont(42);
   crystal_23->GetZaxis()->SetLabelSize(0.035);
   crystal_23->GetZaxis()->SetTitleSize(0.035);
   crystal_23->GetZaxis()->SetTitleFont(42);
   crystal_23->Draw("");
   
   TPaveText *pt = new TPaveText(0.3305172,0.9339831,0.6694828,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("Energy Distribution");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
