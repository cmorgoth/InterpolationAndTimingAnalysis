#include <iostream>
#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLatex.h>

int main ( int argc, char** argv)
{
  gROOT->Reset();

  TFile* f = new TFile("TBAnalysis.root","READ");
  TTree* t = (TTree*)f->Get("t_tree");
  t->Draw("ch2_t0-ch1_t0>>tmp1(100,-20,20)", "ch1Amp>0");
  TH1F* h = (TH1F*)gDirectory->Get("tmp1");

  //Define Modifed Gaussian
  TF1* f1 = new TF1("f1","[0]*([1]/2.*exp([1]/2.*(2.*[2]+[1]*[3]*[3]-2.*x)))*TMath::Erfc(([2]+[1]*[3]*[3]-x)/(TMath::Sqrt(2.)*[3]))", -5,5);
  f1->SetParameter(0,2000);
  f1->SetParameter(1,5);
  f1->SetParameter(2,0);
  f1->SetParameter(3,1.5);
  TCanvas* c = new TCanvas("c", "c", 800, 600);
  h->Fit("f1", "R");
  double max = f1->GetMaximum();
  double step = 0.01;
  double t_low = -20;
  double t_high = 20;
  double t_fwhm_low = 0;
  double t_fwhm_high = 0;
  for ( unsigned int i = 0; i < 10000; i++ )
  {
    double t = t_low + double(i)*step;
    if ( f1->Eval(t) > max/2.0)
    {
      t_fwhm_low = t;
      break;
    }
  }

  for ( unsigned int i = 0; i < 10000; i++ )
  {
    double t = t_high - double(i)*step;
    if ( f1->Eval(t) > max/2.0)
    {
      t_fwhm_high = t;
      break;
    }
  }

  std::cout << t_fwhm_low << " " << t_fwhm_high << " " << t_fwhm_high-t_fwhm_low << std::endl;

  h->SetStats(0);
  h->SetTitle("");
  h->SetXTitle("Delay [ps]");
  h->SetYTitle("Events / (0.4 ps)");
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlue);
  h->SetMarkerSize(1.5);
  //c->SetLogy();
  h->Draw("P");
  f1->Draw("same");
  TLatex* tex = new TLatex(0.9,0.92,Form("FWHM = %.2f [ps]", t_fwhm_high-t_fwhm_low));
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(41);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.2,0.92, "532 nm");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(41);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex->Draw();
  c->SaveAs("deltaT.pdf");

  f->Close();
  return 0;
}
