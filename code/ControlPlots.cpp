#include <string>
#include <iostream>
#include <getopt.h>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
#include "TDatime.h"
#include "TMath.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TApplication.h"

using namespace std;

static
int g_counter = 9;

TH1D*
get_hist_1d(TFile *file, const char *name, int first, int last)
{
  TH2D *h2 = (TH2D*)file->Get(name);
  if (h2 == NULL) {
    std::cerr << "Error: failed to extract histogram \"" << name << '"' << std::endl;
    exit(EXIT_FAILURE);
  }
  TH1D *h1 = h2->ProjectionX(Form("%s-projx_%d", name, g_counter++), first, last);
  h1->SetTitle(h2->GetTitle());
  h1->SetLineColor(h2->GetLineColor());
  return h1;
}

static
int g_bin1 = 1, g_bin2 = 3;

//--->START MAIN PROGRAM
//________________________________________________________________________________
void ControlPlots(const char* data, const char* mc, const char* out)
{
  //gROOT->SetBatch(true);
  //  gROOT->cd();
  TDatime now;                                          //Set time in Root
  now.Print();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.07);
  gStyle->SetTitleSize(0.2);
  gStyle->SetTitleOffset(0.25);

  /////Read files --1-->
  TFile *fdata = new TFile(data);  
  TFile *fmc   = new TFile(mc); 

  //TFile *fdata = new TFile("040506e.root");  
  //TFile *fmc   = new TFile("040506e_MC.root"); 
  /////END Read files <--1--

  Bool_t REBIN = 0;
  Int_t binXX = 5;

  ////Define Histograms --2-->

  //Event Quantities
  TH1D* hdataVertexZ   = get_hist_1d(fdata, "h2VertexZ", g_bin1, g_bin2);
  TH1D* hdataQ2        = get_hist_1d(fdata, "h2Q2", g_bin1, g_bin2);
  TH1D* hdatax         = get_hist_1d(fdata, "h2x", g_bin1, g_bin2);
  TH1D* hdataEmpz      = get_hist_1d(fdata, "h2Empz", g_bin1, g_bin2);
  TH1D* hdataGamma     = get_hist_1d(fdata, "h2Gamma", g_bin1, g_bin2);
  TH1D* hdataPtEt      = get_hist_1d(fdata, "h2PtEt", g_bin1, g_bin2);
  TH1D* hdataEtamax    = get_hist_1d(fdata, "h2Etamax", g_bin1, g_bin2);

  TH1D* hmcVertexZ     = get_hist_1d(fmc, "h2VertexZ", g_bin1, g_bin2);
  TH1D* hmcQ2          = get_hist_1d(fmc, "h2Q2", g_bin1, g_bin2);
  TH1D* hmcx           = get_hist_1d(fmc, "h2x", g_bin1, g_bin2);
  TH1D* hmcEmpz        = get_hist_1d(fmc, "h2Empz", g_bin1, g_bin2);
  TH1D* hmcGamma       = get_hist_1d(fmc, "h2Gamma", g_bin1, g_bin2);
  TH1D* hmcPtEt        = get_hist_1d(fmc, "h2PtEt", g_bin1, g_bin2);
  TH1D* hmcEtamax      = get_hist_1d(fmc, "h2Etamax", g_bin1, g_bin2);

  //Electron Quantities
  TH1D* hdataElecTheta = get_hist_1d(fdata, "h2ElecTheta", g_bin1, g_bin2);
  TH1D* hdataElecPhi   = get_hist_1d(fdata, "h2ElecPhi", g_bin1, g_bin2);
  TH1D* hdataElecE     = get_hist_1d(fdata, "h2ElecE", g_bin1, g_bin2);
  TH1D* hdataElecProb  = get_hist_1d(fdata, "h2ElecProb", g_bin1, g_bin2);
  TH1D* hdataElecy     = get_hist_1d(fdata, "h2Elecy", g_bin1, g_bin2);

  TH1D* hmcElecTheta   = get_hist_1d(fmc, "h2ElecTheta", g_bin1, g_bin2);
  TH1D* hmcElecPhi     = get_hist_1d(fmc, "h2ElecPhi", g_bin1, g_bin2);
  TH1D* hmcElecE       = get_hist_1d(fmc, "h2ElecE", g_bin1, g_bin2);
  TH1D* hmcElecProb    = get_hist_1d(fmc, "h2ElecProb", g_bin1, g_bin2);
  TH1D* hmcElecy       = get_hist_1d(fmc, "h2Elecy", g_bin1, g_bin2);

  //Jets
  TH1D* hdataJetEt     = get_hist_1d(fdata, "h2JetEt", g_bin1, g_bin2);
  TH1D* hdataJetPt     = get_hist_1d(fdata, "h2JetPt", g_bin1, g_bin2);
  TH1D* hdataJetMass   = get_hist_1d(fdata, "h2JetMass", g_bin1, g_bin2);
  TH1D* hdataJetEta    = get_hist_1d(fdata, "h2JetEta", g_bin1, g_bin2);
  TH1D* hdataJetPhi    = get_hist_1d(fdata, "h2JetPhi", g_bin1, g_bin2);

  TH1D* hmcJetEt       = get_hist_1d(fmc, "h2JetEt", g_bin1, g_bin2);
  TH1D* hmcJetPt       = get_hist_1d(fmc, "h2JetPt", g_bin1, g_bin2);
  TH1D* hmcJetMass     = get_hist_1d(fmc, "h2JetMass", g_bin1, g_bin2);
  TH1D* hmcJetEta      = get_hist_1d(fmc, "h2JetEta", g_bin1, g_bin2);
  TH1D* hmcJetPhi      = get_hist_1d(fmc, "h2JetPhi", g_bin1, g_bin2);

  //Decorrelation
  TH1D* hdataDecorrPhi[5];
  TH2D* h2dataPtDecorrPhi[5];
  TH1D* hmcDecorrPhi[5];
  TH2D* h2mcPtDecorrPhi[5];
  for(Int_t ijet =0; ijet < 5; ijet++){  
    h2dataPtDecorrPhi[ijet] = (TH2D*)fdata->Get(Form("h2PtDecorrPhi_%d",ijet));
    h2mcPtDecorrPhi[ijet]   = (TH2D*)fmc  ->Get(Form("h2PtDecorrPhi_%d" ,ijet));

    hdataDecorrPhi[ijet] =
      h2dataPtDecorrPhi[ijet]->ProjectionX(Form("hdataDecorrPhi_%d", ijet));
    hmcDecorrPhi[ijet] =
      h2mcPtDecorrPhi[ijet]->ProjectionX(Form("hmcDecorrPhi_%d", ijet));
  }
  ////End Define Histograms <--2--

  cout << Form( "Number of Events Data: %.0f",hdataVertexZ->GetEntries()) << endl;
  cout << Form( "Number of Events MC  : %.0f",hmcVertexZ->GetEntries()) << endl;

  //////Event Plots --3-->
  TCanvas *c0 = new TCanvas("c0","Event",0,0,1300*2,900*2);
  c0->cd();
  c0->Divide(3,2);

  c0->cd(1);
  //c0->cd(2)->SetLogy();
  TPad *pad000 = new TPad("pad000","pad000", 0.0,0.3,1.0,1.0);
  //pad000->SetLogy();
  pad000->SetBottomMargin(0);
  pad000->Draw();
  pad000->cd();
  if(REBIN) hdataVertexZ->Rebin(binXX);
  hdataVertexZ->SetTitle("Vertex Z");
  hdataVertexZ->GetXaxis()->SetTitle("z (cm))");
  hdataVertexZ->SetLineWidth(2); hdataVertexZ->SetLineColor(4);
  //hdataVertexZ->Scale(1/hdataVertexZ->GetEntries());
  //hdataVertexZ->SetMaximum(470000);
  //hdataVertexZ->SetMinimum(0);
  hdataVertexZ->Draw();
  TLegend *legend0 = new TLegend(0.7,0.70,0.85,0.85);
  legend0->AddEntry(hdataVertexZ, "Data", "l");
  //TLatex* latex1 = new TLatex();
  //latex1->SetTextFont(22);
  //latex1->SetTextSize(0.1);
  //latex1->DrawLatexNDC(0.6,0.60,trigger);
  //if(ijet==0)  latex1->DrawLatexNDC(0.6,0.50,"High jet");
  //else  latex1->DrawLatexNDC(0.6,0.50,"Low jet");
  //
  if(REBIN) hmcVertexZ->Rebin(binXX);
  hmcVertexZ->SetLineColor(2);
  hmcVertexZ->Scale(hdataVertexZ->GetEntries()/hmcVertexZ->GetEntries());
  hmcVertexZ->Draw("same");
  legend0->AddEntry(hmcVertexZ,"MC", "l");
  legend0->SetTextSize(0.05);
  legend0->SetBorderSize(0);
  legend0->SetFillStyle(0);
  legend0->Draw("hist p");
  /////
  c0->cd(1);
  TPad *pad001 = new TPad("pad001","pad001", 0,0.05,1.0,0.3);
  pad001->SetTopMargin(0);
  pad001->SetBottomMargin(0.2);
  pad001->Draw();
  pad001->cd();
  TH1D *hdataVertexZ_R= (TH1D*) hdataVertexZ->Clone();
  //hdataVertexZ_R->Add(hmcVertexZ,-1);
  hdataVertexZ_R->Divide(hmcVertexZ);
  //hdataVertexZ_R->Rebin(5); 
  hdataVertexZ_R->SetBit(TH1::kNoTitle);
  hdataVertexZ_R->SetLineColor(1); 
  hdataVertexZ_R->GetYaxis()->SetTitle("Data/MC");
  hdataVertexZ_R->GetYaxis()->CenterTitle();
  hdataVertexZ_R->GetYaxis()->SetTitleSize(0.07);
  hdataVertexZ_R->GetYaxis()->SetTitleOffset(0.38);
  hdataVertexZ_R->GetYaxis()->SetRangeUser(0.8,1.2);
  //hdataVertexZ_R->SetLineWidth(2);
  hdataVertexZ_R->GetYaxis()->SetLabelSize(0.1);
  hdataVertexZ_R->GetXaxis()->SetLabelSize(0.1);
  hdataVertexZ_R->GetXaxis()->SetTitleSize(0.07);
  hdataVertexZ_R->GetYaxis()->SetNdivisions(5); 
  hdataVertexZ_R->SetMarkerStyle(2);
  hdataVertexZ_R->Draw("E");



  c0->cd(2);
  //c0->cd(2)->SetLogy();
  TPad *pad010 = new TPad("pad010","pad010", 0.0,0.3,1.0,1.0);
  pad010->SetLogy();
  pad010->SetBottomMargin(0);
  pad010->Draw();
  pad010->cd();
  if(REBIN) hdataQ2->Rebin(binXX);
  hdataQ2->SetTitle("Q^{2}");
  hdataQ2->GetXaxis()->SetTitle("Q^{2} (GeV)");
  hdataQ2->SetLineWidth(2);
  hdataQ2->SetLineColor(4);
  hdataQ2->Draw();
  if(REBIN) hmcQ2->Rebin(binXX);
  hmcQ2->SetLineColor(2);
  hmcQ2->Scale(hdataQ2->GetEntries()/hmcQ2->GetEntries());
  hmcQ2->Draw("same");
  legend0->Draw("hist p");
  /////
  c0->cd(2);
  TPad *pad011 = new TPad("pad011","pad011", 0,0.05,1.0,0.3);
  pad011->SetTopMargin(0);
  pad011->SetBottomMargin(0.2);
  pad011->Draw();
  pad011->cd();
  TH1D *hdataQ2_R= (TH1D*) hdataQ2->Clone();
  //hdataQ2_R->Add(hmcQ2,-1);
  hdataQ2_R->Divide(hmcQ2);
  //hdataQ2_R->Rebin(5); 
  hdataQ2_R->SetBit(TH1::kNoTitle);
  hdataQ2_R->SetLineColor(1); 
  hdataQ2_R->GetYaxis()->SetTitle("Data/MC");
  hdataQ2_R->GetYaxis()->CenterTitle();
  hdataQ2_R->GetYaxis()->SetTitleSize(0.07);
  hdataQ2_R->GetYaxis()->SetTitleOffset(0.38);
  hdataQ2_R->GetYaxis()->SetRangeUser(0.8,1.2);
  //hdataQ2_R->SetLineWidth(2);
  hdataQ2_R->GetYaxis()->SetLabelSize(0.1);
  hdataQ2_R->GetXaxis()->SetLabelSize(0.1);
  hdataQ2_R->GetXaxis()->SetTitleSize(0.07);
  hdataQ2_R->GetYaxis()->SetNdivisions(5); 
  hdataQ2_R->SetMarkerStyle(2);
  hdataQ2_R->Draw("E");



  c0->cd(3);
  c0->cd(3)->SetLogy();
  TPad *pad030 = new TPad("pad030","pad030", 0.0,0.3,1.0,1.0);
  pad030->SetLogy();
  pad030->SetBottomMargin(0);
  pad030->Draw();
  pad030->cd();
  if(REBIN) hdatax->Rebin(binXX);
  hdatax->SetTitle("x_{Bj}");
  hdatax->GetXaxis()->SetTitle("x");
  hdatax->SetLineWidth(2); hdatax->SetLineColor(4);
  //hdatax->Scale(1/h2datax->GetEntries());
  //hdatax->SetMaximum(10000000);
  //hdatax->SetMinimum(100);
  hdatax->Draw();
  //
  if(REBIN) hmcx->Rebin(binXX);
  hmcx->SetLineColor(2);
  hmcx->Scale(hdatax->GetEntries()/hmcx->GetEntries());
  hmcx->Draw("same");
  legend0->Draw("hist p");
  //
  c0->cd(3);
  TPad *pad031 = new TPad("pad031","pad031", 0,0.05,1.0,0.3);
  pad031->SetTopMargin(0);
  pad031->SetBottomMargin(0.2);
  pad031->Draw();
  pad031->cd();
  TH1D *hdatax_R= (TH1D*) hdatax->Clone();
  //hdatax_R->Add(hmcx,-1);
  hdatax_R->Divide(hmcx);
  //hdatax_R->Rebin(5); 
  hdatax_R->SetBit(TH1::kNoTitle);
  hdatax_R->SetLineColor(1); 
  hdatax_R->GetYaxis()->SetTitle("Data/MC");
  hdatax_R->GetYaxis()->CenterTitle();
  hdatax_R->GetYaxis()->SetTitleSize(0.07);
  hdatax_R->GetYaxis()->SetTitleOffset(0.38);
  hdatax_R->GetYaxis()->SetRangeUser(0.8,1.2);
  //hdatax_R->SetLineWidth(2);
  hdatax_R->GetYaxis()->SetLabelSize(0.1);
  hdatax_R->GetXaxis()->SetLabelSize(0.1);
  hdatax_R->GetXaxis()->SetTitleSize(0.07);
  hdatax_R->GetYaxis()->SetNdivisions(5); 
  hdatax_R->SetMarkerStyle(2);
  hdatax_R->Draw("E");



  c0->cd(4);
  TPad *pad040 = new TPad("pad040","pad040", 0.0,0.3,1.0,1.0);
  pad040->SetLogy();
  pad040->SetBottomMargin(0);
  pad040->Draw();
  pad040->cd();
  if(REBIN) hdataPtEt->Rebin(binXX);
  hdataPtEt->SetTitle("P_{T} / #sqrt{E_{T}}");
  hdataPtEt->GetXaxis()->SetTitle("(GeV)^{1/2}");
  //hdataPtEt->GetYaxis()->SetTitle("Q^{2}");
  //h2dataPtEt->SetLineWidth(2); h2dataPtEt->SetLineColor(4);
  //h2dataPtEt->Scale(1/h2dataPtEt->GetEntries());
  //h2dataPtEt->SetMaximum(900000);
  //h2dataPtEt->SetMinimum(0);
  hdataPtEt->Draw("");
  if(REBIN) hmcPtEt->Rebin(binXX);
  hmcPtEt->SetLineColor(2);
  hmcPtEt->Scale(hdataPtEt->GetEntries()/hmcPtEt->GetEntries());
  hmcPtEt->Draw("same");
  legend0->Draw("hist p");
  /////
  c0->cd(4);
  TPad *pad041 = new TPad("pad041","pad041", 0,0.05,1.0,0.3);
  pad041->SetTopMargin(0);
  pad041->SetBottomMargin(0.2);
  pad041->Draw();
  pad041->cd();
  TH1D *hdataPtEt_R= (TH1D*) hdataPtEt->Clone();
  //hdataPtEt_R->Add(hmcPtEt,-1);
  hdataPtEt_R->Divide(hmcPtEt);
  //hdataPtEt_R->Rebin(5); 
  hdataPtEt_R->SetBit(TH1::kNoTitle);
  hdataPtEt_R->SetLineColor(1); 
  hdataPtEt_R->GetYaxis()->SetTitle("Data/MC");
  hdataPtEt_R->GetYaxis()->CenterTitle();
  hdataPtEt_R->GetYaxis()->SetTitleSize(0.07);
  hdataPtEt_R->GetYaxis()->SetTitleOffset(0.38);
  hdataPtEt_R->GetYaxis()->SetRangeUser(-1.,3.);
  //hdataPtEt_R->SetLineWidth(2);
  hdataPtEt_R->GetYaxis()->SetLabelSize(0.1);
  hdataPtEt_R->GetXaxis()->SetLabelSize(0.1);
  hdataPtEt_R->GetXaxis()->SetTitleSize(0.07);
  hdataPtEt_R->GetYaxis()->SetNdivisions(5); 
  hdataPtEt_R->SetMarkerStyle(2);
  hdataPtEt_R->Draw("E");


  bool do_plot_etamax = true;
  c0->cd(5);
  TPad *pad050 = new TPad("pad050","pad050", 0.0,0.3,1.0,1.0);
  pad050->SetLogy();
  pad050->SetBottomMargin(0);
  pad050->Draw();
  pad050->cd();
  if (do_plot_etamax) {
    if(REBIN) hdataEtamax->Rebin(binXX);
    hdataEtamax->SetTitle("eta max.");
    hdataEtamax->GetXaxis()->SetTitle("eta max.");
    hdataEtamax->GetXaxis()->SetRangeUser(-2, 5);
    hdataEtamax->Draw("");
    if(REBIN) hmcEtamax->Rebin(binXX);
    hmcEtamax->SetLineColor(2);
    hmcEtamax->Scale(hdataEtamax->GetEntries()/hmcEtamax->GetEntries());
    hmcEtamax->Draw("same");
    legend0->Draw("hist p");
  } else {
    if(REBIN) hdataGamma->Rebin(binXX);
    hdataGamma->SetTitle("Hadronic angle");
    hdataGamma->GetXaxis()->SetTitle("cos(#gamma_{h})");
    hdataGamma->Draw("");
    if(REBIN) hmcGamma->Rebin(binXX);
    hmcGamma->SetLineColor(2);
    hmcGamma->Scale(hdataGamma->GetEntries()/hmcGamma->GetEntries());
    hmcGamma->Draw("same");
    legend0->Draw("hist p");
  }
  c0->cd(5);
  TPad *pad051 = new TPad("pad051","pad051", 0,0.05,1.0,0.3);
  pad051->SetTopMargin(0);
  pad051->SetBottomMargin(0.2);
  pad051->Draw();
  pad051->cd();
  if (do_plot_etamax) {
    TH1D *hdataEtamax_R = (TH1D*) hdataEtamax->Clone();
    hdataEtamax_R->Divide(hmcEtamax);
    hdataEtamax_R->SetBit(TH1::kNoTitle);
    hdataEtamax_R->SetLineColor(1); 
    hdataEtamax_R->GetYaxis()->SetTitle("Data/MC");
    hdataEtamax_R->GetYaxis()->CenterTitle();
    hdataEtamax_R->GetYaxis()->SetTitleSize(0.07);
    hdataEtamax_R->GetYaxis()->SetTitleOffset(0.38);
    hdataEtamax_R->GetYaxis()->SetRangeUser(-1.,3.);
    hdataEtamax_R->GetYaxis()->SetLabelSize(0.1);
    hdataEtamax_R->GetXaxis()->SetLabelSize(0.1);
    hdataEtamax_R->GetXaxis()->SetTitleSize(0.07);
    hdataEtamax_R->GetYaxis()->SetNdivisions(5); 
    hdataEtamax_R->GetXaxis()->SetRangeUser(-2, 5);
    hdataEtamax_R->SetMarkerStyle(2);
    hdataEtamax_R->Draw("E");
  } else {
    TH1D *hdataGamma_R= (TH1D*) hdataGamma->Clone();
    hdataGamma_R->Divide(hmcGamma);
    hdataGamma_R->SetBit(TH1::kNoTitle);
    hdataGamma_R->SetLineColor(1); 
    hdataGamma_R->GetYaxis()->SetTitle("Data/MC");
    hdataGamma_R->GetYaxis()->CenterTitle();
    hdataGamma_R->GetYaxis()->SetTitleSize(0.07);
    hdataGamma_R->GetYaxis()->SetTitleOffset(0.38);
    hdataGamma_R->GetYaxis()->SetRangeUser(-1.,3.);
    hdataGamma_R->GetYaxis()->SetLabelSize(0.1);
    hdataGamma_R->GetXaxis()->SetLabelSize(0.1);
    hdataGamma_R->GetXaxis()->SetTitleSize(0.07);
    hdataGamma_R->GetYaxis()->SetNdivisions(5); 
    hdataGamma_R->SetMarkerStyle(2);
    hdataGamma_R->Draw("E");
  }


  c0->cd(6);
  //c0->cd(6)->SetLogz();
  TPad *pad060 = new TPad("pad060","pad060", 0.0,0.3,1.0,1.0);
  pad060->SetLogy();
  pad060->SetBottomMargin(0);
  pad060->Draw();
  pad060->cd();
  if(REBIN) hdataEmpz->Rebin(binXX);
  hdataEmpz->SetTitle("E - p_{z}");
  hdataEmpz->GetXaxis()->SetTitle("E - p_{z}");
  //hdataEmpz->GetYaxis()->SetTitle("Q^{2}");
  //h2dataEmpz->SetLineWidth(2); h2dataEmpz->SetLineColor(4);
  //h2dataEmpz->Scale(1/h2dataEmpz->GetEntries());
  //h2dataEmpz->SetMaximum(900000);
  //h2dataEmpz->SetMinimum(0);
  hdataEmpz->Draw("");
  if(REBIN) hmcEmpz->Rebin(binXX);
  hmcEmpz->SetLineColor(2);
  hmcEmpz->Scale(hdataEmpz->GetEntries()/hmcEmpz->GetEntries());
  hmcEmpz->Draw("same");
  legend0->Draw("hist p");
  /////
  c0->cd(6);
  TPad *pad061 = new TPad("pad061","pad061", 0,0.05,1.0,0.3);
  pad061->SetTopMargin(0);
  pad061->SetBottomMargin(0.2);
  pad061->Draw();
  pad061->cd();
  TH1D *hdataEmpz_R= (TH1D*) hdataEmpz->Clone();
  //hdataEmpz_R->Add(hmcEmpz,-1);
  hdataEmpz_R->Divide(hmcEmpz);
  //hdataEmpz_R->Rebin(5); 
  hdataEmpz_R->SetBit(TH1::kNoTitle);
  hdataEmpz_R->SetLineColor(1); 
  hdataEmpz_R->GetYaxis()->SetTitle("Data/MC");
  hdataEmpz_R->GetYaxis()->CenterTitle();
  hdataEmpz_R->GetYaxis()->SetTitleSize(0.07);
  hdataEmpz_R->GetYaxis()->SetTitleOffset(0.38);
  hdataEmpz_R->GetYaxis()->SetRangeUser(-1.,3.);
  //hdataEmpz_R->SetLineWidth(2);
  hdataEmpz_R->GetYaxis()->SetLabelSize(0.1);
  hdataEmpz_R->GetXaxis()->SetLabelSize(0.1);
  hdataEmpz_R->GetXaxis()->SetTitleSize(0.07);
  hdataEmpz_R->GetYaxis()->SetNdivisions(5); 
  hdataEmpz_R->SetMarkerStyle(2);
  hdataEmpz_R->Draw("E");

  //////ENDEvent Plots <--3--



  //////Electron Plots --4-->
  //TCanvas *ce = new TCanvas("ce","Lepton",0,0,1300*2,900*2);
  TCanvas *ce = new TCanvas("ce","Lepton",0,0,950*2,900*2);
  ce->cd();
  //ce->Divide(3,2);
  ce->Divide(2,2);

  ce->cd(1);
  TPad *pad020 = new TPad("pad020","pad020", 0.0,0.3,1.0,1.0);
  pad020->SetBottomMargin(0);
  pad020->Draw();
  pad020->cd();
  if(REBIN) hdataElecTheta->Rebin(binXX);
  hdataElecTheta->SetTitle("Lepton (#theta) angle");
  hdataElecTheta->GetXaxis()->SetTitle("(degrees)");
  hdataElecTheta->SetLineWidth(2); hdataElecTheta->SetLineColor(4);
  hdataElecTheta->Draw();
  //
  if(REBIN) hmcElecTheta->Rebin(binXX);
  hmcElecTheta->SetLineColor(2);
  hmcElecTheta->Scale(hdataElecTheta->GetEntries()/hmcElecTheta->GetEntries());
  hmcElecTheta->Draw("same");
  legend0->Draw("hist p");
  ////
  ce->cd(1);
  TPad *pad021 = new TPad("pad021","pad021", 0,0.05,1.0,0.3);
  pad021->SetTopMargin(0);
  pad021->SetBottomMargin(0.2);
  pad021->Draw();
  pad021->cd();
  TH1D *hdataElecTheta_R= (TH1D*) hdataElecTheta->Clone();
  hdataElecTheta_R->Divide(hmcElecTheta);
  hdataElecTheta_R->SetBit(TH1::kNoTitle);
  hdataElecTheta_R->SetLineColor(1); 
  hdataElecTheta_R->GetYaxis()->SetTitle("Data/MC");
  hdataElecTheta_R->GetYaxis()->CenterTitle();
  hdataElecTheta_R->GetYaxis()->SetTitleSize(0.07);
  hdataElecTheta_R->GetYaxis()->SetTitleOffset(0.38);
  hdataElecTheta_R->GetYaxis()->SetRangeUser(0.8,1.2);
  hdataElecTheta_R->GetYaxis()->SetLabelSize(0.1);
  hdataElecTheta_R->GetXaxis()->SetLabelSize(0.1);
  hdataElecTheta_R->GetXaxis()->SetTitleSize(0.07);
  hdataElecTheta_R->GetYaxis()->SetNdivisions(5); 
  hdataElecTheta_R->SetMarkerStyle(2);
  hdataElecTheta_R->Draw("E");



  ce->cd(2);
  TPad *pad220 = new TPad("pad220","pad220", 0.0,0.3,1.0,1.0);
  pad220->SetBottomMargin(0);
  pad220->Draw();
  pad220->cd();
  if(REBIN) hdataElecE->Rebin(binXX);
  hdataElecE->SetTitle("Lepton energy");
  hdataElecE->GetXaxis()->SetTitle("E (GeV)");
  hdataElecE->SetLineWidth(2); hdataElecE->SetLineColor(4);
  hdataElecE->Draw();
  //
  if(REBIN) hmcElecE->Rebin(binXX);
  hmcElecE->SetLineColor(2);
  hmcElecE->Scale(hdataElecE->GetEntries()/hmcElecE->GetEntries());
  hmcElecE->Draw("same");
  legend0->Draw("hist p");
  ////
  ce->cd(2);
  TPad *pad221 = new TPad("pad221","pad221", 0,0.05,1.0,0.3);
  pad221->SetTopMargin(0);
  pad221->SetBottomMargin(0.2);
  pad221->Draw();
  pad221->cd();
  TH1D *hdataElecE_R= (TH1D*) hdataElecE->Clone();
  //hdataElecE_R->Add(hmcElecE,-1);
  hdataElecE_R->Divide(hmcElecE);
  //hdataElecE_R->Rebin(5); 
  hdataElecE_R->SetBit(TH1::kNoTitle);
  hdataElecE_R->SetLineColor(1); 
  hdataElecE_R->GetYaxis()->SetTitle("Data/MC");
  hdataElecE_R->GetYaxis()->CenterTitle();
  hdataElecE_R->GetYaxis()->SetTitleSize(0.07);
  hdataElecE_R->GetYaxis()->SetTitleOffset(0.38);
  hdataElecE_R->GetYaxis()->SetRangeUser(0.,2.);
  //hdataElecE_R->SetLineWidth(2);
  hdataElecE_R->GetYaxis()->SetLabelSize(0.1);
  hdataElecE_R->GetXaxis()->SetLabelSize(0.1);
  hdataElecE_R->GetXaxis()->SetTitleSize(0.07);
  hdataElecE_R->GetYaxis()->SetNdivisions(5); 
  hdataElecE_R->SetMarkerStyle(2);
  hdataElecE_R->Draw("E");




  //ce->cd(5);
  //TPad *pad520 = new TPad("pad520","pad520", 0.0,0.3,1.0,1.0);
  //pad520->SetLogy();
  //pad520->SetBottomMargin(0);
  //pad520->Draw();
  //pad520->cd();
  //if(REBIN) hdataElecProb->Rebin(binXX);
  //hdataElecProb->SetTitle("Probability of being DIS lepton");
  //hdataElecProb->SetLineWidth(2); hdataElecProb->SetLineColor(4);
  //hdataElecProb->Draw();
  ////
  //if(REBIN) hmcElecProb->Rebin(binXX);
  //hmcElecProb->SetLineColor(2);
  //hmcElecProb->Scale(hdataElecProb->GetEntries()/hmcElecProb->GetEntries());
  //hmcElecProb->Draw("same");
  //legend0->Draw("hist p");
  //////
  //ce->cd(5);
  //TPad *pad521 = new TPad("pad521","pad521", 0,0.05,1.0,0.3);
  //pad521->SetTopMargin(0);
  //pad521->SetBottomMargin(0.2);
  //pad521->Draw();
  //pad521->cd();
  //TH1D *hdataElecProb_R= (TH1D*) hdataElecProb->Clone();
  //hdataElecProb_R->Divide(hmcElecProb);
  //hdataElecProb_R->SetBit(TH1::kNoTitle);
  //hdataElecProb_R->SetLineColor(1); 
  //hdataElecProb_R->GetYaxis()->SetTitle("Data/MC");
  //hdataElecProb_R->GetYaxis()->CenterTitle();
  //hdataElecProb_R->GetYaxis()->SetTitleSize(0.07);
  //hdataElecProb_R->GetYaxis()->SetTitleOffset(0.38);
  //hdataElecProb_R->GetYaxis()->SetRangeUser(-1.,3.);
  //hdataElecProb_R->GetYaxis()->SetLabelSize(0.1);
  //hdataElecProb_R->GetXaxis()->SetLabelSize(0.1);
  //hdataElecProb_R->GetXaxis()->SetTitleSize(0.07);
  //hdataElecProb_R->GetYaxis()->SetNdivisions(5); 
  //hdataElecProb_R->SetMarkerStyle(2);
  //hdataElecProb_R->Draw("E");


  //ce->cd(3);
  //ce->cd(3)->SetLogz();

  ce->cd(4);
  //  ce->cd(5)->SetLogy();
  TPad *pad620 = new TPad("pad620","pad620", 0.0,0.3,1.0,1.0);
  //pad620->SetLogy();
  pad620->SetBottomMargin(0);
  pad620->Draw();
  pad620->cd();
  if(REBIN) hdataElecy->Rebin(binXX);
  hdataElecy->SetTitle("Lepton Inelasticity");
  hdataElecy->GetXaxis()->SetTitle("y");
  hdataElecy->SetLineWidth(2); hdataElecy->SetLineColor(4);
  hdataElecy->Draw();
  //
  if(REBIN) hmcElecy->Rebin(binXX);
  hmcElecy->SetLineColor(2);
  hmcElecy->Scale(hdataElecy->GetEntries()/hmcElecy->GetEntries());
  hmcElecy->Draw("same");
  legend0->Draw("hist p");
  ////
  ce->cd(4);
  TPad *pad621 = new TPad("pad621","pad621", 0,0.05,1.0,0.3);
  pad621->SetTopMargin(0);
  pad621->SetBottomMargin(0.2);
  pad621->Draw();
  pad621->cd();
  TH1D *hdataElecy_R= (TH1D*) hdataElecy->Clone();
  hdataElecy_R->Divide(hmcElecy);
  hdataElecy_R->SetBit(TH1::kNoTitle);
  hdataElecy_R->SetLineColor(1); 
  hdataElecy_R->GetYaxis()->SetTitle("Data/MC");
  hdataElecy_R->GetYaxis()->CenterTitle();
  hdataElecy_R->GetYaxis()->SetTitleSize(0.07);
  hdataElecy_R->GetYaxis()->SetTitleOffset(0.38);
  hdataElecy_R->GetYaxis()->SetRangeUser(0.5,1.5);
  hdataElecy_R->GetYaxis()->SetLabelSize(0.1);
  hdataElecy_R->GetXaxis()->SetLabelSize(0.1);
  hdataElecy_R->GetXaxis()->SetTitleSize(0.07);
  hdataElecy_R->GetYaxis()->SetNdivisions(5); 
  hdataElecy_R->SetMarkerStyle(2);
  hdataElecy_R->Draw("E");

  //////ENDElectron Plots <--4--



  //////Jets plots --5-->
  TCanvas *c1 = new TCanvas("c1","Jets",0,0,950*2,900*2);
  c1->cd();
  c1->Divide(1,2);

  c1->cd(1);
  if(REBIN) hdataElecPhi->Rebin(binXX);
  TPad *pad420 = new TPad("pad420","pad420", 0.0,0.3,1.0,1.0);
  pad420->SetBottomMargin(0);
  pad420->Draw();
  pad420->cd();
  hdataElecPhi->SetTitle("Azimuthal angle (#phi) of the lepton");
  hdataElecPhi->GetXaxis()->SetTitle("(rad)");
  hdataElecPhi->SetLineWidth(2); hdataElecPhi->SetLineColor(4);
  hdataElecPhi->Draw();
  //
  if(REBIN) hmcElecPhi->Rebin(binXX);
  hmcElecPhi->SetLineColor(2);
  hmcElecPhi->Scale(hdataElecPhi->GetEntries()/hmcElecPhi->GetEntries());
  hmcElecPhi->Draw("same");
  legend0->Draw("hist p");
  ////
  c1->cd(1);
  TPad *pad421 = new TPad("pad421","pad421", 0,0.05,1.0,0.3);
  pad421->SetTopMargin(0);
  pad421->SetBottomMargin(0.2);
  pad421->Draw();
  pad421->cd();
  TH1D *hdataElecPhi_R= (TH1D*) hdataElecPhi->Clone();
  hdataElecPhi_R->Divide(hmcElecPhi);
  hdataElecPhi_R->SetBit(TH1::kNoTitle);
  hdataElecPhi_R->SetLineColor(1); 
  hdataElecPhi_R->GetYaxis()->SetTitle("Data/MC");
  hdataElecPhi_R->GetYaxis()->CenterTitle();
  hdataElecPhi_R->GetYaxis()->SetTitleSize(0.07);
  hdataElecPhi_R->GetYaxis()->SetTitleOffset(0.38);
  hdataElecPhi_R->GetYaxis()->SetRangeUser(0.8,1.2);
  hdataElecPhi_R->GetYaxis()->SetLabelSize(0.1);
  hdataElecPhi_R->GetXaxis()->SetLabelSize(0.1);
  hdataElecPhi_R->GetXaxis()->SetTitleSize(0.07);
  hdataElecPhi_R->GetYaxis()->SetNdivisions(5); 
  hdataElecPhi_R->SetMarkerStyle(2);
  hdataElecPhi_R->Draw("E");


  //c1->cd(1);
  //TPad *padJ30 = new TPad("padJ30","padJ30", 0.0,0.3,1.0,1.0);
  //padJ30->SetLogy();
  //padJ30->SetBottomMargin(0);
  //padJ30->Draw();
  //padJ30->cd();
  //if(REBIN) hdataJetEt->Rebin(binXX);
  //hdataJetEt->SetTitle("Jets transverse Energy");
  //hdataJetEt->GetXaxis()->SetTitle("E_{T} (GeV)");
  //hdataJetEt->SetLineWidth(2); hdataJetEt->SetLineColor(4);
  //hdataJetEt->SetMaximum(10000000);
  //hdataJetEt->SetMinimum(1);
  //hdataJetEt->Draw();
  ////
  //if(REBIN) hmcJetEt->Rebin(binXX);
  //hmcJetEt->SetLineColor(2);
  //hmcJetEt->Scale(hdataJetEt->GetEntries()/hmcJetEt->GetEntries());
  //hmcJetEt->Draw("same");
  //legend0->Draw("hist p");
  ////
  //c1->cd(1);
  //TPad *padJ31 = new TPad("padJ31","padJ31", 0,0.05,1.0,0.3);
  //padJ31->SetTopMargin(0);
  //padJ31->SetBottomMargin(0.2);
  //padJ31->Draw();
  //padJ31->cd();
  //TH1D *hdataJetEt_R= (TH1D*) hdataJetEt->Clone();
  //hdataJetEt_R->Divide(hmcJetEt);
  //hdataJetEt_R->SetBit(TH1::kNoTitle);
  //hdataJetEt_R->SetLineColor(1); 
  //hdataJetEt_R->GetYaxis()->SetTitle("Data/MC");
  //hdataJetEt_R->GetYaxis()->CenterTitle();
  //hdataJetEt_R->GetYaxis()->SetTitleSize(0.07);
  //hdataJetEt_R->GetYaxis()->SetTitleOffset(0.38);
  //hdataJetEt_R->GetYaxis()->SetRangeUser(0.,2.);
  ////hdataJetEt_R->SetLineWidth(2);
  //hdataJetEt_R->GetYaxis()->SetLabelSize(0.1);
  //hdataJetEt_R->GetXaxis()->SetLabelSize(0.1);
  //hdataJetEt_R->GetXaxis()->SetTitleSize(0.07);
  //hdataJetEt_R->GetYaxis()->SetNdivisions(5); 
  //hdataJetEt_R->SetMarkerStyle(2);
  //hdataJetEt_R->Draw("E");




  //c1->cd(2);
  //TPad *pad230 = new TPad("pad230","pad230", 0.0,0.3,1.0,1.0);
  //pad230->SetLogy();
  //pad230->SetBottomMargin(0);
  //pad230->Draw();
  //pad230->cd();
  //if(REBIN) hdataJetPt->Rebin(binXX);
  //hdataJetPt->SetTitle("Transverse Momentun of the Jets");
  //hdataJetPt->GetXaxis()->SetTitle("Pt(GeV/c)");
  //hdataJetPt->SetLineWidth(2); hdataJetPt->SetLineColor(4);
  //hdataJetPt->SetMinimum(1);
  //hdataJetPt->Draw();
  ////
  //if(REBIN) hmcJetPt->Rebin(binXX);
  //hmcJetPt->SetLineColor(2);
  //hmcJetPt->Scale(hdataJetPt->GetEntries()/hmcJetPt->GetEntries());
  //hmcJetPt->Draw("same");
  //legend0->Draw("hist p");
  //c1->cd(2);
  //TPad *pad231 = new TPad("pad231","pad231", 0,0.05,1.0,0.3);
  //pad231->SetTopMargin(0);
  //pad231->SetBottomMargin(0.2);
  //pad231->Draw();
  //pad231->cd();
  //TH1D *hdataJetPt_R= (TH1D*) hdataJetPt->Clone();
  //hdataJetPt_R->Divide(hmcJetPt);
  //hdataJetPt_R->SetBit(TH1::kNoTitle);
  //hdataJetPt_R->SetLineColor(1); 
  //hdataJetPt_R->GetYaxis()->SetTitle("Data/MC");
  //hdataJetPt_R->GetYaxis()->CenterTitle();
  //hdataJetPt_R->GetYaxis()->SetTitleSize(0.07);
  //hdataJetPt_R->GetYaxis()->SetTitleOffset(0.38);
  //hdataJetPt_R->GetYaxis()->SetRangeUser(0.,2.);
  //hdataJetPt_R->GetYaxis()->SetLabelSize(0.1);
  //hdataJetPt_R->GetXaxis()->SetLabelSize(0.1);
  //hdataJetPt_R->GetXaxis()->SetTitleSize(0.07);
  //hdataJetPt_R->GetYaxis()->SetNdivisions(5); 
  //hdataJetPt_R->SetMarkerStyle(2);
  //hdataJetPt_R->Draw("E");


  //c1->cd(3);
  //TPad *pad330 = new TPad("pad330","pad330", 0.0,0.3,1.0,1.0);
  //pad330->SetBottomMargin(0);
  //pad330->Draw();
  //pad330->cd();
  //if(REBIN) hdataJetEta->Rebin(binXX);
  //hdataJetEta->SetTitle("Jet pseudo-rapidity");
  //hdataJetEta->GetXaxis()->SetTitle("#eta");
  //hdataJetEta->SetLineWidth(2); hdataJetEta->SetLineColor(4);
  //hdataJetEta->Draw();
  ////
  //if(REBIN) hmcJetEta->Rebin(binXX);
  //hmcJetEta->SetLineColor(2);
  //hmcJetEta->Scale(hdataJetEta->GetEntries()/hmcJetEta->GetEntries());
  //hmcJetEta->Draw("same");
  //legend0->Draw("hist p");

  //c1->cd(3);
  //TPad *pad331 = new TPad("pad331","pad331", 0,0.05,1.0,0.3);
  //pad331->SetTopMargin(0);
  //pad331->SetBottomMargin(0.2);
  //pad331->Draw();
  //pad331->cd();
  //TH1D *hdataJetEta_R= (TH1D*) hdataJetEta->Clone();
  //hdataJetEta_R->Divide(hmcJetEta);
  //hdataJetEta_R->SetBit(TH1::kNoTitle);
  //hdataJetEta_R->SetLineColor(1); 
  //hdataJetEta_R->GetYaxis()->SetTitle("Data/MC");
  //hdataJetEta_R->GetYaxis()->CenterTitle();
  //hdataJetEta_R->GetYaxis()->SetTitleSize(0.07);
  //hdataJetEta_R->GetYaxis()->SetTitleOffset(0.38);
  //hdataJetEta_R->GetYaxis()->SetRangeUser(0.9,1.1);
  //hdataJetEta_R->GetYaxis()->SetLabelSize(0.1);
  //hdataJetEta_R->GetXaxis()->SetLabelSize(0.1);
  //hdataJetEta_R->GetXaxis()->SetTitleSize(0.07);
  //hdataJetEta_R->GetYaxis()->SetNdivisions(5); 
  //hdataJetEta_R->SetMarkerStyle(2);
  //hdataJetEta_R->Draw("E");

  c1->cd(2);
  TPad *pad430 = new TPad("pad430","pad430", 0.0,0.3,1.0,1.0);
  pad430->SetBottomMargin(0);
  pad430->Draw();
  pad430->cd();
  if(REBIN) hdataJetPhi->Rebin(binXX);
  hdataJetPhi->SetTitle("Jet azimuth angle");
  hdataJetPhi->GetXaxis()->SetTitle("#phi");
  hdataJetPhi->SetLineWidth(2); hdataJetPhi->SetLineColor(4);
  hdataJetPhi->Draw();
  //
  if(REBIN) hmcJetPhi->Rebin(binXX);
  hmcJetPhi->SetLineColor(2);
  hmcJetPhi->Scale(hdataJetPhi->GetEntries()/hmcJetPhi->GetEntries());
  hmcJetPhi->Draw("same");
  legend0->Draw("hist p");
  c1->cd(2);
  TPad *pad431 = new TPad("pad431","pad431", 0,0.05,1.0,0.3);
  pad431->SetTopMargin(0);
  pad431->SetBottomMargin(0.2);
  pad431->Draw();
  pad431->cd();
  TH1D *hdataJetPhi_R= (TH1D*) hdataJetPhi->Clone();
  hdataJetPhi_R->Divide(hmcJetPhi);
  hdataJetPhi_R->SetBit(TH1::kNoTitle);
  hdataJetPhi_R->SetLineColor(1); 
  hdataJetPhi_R->GetYaxis()->SetTitle("Data/MC");
  hdataJetPhi_R->GetYaxis()->CenterTitle();
  hdataJetPhi_R->GetYaxis()->SetTitleSize(0.07);
  hdataJetPhi_R->GetYaxis()->SetTitleOffset(0.38);
  hdataJetPhi_R->GetYaxis()->SetRangeUser(0.9,1.1);
  hdataJetPhi_R->GetYaxis()->SetLabelSize(0.1);
  hdataJetPhi_R->GetXaxis()->SetLabelSize(0.1);
  hdataJetPhi_R->GetXaxis()->SetTitleSize(0.07);
  hdataJetPhi_R->GetYaxis()->SetNdivisions(5); 
  hdataJetPhi_R->SetMarkerStyle(2);
  hdataJetPhi_R->Draw("E");



  ///////Decorrelation
  TCanvas *cQ = new TCanvas("cQ","Data/MC",0,0,950*2,900*2);
  cQ->cd();
  cQ->Divide(1,2);

  cQ->cd(1);
  if(REBIN) hdataDecorrPhi[0]->Rebin(binXX);
  hdataDecorrPhi[0]->SetTitle("Jet - lepton decorrelation angle");
  hdataDecorrPhi[0]->GetXaxis()->SetTitle("(rad)");
  hdataDecorrPhi[0]->SetLineWidth(2); hdataDecorrPhi[0]->SetLineColor(4);
  hdataDecorrPhi[0]->Draw();
  //
  if(REBIN) hmcDecorrPhi[0]->Rebin(binXX);
  hmcDecorrPhi[0]->SetLineColor(2);
  hmcDecorrPhi[0]->Scale(hdataDecorrPhi[0]->GetEntries()/hmcDecorrPhi[0]->GetEntries());
  hmcDecorrPhi[0]->Draw("same");
  legend0->Draw("hist p");

  cQ->cd(2);
  TH1D *hdataDecorrPhi_R= (TH1D*) hdataDecorrPhi[0]->Clone();
  hdataDecorrPhi_R->Divide(hmcDecorrPhi[0]);
  hdataDecorrPhi_R->SetTitle("Data/MC");
  hdataDecorrPhi_R->GetYaxis()->SetRangeUser(0.9,1.1);
  hdataDecorrPhi_R->Draw("");


  //cQ->cd(2);
  //cQ->cd(2)->SetLogz();
  //h2dataPtDecorrPhi[0]->SetTitle("DATA: Jet - lepton decorrelation angle");
  //h2dataPtDecorrPhi[0]->GetXaxis()->SetTitle("(rad)");
  //h2dataPtDecorrPhi[0]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  //h2dataPtDecorrPhi[0]->SetLineWidth(2); //h2dataPtDecorrPhi[0]->SetLineColor(4);
  //h2dataPtDecorrPhi[0]->Draw("colz");
  //h2dataPtDecorrPhi[0]->ProfileX()->Draw("same");

  //cQ->cd(4);
  //cQ->cd(4)->SetLogy();
  //hdataDecorrPhi[0]->SetTitle("Jet - lepton decorrelation angle");
  //hdataDecorrPhi[0]->GetXaxis()->SetTitle("(rad)");
  //hdataDecorrPhi[0]->SetLineWidth(2); hdataDecorrPhi[0]->SetLineColor(4);
  //hdataDecorrPhi[0]->SetMinimum(10);
  //hdataDecorrPhi[0]->Draw();
  //TLegend *legend1 = new TLegend(0.2,0.7,0.35,0.88);
  //legend1->AddEntry(hdataDecorrPhi[0], "jets #geq 1", "l");
  //for(Int_t ijet=1; ijet < 5 ; ijet++){
    //legend1->AddEntry(hdataDecorrPhi[ijet], Form("jets #geq %d",ijet+1), "l");
    //hdataDecorrPhi[ijet]->SetLineWidth(2); hdataDecorrPhi[ijet]->SetLineColor(4+ijet);
    //hdataDecorrPhi[ijet]->Draw("same");
  //}
  //legend1->SetBorderSize(0);
  //legend1->SetFillStyle(0);
  //legend1->Draw("");

  //////ENDJets Plots <--5--

  if (out) {
    //assert(system("mkdir -vp /tmp/ControlPlots") == 0);
    //c0->Print("/tmp/ControlPlots/c0.png");
    //c1->Print("/tmp/ControlPlots/c1.png");
    //ce->Print("/tmp/ControlPlots/ce.png");
    //cQ->Print("/tmp/ControlPlots/cQ.png");
    //assert(system(Form("convert /tmp/ControlPlots/c0.png /tmp/ControlPlots/c1.png /tmp/ControlPlots/ce.png /tmp/ControlPlots/cQ.png %s", out)) == 0);
    assert(system(Form("mkdir -vp %s", out)) == 0);
    c0->Print(Form("%s/c0.png", out));
    c1->Print(Form("%s/c1.png", out));
    ce->Print(Form("%s/ce.png", out));
    cQ->Print(Form("%s/cQ.png", out));
  }
  TDatime now1;
  now1.Print();


}

int
main(int argc, char **argv)
{
  struct option long_opts[] = {
    { "help", false, NULL, 'h' },
    { "data", true, NULL, 'd' },
    { "mc", true, NULL, 'm' },
    { "output", true, NULL, 'o' },
    { "bin", true, NULL, 'b' },
    { "njet", true, NULL, 'j' },
    { "batch-mode", false, NULL, 0x01FF },
    { 0, 0, 0, 0 }
  };

  std::string dpath, mpath, opath;
  bool batch_mode = false;

  int opt;
  while ((opt = getopt_long(argc, argv, "hd:m:o:b:j:", long_opts, NULL)) > 0) {
    switch (opt) {
      case 'h':
        std::cout << "usage: " << argv[0] << " [OPTIONS]" << std::endl
                  << "OPTIONS:" << std::endl
                  << "  --help -h         Show this message and exit." << std::endl
                  << "  --data -d <path>  Specify path to the data." << std::endl
                  << "  --mc   -m <path>  Specify path to the MC (detector level)." << std::endl
                  << "  --bin  -b <#bin>  Number of the Pt/Q2-bin to use." << std::endl
                  << "  --njet -j <#jet>  Multiplicity to use (specifies data and MC files)." << std::endl
                  << "  --batch-mode      Run in batch mode." << std::endl
                  ;
        return EXIT_SUCCESS;

      case 'd':
        dpath = optarg;
        break;

      case 'm':
        mpath = optarg;
        break;

      case 'o':
        opath = optarg;
        break;

      case 'b':
      {
        int bin = atoi(optarg);
        if (bin < 1 || bin > 3) {
          std::cerr << "error: invalid bin number (can only be 1, 2, or 3)" << std::endl;
          return EXIT_FAILURE;
        }
        g_bin1 = g_bin2 = bin;
      }
      break;

      case 'j':
      {
        int n = atoi(optarg);
        if (n < 1 || n > 4) {
          std::cerr << "error: invalid multiplicity, " << n << " (can only be 1, 2, 3, or 4)" << std::endl;
          return EXIT_FAILURE;
        }
        dpath = Form("./Output/abs_eta_le_1/central/pt-mult-bins/data-%djet.root", n);
        mpath = Form("./Output/abs_eta_le_1/central/pt-mult-bins/mcincl-%djet.root", n);
      }
      break;

      case 0x01FF:
        batch_mode = true;
        break;

      default:
        std::cerr << "error: undefined command-line option" << std::endl;
        return EXIT_FAILURE;
    }
  }

  std::clog << "data: " << '"' << dpath << '"' << std::endl;
  std::clog << "MC:   " << '"' << mpath << '"' << std::endl;

  TApplication app { "ControlPlots", &argc, argv };

  ControlPlots(dpath.c_str(), mpath.c_str(), opath.empty() ? NULL : opath.c_str());

  if (not batch_mode)
    app.Run();

  return EXIT_SUCCESS;
}
