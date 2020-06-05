/////Amilkar Quintero
////Temple University
////October 2018
//
//v2, 2 Nov 2018
//Fill the histograms per event to apply the weight
//
//v3, 5 Nov 2018
//- Better way to apply the weight (No weight applied now)
//- New binning and ranges 
//- Include Electron phi histogram
//- Include decorrelation histograms
//
//v4, 10 Nov 2018
//- Loop over each event
//- Modified the orange.h (done with MakeClass()) to accomodate data and MC
//
//v5, 26 Nov 2018
//- Change names of orange to JetDecorr2018
//- Use the cuts suggested in the email of 13 November 2018
//
//v5_1, 11 Dec 2018
//- Changes due comments after the meeting today 
//- Apply "chimney" cut, -1.5 < eta < 1.8, yJB>0.04 and Pt/sqrt(Et)<2.5
//- Remove Q2 high limit.Sinistra might only go to 10000
//- Calculate E-pz from Zufo.
//
//v6, 10 June 2019
//- Cuts used for preliminary
//- Include and modify some histograms
//- Include weights calculated with phi angles
//
//v7, 19 June 2019 (after first prel. presentation)
//- Remove Reweighting
//- y_el < 0.7
//- Include 2D histograms for DecorrPhi detector vs hadron jet, 
//  and detector - hadron jet vs Q2 and Pt
//- Use highest E_T jet for true MC
//
//v8, 21 August 2019 (looking for cross section preliminary)
//- Clean several unused histograms
//- Apply matching of detector and hadron level tracks
//

#define JetOrange2018_cxx
#include <stdlib.h>
#include <string>
#include "Riostream.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
#include "TDatime.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "JetOrange2018.h"

using namespace std;


//--->START MAIN PROGRAM
//________________________________________________________________________________
void MakeHist_v8(
		const Char_t *eachfile = "/pnfs/desy.de/dphep/online/zeus/z/ntup/07p/v08b/mc/root/zeusmc.hfix627.h14207.djangoh.cdm.pos.q24.gnor2515.root",
    const Char_t *outdir="./")
{

  //gROOT->cd();
  TDatime now;                                          //Set time in Root
  now.Print();
  //gErrorIgnoreLevel=kError;       //This gets rid of the warnings
  gErrorIgnoreLevel=kFatal;       //This gets rid of the errors
  //There are some unused branched for data or MC

  /////////////////////Read all files --1-->
  TFile f(eachfile); if (f.IsZombie())  {cout<< "Problems with file: " << eachfile << endl; return;}
  string eachfile2 = eachfile;
  size_t found = eachfile2.find("/data_");
  Bool_t isdata = found != std::string::npos ? true : false;
  Int_t runnumber = 0;
  Char_t runnumber_size[10];
  Char_t sample[50];
  Float_t Weight = 1.0;
  if(isdata){
    cout << "Data" << endl;
    sscanf(&eachfile2[found],"/data_%[^_]_%d_%[^.].root",sample,&runnumber,runnumber_size);
    cout<< "Sample: " << sample << endl;
    cout << "Run Number: " << runnumber << endl;
    //cout<< size << endl;
  } else{
    cout << "MC" << endl;
    found = eachfile2.find("/zeusmc.");
    if (found == std::string::npos) {cout << "Problem finding the MC file" << endl; return;}
    sscanf(&eachfile2[found],"/zeusmc.%s",sample);
    cout << "Sample MC: " << sample << endl;
  }

  TTree *firstJet = (TTree*)f.Get("orange");
  if (!firstJet) {cout<< "Problems with file: " << eachfile << endl; return;}
  JetOrange2018 *JetOrange = new JetOrange2018();
  JetOrange->Init(firstJet);
  JetOrange->SetData(isdata);

  TString period;
  if(runnumber>=47010 && runnumber <=51245)  period = "04p";
  if(runnumber>=52258 && runnumber <=57123)  period = "0405e";
  if(runnumber>=58207 && runnumber <=59947)  period = "06e";
  if(runnumber>=60005 && runnumber <= 62639) period = "0607p";
  if(!isdata) period=sample;
  /////END Reading all files <--1--


  /////////////////////Define Histograms --2-->  
  Int_t q2start = 10;
  Int_t q2end   = 350;
  //Event quantities
  TH1D* hVertexX   = new TH1D("hVertexX", "Vertex X"            , 4,-4,200);
  TH1D* hVertexY   = new TH1D("hVertexY", "Vertex Y"            , 4,-4,200);
  TH1D* hVertexXmY = new TH1D("hVertexXmY", "Vertex X - Y"      , 4,-4,200);
  TH1D* hVertexZ   = new TH1D("hVertexZ", "Vertex Z"            , 40,-40,40);

  TH2D *h2Q2x      = new TH2D("h2Q2x"   , "Q^{2} vs x"          , 100, 2.e-5, 2.e-2, q2end-q2start, q2start, q2end);
  TH1D* hQ2        = new TH1D("hQ2"     , "Q^{2}"               , q2end-q2start, q2start, q2end);
  TH1D* hx         = new TH1D("hx"      , "Momentum fraction x" , 100, 2.e-5, 2.e-2);
  TH1D* hJetMult   = new TH1D("hJetMult", "Jets Multiplicity "  , 11,0,11);
  TH1D* hEmpz      = new TH1D("hEmpz"   , "E_{M} - p_{Z}"       , 80,30,70);
  TH1D* hCalEmpz   = new TH1D("hCalEmpz", "E_{M} - p_{Z} (cal)" , 80,30,70);
  TH1D* hGamma     = new TH1D("hGamma"  , "cos(#gamma_{h})"     , 120,-1,-0.6);
  TH1D* hPtEt      = new TH1D("hPtEt"   , "P_{T} / #sqrt(E_{T}): from Zufo", 90,0,3);
  TH1D* hDiffEmpz  = new TH1D("hDiffEmpz", "E - p_{Z}: Zufo - Cal", 200,-20,20);
  TH1D* hEtamax    = new TH1D("hEtamax", "Eta max of all cell", 200, -2, 5);
  TH1D* hEtamax400 = new TH1D("hEtamax400", "Eta max of all cell > 400 GeV", 200,-20,20);
  TH1D* hEtamax_zu = new TH1D("hEtamax_zu", "Eta max from ZUFOs", 200,-20,20);
  TH1D* hEtamax_zu4 = new TH1D("hEtamax_zu4", "Eta max from ZUFOs (E > 400 MeV)", 200,-20,20);
  //Electron quantities
  TH1D* hElecTheta = new TH1D("hElecTheta", "Angle #theta of the DIS electron", 120, 130., 190.);
  TH1D* hElecPhi   = new TH1D("hElecPhi"  , "Angle #phi of the DIS electron"  , 60,-TMath::Pi(),TMath::Pi());
  TH1D* hElecE     = new TH1D("hElecE"    , "Energy of the DIS electron"      , 150, 10., 25.); //60.);
  TH1D* hElecProb  = new TH1D("hElecProb" , "Sinistra probability"            , 200, 0.89, 1.);
  TH1D* hElecy     = new TH1D("hElecy"    , "Inelasticity JB method"          , 220, -0.1, 1.);
  TH2D* h2ElecPos  = new TH2D("h2ElecPos" , "Sinistra probability"            , 800,-200,200,800,-200,200);
  //Jets
  TH1D* hJetEt   = new TH1D("hJetEt"  ,"Jets Transverse Energy"  , 110,0,55);
  TH1D* hJetMass = new TH1D("hJetMass","Jets Mass"               , 110,0,55);
  TH1D* hJetPt   = new TH1D("hJetPt"  ,"Jets Transverse Momentum", 110,0,55);
  TH1D* hJetEta  = new TH1D("hJetEta" ,"Jets Eta"                , 100,-2.5,2.5);
  TH1D* hJetPhi  = new TH1D("hJetPhi" ,"Jets Phi"                , 60,0,2*TMath::Pi());
  // Resolution
  TH1D* hElecEDiff = new TH1D("hElecEDiff", "Resolution of the electron energy", 150, -10, 10); //60.);
  TH1D* hJetEtDiff = new TH1D("hJetEtDiff", "Resolution of the Jet Transverse Energy", 150, -20, 20);
  TH1D* hPtEtDiff  = new TH1D("hPtEtDiff" , "Resolution of P_{T} / #sqrt(E_{T})", 150, -2, 2);
  TH1D* hElecYjbDiff = new TH1D("hElecYjbDiff", "Resolution of Y_jb", 200, -0.5, 0.5);
  TH1D* hElecYelDiff = new TH1D("hElecYelDiff", "Resolution of Y_el", 200, -0.5, 0.5);
  //Decorrelation
  TH1D* hDecorrPhi[5];   TH2D* h2PtDecorrPhi[5]; TH2D* h2Q2DecorrPhi[5];
  TH2D* h2PtDecorrPhi_true[5]; TH2D* h2Q2DecorrPhi_true[5]; TH2D* h2DecorrPhiHadDec[5]; //TH2D* h2DecorrPhiHadDecMat[5];
  TH3D* h3DecorrPhiHadDecPt[5];  TH3D* h3DecorrPhiHadDecQ2[5];
  TH3D* h3DecorrPhiHadDecPtMat[5];  TH3D* h3DecorrPhiHadDecQ2Mat[5];
  TH2D* h2ElecPhiHadDec[5];
  TH2D* h2JetPhiHadDec[5];
  TH2D* h2JetPtHadDec[5];
  TH2D* h2Q2HadDec[5];
  for(Int_t ijet =0; ijet < 5; ijet++){ 
    hDecorrPhi[ijet]    = new TH1D(Form("hDecorrPhi_%d"   ,ijet),
        Form("Jet - electron decorrelation angle, %d jets",ijet+1), 180,0,TMath::Pi());
    h2PtDecorrPhi[ijet] = new TH2D(Form("h2PtDecorrPhi_%d",ijet),
        Form("P_{T} of the Jet vs angle, %d jets",ijet+1), 180,0,TMath::Pi(),110,0,55);
    h2Q2DecorrPhi[ijet] = new TH2D(Form("h2Q2DecorrPhi_%d",ijet),
        Form("Q^{2} vs angle, %d jets",ijet+1), 180,0,TMath::Pi(),q2end-q2start, q2start, q2end);
    h2PtDecorrPhi_true[ijet] = new TH2D(Form("h2PtDecorrPhi_true_%d",ijet),
        Form("P_{T} of the Jet vs angle, %d jets",ijet+1), 30,TMath::Pi()/2,TMath::Pi(),110,0,55);
    h2Q2DecorrPhi_true[ijet] = new TH2D(Form("h2Q2DecorrPhi_true_%d",ijet),
        Form("Q^{2} vs angle, %d jets",ijet+1), 30,TMath::Pi()/2,TMath::Pi(),q2end-q2start, q2start, q2end);
    h3DecorrPhiHadDecPtMat[ijet] = new TH3D(Form("h3DecorrPhiHadDecPtMat_%d",ijet),
        Form("Matched Hadron vs detector #Delta#phi, %d jets",ijet+1), 180,0,TMath::Pi(), 60,0,TMath::Pi(),110,0,55);
    h3DecorrPhiHadDecQ2Mat[ijet] = new TH3D(Form("h3DecorrPhiHadDecQ2Mat_%d",ijet),
        Form("Matched Hadron vs detector #Delta#phi, %d jets",ijet+1), 180,0,TMath::Pi(), 60,0,TMath::Pi(),34,q2start,q2end);
    h2DecorrPhiHadDec[ijet]    = new TH2D(Form("h2DecorrPhiHadDec_%d",ijet),
        Form("Hadron vs detector #Delta#phi, %d jets",ijet+1), 180,0,TMath::Pi(), 60,0,TMath::Pi());  
    h3DecorrPhiHadDecPt[ijet]  = new TH3D(Form("h3DecorrPhiHadDecPt_%d",ijet),
        Form("Hadron vs detector #Delta#phi, %d jets",ijet+1), 180,0,TMath::Pi(), 60,0,TMath::Pi(),110,0,55);  
    h3DecorrPhiHadDecQ2[ijet]  = new TH3D(Form("h3DecorrPhiHadDecQ2_%d",ijet),
        Form("Hadron vs detector #Delta#phi, %d jets",ijet+1), 180,0,TMath::Pi(), 60,0,TMath::Pi(),34,q2start,q2end);  
    h2ElecPhiHadDec[ijet]   = new TH2D(Form("h2ElecPhiHadDec_%d",ijet),"Hadron vs detector Electron_{#phi}", 60,0,2*TMath::Pi(), 60,0,2*TMath::Pi());
    h2JetPhiHadDec[ijet]    = new TH2D(Form("h2JetPhiHadDec_%d" ,ijet),"Hadron vs detector Jet_{#phi}"     , 60,0,2*TMath::Pi(), 60,0,2*TMath::Pi());
    h2JetPtHadDec[ijet]     = new TH2D(Form("h2JetPtHadDec_%d"  ,ijet),"Hadron vs detector Jet p_{T}"      ,110,0,55,110,0,55);
    h2Q2HadDec[ijet]        = new TH2D(Form("h2Q2HadDec_%d"     ,ijet),"Q2 Hadron vs detector" ,  q2end-q2start, q2start, q2end,  q2end-q2start, q2start, q2end);  
  }
  TH2D* h2JetElecPhi  = new TH2D("h2JetElecPhi" ,"Jet vs electron #phi angle", 60,-TMath::Pi(),TMath::Pi(), 60,0,2*TMath::Pi());
  ////

  ////END Define Histograms <--2--


  ///////Start filling the histograms from the trees  --3-->
  Int_t count = 0;
  Long64_t nentries = firstJet->GetEntriesFast();
  cout << "Number of events before the cuts: " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    JetOrange->GetEntry(jentry);

    /////////////////////Define Cuts --3.1-->
    /////Using only the electron with the highest probability Siq2el[0]
    /////Using only the highest jet in the event: Kt_etajet_b[0]
    /////Phase Space 
    if(JetOrange->Siq2el[0]<q2start || JetOrange->Siq2el[0]>q2end) continue;
    /////Cleanning cuts
    if(JetOrange->Zvtx<-40 || JetOrange->Zvtx>40 || JetOrange->Zvtx==0)    continue;
    if(JetOrange->Cal_empz<45 || JetOrange->Cal_empz>65) continue;  
    //if(JetOrange->Cal_empz<45+3.456 || JetOrange->Cal_empz>65+3.456) continue; // widereta
    //if(JetOrange->Cal_empz<45-3.45 || JetOrange->Cal_empz>65-3.45) continue;   // widereta
    //Calculate E-pz from zufo
    Float_t Empz = 0;
    for(Int_t zloop=0; zloop<JetOrange->Nzufos; zloop++){
      TLorentzVector v(JetOrange->Zufo[zloop][0], JetOrange->Zufo[zloop][1], 
          JetOrange->Zufo[zloop][2], JetOrange->Zufo[zloop][3]);
      Empz += JetOrange->Zufo[zloop][3] - JetOrange->Zufo[zloop][2];
    }
    //  cout << Form("Calorimeter: %f     Zufo: %f     Diff: %f",JetOrange->Cal_empz,Empz,JetOrange->Cal_empz - Empz) << endl;

    /* - - - - - - - - - - - - - - - - - - *
     * res(empz | widereta) = 2.97         *
     * res(empz | abs(eta) < 1) = 2.93     *
     * - - - - - - - - - - - - - - - - - - */
    //Empz += 2.93; // |eta| < 1
    //Empz -= 2.93; // |eta| < 1
    /* - - - - - - - - - - - - - - - - - - */
    //Empz += 2.97; // widereta
    //Empz -= 2.97; // widereta
    /* - - - - - - - - - - - - - - - - - - */
    if(Empz < 45. || Empz > 65.) continue;   
    /* - - - - - - - - - - - - - - - - - - - - - - - - *
     * res(yel | abs(eta) < 1) = 0.027  (gaus fit)     *
     *                         = 0.06 (std dev)        *
     * res(yjb | abs(eta) < 1) = 9.62E-03  (gaus fit)  *
     *                         = 0.04  (std dev)       *
     * - - - - - - - - - - - - - - - - - - - - - - - - */
    double yel = JetOrange->Siyel[0];
    double yjb = JetOrange->Siyjb[0];
    //yel += 0.027;
    //yel -= 0.027;
    //yjb += 9.62e-3;
    //yjb -= 9.62e-3;
    if(yel > 0.7) continue;  //Change from 0.95 to 0.7 for version 7 of this code
    if(yjb < 0.04) continue;
    if(JetOrange->Cal_pt / TMath::Sqrt(JetOrange->Cal_et) > 2.5) continue;
    /////Electron cuts
    if(JetOrange->Siecorr[0][2] < 10) continue;
    if(JetOrange->Sith[0]*180.0/TMath::Pi() < 140 || JetOrange->Sith[0]*180.0/TMath::Pi() > 180.0) continue;  ///Need to reduce 140 to 60 (future)
    if(JetOrange->Sipos[0][2] < -148. && JetOrange->Sipos[0][0] > -14. 
        && JetOrange->Sipos[0][0] < 12. && JetOrange->Sipos[0][1] > 90.) continue;                //Chimney cut
    if(sqrt(JetOrange->Sipos[0][0]*JetOrange->Sipos[0][0] + JetOrange->Sipos[0][1]*JetOrange->Sipos[0][1]) < 20.0) continue;
    if (JetOrange->Sienin[0] > 0.1*(JetOrange->Siein[0] + JetOrange->Sienin[0]) ) continue; // Energy in cone
    if(JetOrange->Siprob[0] < 0.9) continue;  
    //if(JetOrange->Siprob[0] + 2e-3 < 0.9) continue;  
    //if(JetOrange->Siprob[0] - 2e-3 < 0.9) continue;  
    /////Triggers
    if(period == "0405e" && !(JetOrange->Tltw[2] & (1 << 1)) ) continue;      //SPP02
    if(period == "06e"   && !(JetOrange->Tltw[2] & (1 << 8)) ) continue;      //SPP09
    if(period == "0607p" && !(JetOrange->Tltw[2] & (1 << 8)) ) continue;      //SPP09
    /////Jet selection
    if(JetOrange->Kt_njet_b <= 0) continue;

    if(TMath::Abs(JetOrange->Kt_etajet_b[0]) > 1.0 ) continue;
    //if(JetOrange->Kt_etajet_b[0] < -1.5 || JetOrange->Kt_etajet_b[0] > 1.8) continue; // Gmail, Amilkar, 11 Oct. 2019

    /* Centrail Et cut. */
    if(JetOrange->Kt_etjet_b[0] < 2.5) continue;
    /* Cuts for systematics: if Et < 10 then 5% else 2.5% */
    //float jet_Et = JetOrange->Kt_etjet_b[0];
    //if (jet_Et < 10) {
      ////if (jet_Et < 2.37) continue; // left cut
      //if (jet_Et < 2.6) continue; // right cut
    //} else {
      ////if (jet_Et < 2.437) continue; // left cut
      //if (jet_Et < 2.56) continue; // right cut
    //}

    Float_t JetPt = sqrt(JetOrange->Kt_etjet_b[0]*JetOrange->Kt_etjet_b[0]-JetOrange->Kt_masjet_b[0]*JetOrange->Kt_masjet_b[0]) ;
    if(JetPt > 30) continue;
    /////////////////////END Define Cuts <--3.1--


    /////////////////////Fill histograms --3.2-->
    //------------Event
    hVertexX  ->Fill(JetOrange->Xvtx                          ); 
    hVertexY  ->Fill(JetOrange->Yvtx                          ); 
    hVertexXmY->Fill(JetOrange->Xvtx - JetOrange->Yvtx        ); 
    hVertexZ  ->Fill(JetOrange->Zvtx                          ); 

    h2Q2x     ->Fill(JetOrange->Sixel[0],JetOrange->Siq2el[0] );
    hQ2       ->Fill(JetOrange->Siq2el[0]                     );
    hx        ->Fill(JetOrange->Sixel[0]                      );
    hJetMult  ->Fill(JetOrange->Kt_njet_b                     );
    hEmpz     ->Fill(Empz                                     );
    hCalEmpz  ->Fill(JetOrange->Cal_empz                      );
    hGamma    ->Fill(TMath::Cos(JetOrange->Cc_gamma)          );
    hPtEt     ->Fill(JetOrange->Cal_pt/sqrt(JetOrange->Cal_et));
    hDiffEmpz ->Fill(Empz - JetOrange->Cal_empz               );
    hEtamax   ->Fill(JetOrange->Etamax_ce                     );    //Include to check if there are diffractive processes
    hEtamax400->Fill(JetOrange->Etamax_ce4                    );    // Amilkar 3 Sep 2019
    hEtamax_zu->Fill(JetOrange->Etamax_zu);
    hEtamax_zu4->Fill(JetOrange->Etamax_zu4);
    //------------Electron 
    hElecTheta->Fill(JetOrange->Sith[0]*180.0/TMath::Pi()         );
    hElecPhi  ->Fill(JetOrange->Siph[0]                           );
    hElecE    ->Fill(JetOrange->Siecorr[0][2]                     );
    hElecProb ->Fill(JetOrange->Siprob[0]                         );
    hElecy    ->Fill(JetOrange->Siyel[0]                          );
    h2ElecPos ->Fill(JetOrange->Sipos[0][0],JetOrange->Sipos[0][1]);
    //------------Jet
    hJetEt  ->Fill(JetOrange->Kt_etjet_b[0] );
    hJetMass->Fill(JetOrange->Kt_masjet_b[0]);
    hJetPt  ->Fill(JetPt                    );
    hJetEta ->Fill(JetOrange->Kt_etajet_b[0]);
    hJetPhi ->Fill(JetOrange->Kt_phijet_b[0]);

    //------------Correlation
    Float_t ElectronPhi = JetOrange->Siph[0] > 0. ? JetOrange->Siph[0] : 2*TMath::Pi() + JetOrange->Siph[0];  
    Float_t DecorrPhi   = TMath::Abs( JetOrange->Kt_phijet_b[0] - ElectronPhi) ;
    if(DecorrPhi > TMath::Pi()) DecorrPhi = 2*TMath::Pi() - DecorrPhi ;

    for(Int_t ijet =0; ijet < 5; ijet++){ 
      if(JetOrange->Kt_njet_b > ijet) {
        hDecorrPhi[ijet]   ->Fill(DecorrPhi       );
        h2PtDecorrPhi[ijet]->Fill(DecorrPhi, JetPt);
        h2Q2DecorrPhi[ijet]->Fill(DecorrPhi, JetOrange->Siq2el[0]);
      }
    }

    h2JetElecPhi->Fill(JetOrange->Siph[0],JetOrange->Kt_phijet_b[0]);

    ///////////////MC information --3.2.1-->
    if(!isdata){
      Int_t highEt_jet = 0;
      //if(JetOrange->Nhmjets <= 0) continue;
      TLorentzVector jet_initial, jet_compare;
      jet_initial.SetPxPyPzE( 0,0,0,0);
      for(int high = 0; high < JetOrange->Nhmjets; high++){	
        jet_compare.SetPxPyPzE( JetOrange->Pxhmjet[high], JetOrange->Pyhmjet[high], 
            JetOrange->Pzhmjet[high], JetOrange-> Ehmjet[high]);
        if(jet_compare.Et() > jet_initial.Et()) {
          highEt_jet = high;
          jet_initial.SetPxPyPzE(JetOrange->Pxhmjet[high], JetOrange->Pyhmjet[high], 
              JetOrange->Pzhmjet[high], JetOrange-> Ehmjet[high]);
        }
      }
      TLorentzVector jet_hadron(JetOrange->Pxhmjet[highEt_jet], JetOrange->Pyhmjet[highEt_jet], 
          JetOrange->Pzhmjet[highEt_jet], JetOrange-> Ehmjet[highEt_jet]);
      TLorentzVector lep_true(JetOrange->Mc_pfsl[0],JetOrange->Mc_pfsl[1],JetOrange->Mc_pfsl[2],JetOrange->Mc_pfsl[3]);
      Float_t DecorrPhi_true = TMath::Abs(jet_hadron.Phi() - lep_true.Phi());
      if(DecorrPhi_true > TMath::Pi()) DecorrPhi_true = 2*TMath::Pi() - DecorrPhi_true ;
      if(JetOrange->Nhmjets <= 0) continue; //DecorrPhi_true = 0;

      //To match the reco jet with the true jet and reco lepton with true lepton
      Float_t jet_deltaR = sqrt(  (jet_hadron.Eta()-JetOrange->Kt_etajet_b[0])*(jet_hadron.Eta()-JetOrange->Kt_etajet_b[0]) 
          + (jet_hadron.Phi()-JetOrange->Kt_phijet_b[0])*(jet_hadron.Phi()-JetOrange->Kt_phijet_b[0]) );
      Float_t lep_eta = (-1)*TMath::Log( TMath::Tan(JetOrange->Sith[0]/2) );
      Float_t lep_deltaR = sqrt(  (lep_true.Eta()-lep_eta)*(lep_true.Eta()-lep_eta) 
          + (lep_true.Phi()-ElectronPhi)*(lep_true.Phi()-ElectronPhi) );

      // if(jet_deltaR<0.5 && lep_deltaR<0.5){ 
      //cout << lep_deltaR << endl;
      Float_t ElectronPhi_true = lep_true.Phi() > 0. ? lep_true.Phi() : 2*TMath::Pi() + lep_true.Phi();  
      Float_t JetPhi_true = JetOrange->Kt_phijet_b[0] > 0. ? JetOrange->Kt_phijet_b[0] : 2*TMath::Pi() + JetOrange->Kt_phijet_b[0];  
      //}

      for(Int_t ijet =0; ijet < 5; ijet++){ 
        if(JetOrange->Kt_njet_b > ijet) {
          if(jet_deltaR<0.5 && lep_deltaR<0.5){
            h3DecorrPhiHadDecPtMat[ijet]->Fill(DecorrPhi,DecorrPhi_true,JetPt);  
            h3DecorrPhiHadDecQ2Mat[ijet]->Fill(DecorrPhi,DecorrPhi_true,JetOrange->Siq2el[0]);  
          }	  
          h2DecorrPhiHadDec[ijet]->Fill(DecorrPhi,DecorrPhi_true);  
          h3DecorrPhiHadDecPt[ijet]->Fill(DecorrPhi,DecorrPhi_true,JetPt);  
          h3DecorrPhiHadDecQ2[ijet]->Fill(DecorrPhi,DecorrPhi_true,JetOrange->Siq2el[0]);  

          h2PtDecorrPhi_true[ijet]->Fill(DecorrPhi_true, jet_hadron.Pt());
          h2Q2DecorrPhi_true[ijet]->Fill(DecorrPhi_true, JetOrange->Mc_q2);

          h2ElecPhiHadDec[ijet]->Fill(ElectronPhi,ElectronPhi_true);  
          h2JetPhiHadDec[ijet] ->Fill(JetOrange->Kt_phijet_b[0],jet_hadron.Phi());  
          h2JetPtHadDec[ijet]  ->Fill(JetPt,jet_hadron.Pt());  
          h2Q2HadDec[ijet]     ->Fill(JetOrange->Siq2el[0], JetOrange->Mc_q2);
        }
      }

      hElecEDiff->Fill(lep_true.E() - JetOrange->Siecorr[0][2]);
      hJetEtDiff->Fill(jet_hadron.Et() - JetOrange->Kt_etjet_b[0]);
      double PtEt_det = JetOrange->Cal_pt / sqrt(JetOrange->Cal_et);
      double PtEt_had = lep_true.Pt() / sqrt(lep_true.Et());
      hPtEtDiff->Fill(PtEt_det - PtEt_had);


      double y_jb = JetOrange->Siyjb[0];
      double y_el = JetOrange->Siyel[0];
      double y_mc = JetOrange->Mc_y;

      double dy_jb = 0.04 * 0.15;
      double dy_el = 0.7  * 0.15;

      if (0.04 - dy_jb < y_jb && y_jb < 0.04 + dy_jb)
        hElecYjbDiff->Fill(y_jb - y_mc);

      if (0.7 - dy_el < y_el && y_el < 0.7 + dy_el)
        hElecYelDiff->Fill(y_el - y_mc);
    }
    ///////////////END MC information <--3.2.1--


    /////////////////////END Fill histograms <--3.2--
    count++;
  }
  cout << "Number of events after the cuts : " << count << endl;
  ///END Start filling the histograms from the trees  <--3--


  /////////////////////Save histograms in a file --5-->    
  if(hVertexZ->GetEntries() == 0) {cout << "No events in: " << runnumber_size << endl; return ;}
  TString outdir1(outdir);
  TString outfile(outdir1+"/Hist_"+period+"_"+runnumber_size+".root");
  TString outfileMC(outdir1+"/Hist_"+sample);
  if(Weight!=1) outfileMC(outdir1+"/Hist_rw"+sample); 
  TFile fout(isdata ? outfile : outfileMC,"RECREATE");
  //Event
  hVertexX->Write();
  hVertexY->Write();
  hVertexXmY->Write();
  hVertexZ->Write(); //h2Q2x->Write(); 
  hQ2->Write(); hx->Write(); 
  hJetMult->Write(); 
  hCalEmpz->Write();
  hEmpz->Write(); hGamma->Write();  hPtEt->Write(); hDiffEmpz->Write();
  hEtamax->Write();
  hEtamax400->Write();
  hEtamax_zu->Write();
  hEtamax_zu4->Write();
  //Electron
  hElecTheta->Write(); hElecPhi->Write(); hElecE->Write(); hElecProb->Write(); 
  hElecy->Write(); h2ElecPos->Write();
  //Jet
  hJetEt->Write(); hJetPt->Write(); hJetMass->Write(); hJetEta->Write(); hJetPhi->Write();
  // Resolution
  hElecEDiff->Write();
  hJetEtDiff->Write();
  hPtEtDiff->Write();
  hElecYjbDiff->Write();
  hElecYelDiff->Write();
  //Decorrelation
  for(Int_t ijet =0; ijet < 5; ijet++){ 
    //hDecorrPhi[ijet]->Write();   
    h2PtDecorrPhi[ijet]->Write();
    h2Q2DecorrPhi[ijet]->Write();
    //h2DecorrPhiHadDec[ijet]->Write(); 
    if(!isdata) {   
      h2PtDecorrPhi_true[ijet]->Write();h2Q2DecorrPhi_true[ijet]->Write();
      //h2DecorrPhiHadDecMat[ijet]->Write(); 
      h3DecorrPhiHadDecPt[ijet]->Write(); h3DecorrPhiHadDecQ2[ijet]->Write(); 
      h3DecorrPhiHadDecPtMat[ijet]->Write(); h3DecorrPhiHadDecQ2Mat[ijet]->Write(); 
      h2ElecPhiHadDec[ijet]->Write();   h2JetPhiHadDec[ijet]->Write();   
      h2JetPtHadDec[ijet]->Write(); h2Q2HadDec[ijet]->Write();
    }
  }
  //h2JetElecPhi->Write();
  fout.Close();
  /////END Save histograms in a file <--5-- 


  TDatime now1;
  now1.Print();
}



// A function to implement bubble sort 
/*void bubbleSort(int arr[], int n) { 
  int i, j; 
  for (i = 0; i < n-1; i++)       

// Last i elements are already in place    
for (j = 0; j < n-i-1; j++)  
if (arr[j] > arr[j+1]) 
swap(&arr[j], &arr[j+1]); 
} 
*/
