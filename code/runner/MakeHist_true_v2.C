/////Amilkar Quintero
////Temple University
////June 2019
//
// True MC data only
//
//v2, 8 July 2019
//Use the decorrelation angle of the leading E_T jet and the lepton
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
void MakeHist_true_v2(const Char_t *eachfile= "~/Desktop/zeusmc.hfix627.h1391.0607p.q4.ari_2911.root ", 
    const Char_t *outdir="./"){

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
    //Weight = CalcWeight(sample);
    //cout << "Weight   : " << Weight << endl; 
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
  TH1D* hVertexZ   = new TH1D("hVertexZ", "Vertex Z"            , 40,-40,40);
  TH2D *h2Q2x      = new TH2D("h2Q2x"   , "Q^{2} vs x"          , 100, 2.e-5, 2.e-2, q2end-q2start, q2start, q2end);
  TH1D* hQ2        = new TH1D("hQ2"     , "Q^{2}"               , q2end-q2start, q2start, q2end);
  TH1D* hx         = new TH1D("hx"      , "Momentum fraction x" , 100, 2.e-5, 2.e-2);
  TH1D* hJetMult   = new TH1D("hJetMult", "Jets Multiplicity "  , 11,0,11);
  TH1D* hEmpz      = new TH1D("hEmpz"   , "E_{M} - p_{Z}"       , 80,30,70);
  TH1D* hGamma     = new TH1D("hGamma"  , "cos(#gamma_{h})"     , 120,-1,-0.6);
  TH1D* hPtEt      = new TH1D("hPtEt"   , "P_{T} / #sqrt(E_{T}): from Zufo", 90,0,3);
  TH1D* hDiffEmpz  = new TH1D("hDiffEmpz", "E - p_{Z}: Zufo - Cal", 200,-20,20);
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
  //Decorrelation
  TH1D* hDecorrPhi[5];   TH2D* h2PtDecorrPhi[5]; TH2D* h2ElDecorrPhi[5]; TH2D* h2Q2DecorrPhi[5];
  for(Int_t ijet =0; ijet < 5; ijet++){ 
    hDecorrPhi[ijet]    = new TH1D(Form("hDecorrPhi_%d"   ,ijet),
        Form("Jet - electron decorrelation angle, %d jets",ijet+1), 60,0,TMath::Pi());
    h2PtDecorrPhi[ijet] = new TH2D(Form("h2PtDecorrPhi_%d",ijet),
        Form("P_{T} of the Jet vs angle, %d jets",ijet+1), 60,0,TMath::Pi(),110,0,55);
    h2ElDecorrPhi[ijet] = new TH2D(Form("h2ElDecorrPhi_%d",ijet),
        Form("Energy of the Lepton vs angle, %d jets",ijet+1), 60,0,TMath::Pi(),150,10,25);
    h2Q2DecorrPhi[ijet] = new TH2D(Form("h2Q2DecorrPhi_%d",ijet),
        Form("Q^{2} vs angle, %d jets",ijet+1), 60,0,TMath::Pi(),q2end-q2start, q2start, q2end);
  }
  TH2D* h2JetElecPhi  = new TH2D("h2JetElecPhi" ,"Jet vs electron #phi angle", 60,-TMath::Pi(),TMath::Pi(), 60,0,2*TMath::Pi());

  TH2D* h2JetPt     = new TH2D("h2JetPt" ,"Jet p_{T} true vs MC", 110,0,55, 110,0,55);
  TH2D* h2ElecE  = new TH2D("h2ElecE" ,"Jet energy true vs MC", 150, 10., 25., 150, 10., 25.);
  TH2D* h2ElecPhi   = new TH2D("h2ElecPhi"  , "Angle #phi of the DIS electron true vs MC"  , 60,-TMath::Pi(),TMath::Pi() , 60,-TMath::Pi(),TMath::Pi());
  TH2D* h2JetPhi  = new TH2D("h2JetPhi" ,"Jets Phi true vs MC"     , 60,0,2*TMath::Pi() , 60,0,2*TMath::Pi());
  ////END Define Histograms <--2--


  ///////Start filling the histograms from the trees  --3-->
  Int_t count = 0, nothighEt=0;
  Long64_t nentries = firstJet->GetEntriesFast();
  cout << "Number of events before the cuts: " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    JetOrange->GetEntry(jentry);

    /////////////////////Define Cuts --3.1-->
    /////Using only the electron with the highest probability Siq2el[0]
    /////Using only the highest jet in the event: Kt_etajet_b[0]
    /////Phase Space 
    if(JetOrange->Mc_q2<q2start || JetOrange->Mc_q2>q2end) continue;
    /////Cleanning cuts
    //if(JetOrange->Mc_vtx[2]<-40 || JetOrange->Mc_vtx[2]>40 || JetOrange->Zvtx==0)continue;   //Vertex is detector component Achim 19/Jun/2019
    //if(JetOrange->Cal_empz<35 || JetOrange->Cal_empz>65) continue;  
    //Calculate E-pz from zufo
    Float_t Empz = JetOrange->Cal_empz;
    //for(Int_t zloop=0; zloop<JetOrange->Nzufos; zloop++){
    // TLorentzVector v(JetOrange->Zufo[zloop][0], JetOrange->Zufo[zloop][1], 
    // JetOrange->Zufo[zloop][2], JetOrange->Zufo[zloop][3]);
    // Empz += JetOrange->Zufo[zloop][3] - JetOrange->Zufo[zloop][2];
    // }
    //  cout << Form("Calorimeter: %f     Zufo: %f     Diff: %f",JetOrange->Cal_empz,Empz,JetOrange->Cal_empz - Empz) << endl;
    //  if(Empz < 35. || Empz > 65.) continue;   
    if(JetOrange->Mc_y > 0.7)  continue;   //Change from 0.95 to 0.7 , first prel. meeting
    if(JetOrange->Mc_y < 0.04) continue;
    //if(JetOrange->Cal_pt / TMath::Sqrt(JetOrange->Cal_et) > 2.5) continue;
    /////Electron cuts
    TLorentzVector lepton_final;
    lepton_final.SetPxPyPzE( JetOrange->Mc_pfsl[0], JetOrange->Mc_pfsl[1], JetOrange->Mc_pfsl[2], JetOrange->Mc_pfsl[3]);
    if(lepton_final.E() < 10) continue;  //Remove the upper energy cut
    if(lepton_final.Theta()*180.0/TMath::Pi() < 140 || lepton_final.Theta()*180.0/TMath::Pi() > 180.0) continue;
    //if(JetOrange->Sipos[0][2] < -148. && JetOrange->Sipos[0][0] > -14. 
    //   && JetOrange->Sipos[0][0] < 12. && JetOrange->Sipos[0][1] > 90.) continue;                //Chimney cut
    //if(sqrt(JetOrange->Sipos[0][0]*JetOrange->Sipos[0][0] + JetOrange->Sipos[0][1]*JetOrange->Sipos[0][1]) < 20.0) continue;
    //if (JetOrange->Sienin[0] > 0.1*(JetOrange->Siein[0] + JetOrange->Sienin[0]) ) continue; // Energy in cone                                                                                        
    /////Triggers
    //if(period == "0405e" && !(JetOrange->Tltw[2] & (1 << 1)) ) continue;      //SPP02
    //if(period == "06e"   && !(JetOrange->Tltw[2] & (1 << 8)) ) continue;      //SPP09
    //if(period == "0607p" && !(JetOrange->Tltw[2] & (1 << 8)) ) continue;      //SPP09

    /////Jet selection
    Int_t highEt_jet = 0;
    if(JetOrange->Nhmjets <= 0) continue;
    TLorentzVector jet_initial, jet_compare, jet_final;
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
    jet_final.SetPxPyPzE( JetOrange->Pxhmjet[highEt_jet], JetOrange->Pyhmjet[highEt_jet], 
        JetOrange->Pzhmjet[highEt_jet], JetOrange-> Ehmjet[highEt_jet]);
    //cout << "The maximum E_T jet is: " << highEt_jet << endl;
    if(highEt_jet != 0) nothighEt++;
    if(TMath::Abs(jet_final.Eta()) > 1.0 ) continue;
    //if(JetOrange->Kt_etajet_b[0]<-1.5 || JetOrange->Kt_etajet_b[0]>1.8) continue;
    if(jet_final.Et() < 2.5) continue;
    if(jet_final.Pt() > 30) continue;
    /////////////////////END Define Cuts <--3.1--


    /////////////////////Fill histograms --3.2-->
    //if(!isdata) Weight = CalcWeight(JetOrange->Siph[0],JetOrange->Kt_phijet_b[0]);
    //cout << "Weight   : " << Weight << endl; 
    //if(Weight ==0) continue; 
    //------------Event
    hVertexZ->Fill(JetOrange->Mc_vtx[2]                         ,Weight); 
    h2Q2x   ->Fill(JetOrange->Mc_x,JetOrange->Mc_q2 ,Weight);
    hQ2     ->Fill(JetOrange->Mc_q2                     ,Weight);
    hx      ->Fill(JetOrange->Mc_x                      ,Weight);
    hJetMult->Fill(JetOrange->Nhmjets                     ,Weight);
    //    hEmpz   ->Fill(Empz                                     ,Weight);
    //hEmpz   ->Fill(JetOrange->Cal_empz                      ,Weight);
    //hGamma  ->Fill(TMath::Cos(JetOrange->Cc_gamma)          ,Weight);
    //hPtEt   ->Fill(JetOrange->Cal_pt/sqrt(JetOrange->Cal_et),Weight);
    //hDiffEmpz->Fill(Empz - JetOrange->Cal_empz              ,Weight);
    //------------Electron 
    hElecTheta->Fill(lepton_final.Theta()*180.0/TMath::Pi()         ,Weight);
    hElecPhi  ->Fill(lepton_final.Phi()                          ,Weight);
    hElecE    ->Fill(lepton_final.E()                     ,Weight);
    //hElecProb ->Fill(JetOrange->Siprob[0]                         ,Weight);
    hElecy    ->Fill(JetOrange->Mc_y                          ,Weight);
    //h2ElecPos ->Fill(JetOrange->Sipos[0][0],JetOrange->Sipos[0][1],Weight);
    //------------Jet
    hJetEt  ->Fill(jet_final.Et() ,Weight);
    hJetMass->Fill(jet_final.M(),Weight);
    //Float_t JetPt = sqrt(JetOrange->Kt_etjet_b[0]*JetOrange->Kt_etjet_b[0]-JetOrange->Kt_masjet_b[0]*JetOrange->Kt_masjet_b[0]) ;
    hJetPt  ->Fill(jet_final.Pt()                    ,Weight);
    hJetEta ->Fill(jet_final.Eta(),Weight);
    hJetPhi ->Fill(jet_final.Phi(),Weight);

    //------------Correlation
    Float_t ElectronPhi = lepton_final.Phi() > 0. ? lepton_final.Phi() : 2*TMath::Pi() + lepton_final.Phi();  
    Float_t DecorrPhi   = TMath::Abs( jet_final.Phi() - ElectronPhi) ;
    if(DecorrPhi > TMath::Pi()) DecorrPhi = 2*TMath::Pi() - DecorrPhi ;

    for(Int_t ijet =0; ijet < 5; ijet++){ 
      if(JetOrange->Nhmjets > ijet) {
        hDecorrPhi[ijet]   ->Fill(DecorrPhi       ,Weight);
        h2PtDecorrPhi[ijet]->Fill(DecorrPhi, jet_final.Pt(),Weight);
        h2ElDecorrPhi[ijet]->Fill(DecorrPhi, lepton_final.E(),Weight);
        h2Q2DecorrPhi[ijet]->Fill(DecorrPhi, JetOrange->Mc_q2,Weight);
      }
    }

    h2JetElecPhi->Fill(lepton_final.Phi(),jet_final.Phi(),Weight);
    /////////////////////END Fill histograms <--3.2--
    count++;
  }
  cout << "Number of events after the cuts : " << count << endl;
  cout << "Number of events with high Et != index 0 : " << nothighEt << endl;
  ///END Start filling the histograms from the trees  <--3--


  /////////////////////Save histograms in a file --5-->    
  if(hVertexZ->GetEntries() == 0) {cout << "No events in: " << runnumber_size << endl; return 0;}
  TString outdir1(outdir);
  TString outfile(outdir1+"/Hist_"+period+"_"+runnumber_size+".root");
  TString outfileMC(outdir1+"/Hist_true_"+sample);
  //if(Weight!=1) outfileMC = outdir1+"/Hist_rw"+sample; 
  TFile fout(isdata ? outfile : outfileMC,"RECREATE");
  //Event
  hVertexZ->Write(); h2Q2x->Write(); hQ2->Write(); hx->Write(); hJetMult->Write(); 
  //hEmpz->Write(); hGamma->Write();  hPtEt->Write(); hDiffEmpz->Write();
  //Electron
  hElecTheta->Write(); hElecPhi->Write(); hElecE->Write(); hElecy->Write(); //h2ElecPos->Write(); hElecProb->Write();
  //Jet
  hJetEt->Write(); hJetPt->Write(); hJetMass->Write(); hJetEta->Write(); hJetPhi->Write();
  //Decorrelation
  for(Int_t ijet =0; ijet < 5; ijet++){ 
    hDecorrPhi[ijet]->Write();   h2PtDecorrPhi[ijet]->Write();
    //h2ElDecorrPhi[ijet]->Write();
    h2Q2DecorrPhi[ijet]->Write();}
  h2JetElecPhi->Write();
  h2JetPt->Write(); h2ElecE->Write(); h2ElecPhi->Write(); h2JetPhi->Write();
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
if (arr[j] > arr[j+1]) ////CHANGE here E_T[j] > E_T[j+1]
swap(&arr[j], &arr[j+1]); 
} 
*/

