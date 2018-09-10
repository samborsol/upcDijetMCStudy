#include <ctime>
#include <TMath.h>
#include "../cutsAndBin_bgk.h"
#include <TLorentzVector.h>
#include ".../commonUtility.h"
#include "Riostream.h"
#include <iostream>
#include <algorithm>
#include <vector>

using std::stringstream;
using std::vector;
using std::abs;
     
static const long MAXTREESIZE = 10000000000;

bool arrange(float i, float j) {return (i>j);}
TString getDayAndTime();

void forest2diJetSkim_pp5(
                           TString fname = "/u/user/bekim/758p3_UPC/src/Muon_test/4_UPCTriggers_pp_reco_170427.root",
			   TString outputFname = "upcDiJetSkim180827", 
			   TString trig = "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1",  // "HLT_HIUPCL1SingleEG5NotHF2_v1",    // "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1" //  "HLT_HIUPCL1SingleEG5NotHF2_v1"
			   TString jetCollection = "ak4PFJetAnalyzer", // "akPu5PFJetAnalyzer",
			   float minjPt = 0,
			   int nevt=-1
			   ) {  
  
  using namespace std;
  TFile *f1 = new TFile(fname.Data());

  TFile *f2 = new TFile("h_eff.root");
  TH1F *heff = (TH1F*)f2->Get("h4");
  
  TTree *HltTree = (TTree*)f1->Get("hltanalysis/HltTree");
  TTree *hiTree  = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");
//  TTree *trackTree = (TTree*)f1->Get("anaTrack/trackTree");
  TTree *t = (TTree*)f1->Get(Form("%s/t",jetCollection.Data()));
  TTree *trackTree  = (TTree*)f1->Get("ppTrack/trackTree");
  TTree *hbhe = (TTree*)f1->Get("rechitanalyzer/hbhe");
  TTree *hf = (TTree*)f1->Get("rechitanalyzer/hf");
  TTree *ee = (TTree*)f1->Get("rechitanalyzer/ee");
  TTree *eb = (TTree*)f1->Get("rechitanalyzer/eb");
  //  TTree *akpu5pf = (TTree*)f1->Get("t");
  //  TTree *pfTree  = (TTree*)f1->Get("pfcandAnalyzer/pfTree");
  TString dayTime = getDayAndTime();
  TFile* newfile = new TFile(Form("../data/skimmedFiles/TrkCutsppreco_%s_trig%s_jetCollection%s_minJetPt%d_%s.root",outputFname.Data(), trig.Data(), jetCollection.Data(), (int)minjPt, dayTime.Data()  ),"recreate");
//  TFile* newfile = new TFile("skimmedFiles/test.root", "recreate");  

  t->AddFriend(HltTree);
  t->AddFriend(hiTree);
  t->AddFriend(trackTree);
  t->AddFriend(hbhe);
  t->AddFriend(hf);
  t->AddFriend(ee);
  t->AddFriend(eb);

//  hf->AddFriend(t);
//  hf->AddFriend(hiTree);
   // import the tree to the RooDataSet
   Int_t           trigBit;
   TBranch        *b_trigBit;
   if (trig == "" ) {  cout << " No Trigger selection! " << endl ;}     
   else { 
     HltTree->SetBranchAddress(trig.Data(), &trigBit, &b_trigBit);
   }

  Int_t           Run;
  TBranch        *b_Run;   //!
  HltTree->SetBranchAddress("Run", &Run, &b_Run);
  Int_t           LumiBlock;
  TBranch        *b_LumiBlock;   //!
  HltTree->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
  ULong64_t       Event;
  TBranch        *b_Event;   //!
  HltTree->SetBranchAddress("Event", &Event, &b_Event);

  UInt_t          run;
  ULong64_t       evt;
  UInt_t          lumi;
  TBranch        *b_run;
  TBranch        *b_evt;
  TBranch        *b_lumi;
  hiTree->SetBranchAddress("run", &run, &b_run);
  hiTree->SetBranchAddress("evt", &evt, &b_evt);
  hiTree->SetBranchAddress("lumi", &lumi, &b_lumi);

  /// Photon inputs : 
/*  Int_t           nPFpart;
  vector<int>     *pfId;
  vector<float>   *pfPt;
  vector<float>   *pfEnergy;
  vector<float>   *pfVsPtInitial;
  vector<float>   *pfEta;
  vector<float>   *pfPhi;
  nPFpart = 0;
  pfId = 0;
  pfPt = 0;
  pfEnergy = 0;
  pfVsPtInitial = 0;
  pfEta = 0;
  pfPhi = 0;

  TBranch        *b_nPFpart;   //!
  TBranch        *b_pfId;   //!
  TBranch        *b_pfPt;   //!
  TBranch        *b_pfEnergy;   //!
  TBranch        *b_pfVsPtInitial;   //!
  TBranch        *b_pfEta;   //!
  TBranch        *b_pfPhi;   //!
  
   pfTree->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
   pfTree->SetBranchAddress("pfId", &pfId, &b_pfId);
   pfTree->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
   pfTree->SetBranchAddress("pfEnergy", &pfEnergy, &b_pfEnergy);
   pfTree->SetBranchAddress("pfVsPtInitial", &pfVsPtInitial, &b_pfVsPtInitial);
   pfTree->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
   pfTree->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);
*/
  // Jet inputs :
  Int_t nref;
  Float_t jtpt[200] = {0};
  Float_t rawpt[200] = {0};
  Float_t jteta[200] = {0};
  Float_t jty[200] = {0};
  Float_t jtphi[200] = {0};
  Float_t jtm[200] = {0};
  Float_t jtmass[200] = {0};
  Float_t pt[200] = {0};
  Float_t eta[200] = {0};
  Float_t y[200] = {0};
  Float_t phi[200] = {0};
  Float_t m[200] = {0};
//  Float_t e[200] = {0};
  Float_t mass[200] = {0}; 
  Float_t jtpt2[200] = {0};
  Float_t rawpt2[200] = {0};
  Float_t jteta2[200] = {0};
  Float_t jty2[200] = {0};
  Float_t jtphi2[200] = {0};
  Float_t jtm2[200] = {0};
  Float_t pt2[200] = {0};
  Float_t eta2[200] = {0};
  Float_t y2[200] = {0};
  Float_t phi2[200] = {0};
  Float_t m2[200] = {0};
  Float_t e2[200] = {0};
  Float_t mass2[200] = {0};

  TBranch *b_nref;
  TBranch *b_jtpt;
  TBranch *b_rawpt;
  TBranch *b_jteta;
  TBranch *b_jty;
  TBranch *b_jtphi;
  TBranch *b_jtm;

  t->SetBranchAddress("nref", &nref, &b_nref);
  t->SetBranchAddress("jtpt", jtpt, &b_jtpt);
  t->SetBranchAddress("rawpt", rawpt, &b_rawpt);
  t->SetBranchAddress("jteta", jteta, &b_jteta);
  t->SetBranchAddress("jty", jty, &b_jty);
  t->SetBranchAddress("jtphi", jtphi, &b_jtphi);
  t->SetBranchAddress("jtm", jtm, &b_jtm);

   // HiTree inputs : 
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiHFplusEta4;   //!
   TBranch        *b_hiHFminusEta4;   //!
   hiTree->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   hiTree->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   hiTree->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   hiTree->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   hiTree->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   hiTree->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);

   Int_t           nTrk;
   Int_t           nVtx;
   Float_t         xVtx[50];
   Float_t         yVtx[50];
   Float_t         zVtx[50];
   Float_t         trkEta[33885];
   Float_t         trkPhi[33885];
   Float_t         trkPt[33885];
   Float_t         trkPtError[33885];
   Bool_t          highPurity[33885];
   Float_t         trkDxy1[33885];   //[nTrk]
   Float_t         trkDxyError1[33885];   //[nTrk]
   Float_t         trkDz1[33885];   //[nTrk]
   Float_t         trkDzError1[33885];   //[nTrk]
   UChar_t         trkNHit[33885];
   Float_t         trkChi2[33885];
   UChar_t         trkNdof[33885];
   UChar_t         trkNlayer[33885];
   TBranch        *b_nTrk;   //!
   TBranch        *b_nVtx;
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_trkEta;
   TBranch        *b_trkPhi;
   TBranch        *b_trkPt;
   TBranch        *b_trkPtError;
   TBranch        *b_highPurity;
   TBranch        *b_trkDxy1;
   TBranch        *b_trkDxyError1;
   TBranch        *b_trkDz1;
   TBranch        *b_trkDzError1;
   TBranch        *b_trkNHit;
   TBranch        *b_trkChi2;
   TBranch        *b_trkNdof;
   TBranch        *b_trkNlayer;
   trackTree->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   trackTree->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   trackTree->SetBranchAddress("xVtx", xVtx, &b_xVtx);
   trackTree->SetBranchAddress("yVtx", yVtx, &b_yVtx);
   trackTree->SetBranchAddress("zVtx", zVtx, &b_zVtx);
   trackTree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   trackTree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   trackTree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   trackTree->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   trackTree->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   trackTree->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   trackTree->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   trackTree->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   trackTree->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   trackTree->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   trackTree->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   trackTree->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   trackTree->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);

   Int_t           hben;
   Float_t         hbee[1000];   
   Float_t         hbeet[1000];   
   Float_t         hbeeta[1000];   
   Float_t         hbephi[1000];   
   Float_t         hbeperp[1000];   
   Bool_t          hbeisjet[1000];   
   Int_t           hbedepth[1000];   
   TBranch        *b_hben;   
   TBranch        *b_hbee;   
   TBranch        *b_hbeet;   
   TBranch        *b_hbeeta;   
   TBranch        *b_hbephi;   
   TBranch        *b_hbeperp;   
   TBranch        *b_hbeisjet;   
   TBranch        *b_hbedepth;   
   hbhe->SetBranchAddress("n", &hben, &b_hben);
   hbhe->SetBranchAddress("e", hbee, &b_hbee);
   hbhe->SetBranchAddress("et", hbeet, &b_hbeet);
   hbhe->SetBranchAddress("eta", hbeeta, &b_hbeeta);
   hbhe->SetBranchAddress("phi", hbephi, &b_hbephi);
   hbhe->SetBranchAddress("perp", hbeperp, &b_hbeperp);
   hbhe->SetBranchAddress("isjet", hbeisjet, &b_hbeisjet);
   hbhe->SetBranchAddress("depth", hbedepth, &b_hbedepth); 

   Int_t           hfn;
   Float_t         hfe[1000];
   Float_t         hfet[1000];
   Float_t         hfeta[1000];
   Float_t         hfphi[1000];
   Float_t         hfperp[1000];
   Bool_t          hfisjet[1000];
   Int_t           hfdepth[1000];
   Float_t         totalEtHFplus;
   Float_t         totalEtHFminus;
   Float_t         totalEtHF;
   TBranch        *b_hfn;
   TBranch        *b_hfe;
   TBranch        *b_hfet;
   TBranch        *b_hfeta;
   TBranch        *b_hfphi;
   TBranch        *b_hfperp;
   TBranch        *b_hfisjet;
   TBranch        *b_hfdepth;
   TBranch        *b_totalEtHFplus;
   TBranch        *b_totalEtHFminus;
   TBranch        *b_totalEtHF;
   hf->SetBranchAddress("n", &hfn, &b_hfn);
   hf->SetBranchAddress("e", hfe, &b_hfe);
   hf->SetBranchAddress("et", hfet, &b_hfet);
   hf->SetBranchAddress("eta", hfeta, &b_hfeta);
   hf->SetBranchAddress("phi", hfphi, &b_hfphi);
   hf->SetBranchAddress("perp", hfperp, &b_hfperp);
   hf->SetBranchAddress("isjet", hfisjet, &b_hfisjet);
   hf->SetBranchAddress("depth", hfdepth, &b_hfdepth);
   hf->SetBranchAddress("totalEtHFplus",&totalEtHFplus,&b_totalEtHFplus);
   hf->SetBranchAddress("totalEtHFminus",&totalEtHFminus,&b_totalEtHFminus);
   hf->SetBranchAddress("totalEtHF",&totalEtHF,&b_totalEtHF);

   Int_t           een;
   Float_t         eee[1000];
   Float_t         eeet[1000];
   Float_t         eeeta[1000];
   Float_t         eephi[1000];
   Float_t         eeperp[1000];
   Bool_t          eeisjet[1000];
   TBranch        *b_een;
   TBranch        *b_eee;
   TBranch        *b_eeet;
   TBranch        *b_eeeta;
   TBranch        *b_eephi;
   TBranch        *b_eeperp;
   TBranch        *b_eeisjet;
   ee->SetBranchAddress("n", &een, &b_een);
   ee->SetBranchAddress("e", eee, &b_eee);
   ee->SetBranchAddress("et", eeet, &b_eeet);
   ee->SetBranchAddress("eta", eeeta, &b_eeeta);
   ee->SetBranchAddress("phi", eephi, &b_eephi);
   ee->SetBranchAddress("perp", eeperp, &b_eeperp);
   ee->SetBranchAddress("isjet", eeisjet, &b_eeisjet);

   Int_t           ebn;
   Float_t         ebe[1000];
   Float_t         ebet[1000];
   Float_t         ebeta[1000];
   Float_t         ebphi[1000];
   Float_t         ebperp[1000];
   Bool_t          ebisjet[1000];
   TBranch        *b_ebn;
   TBranch        *b_ebe;
   TBranch        *b_ebet;
   TBranch        *b_ebeta;
   TBranch        *b_ebphi;
   TBranch        *b_ebperp;
   TBranch        *b_ebisjet;
   eb->SetBranchAddress("n", &ebn, &b_ebn);
   eb->SetBranchAddress("e", ebe, &b_ebe);
   eb->SetBranchAddress("et", ebet, &b_ebet);
   eb->SetBranchAddress("eta", ebeta, &b_ebeta);
   eb->SetBranchAddress("phi", ebphi, &b_ebphi);
   eb->SetBranchAddress("perp", ebperp, &b_ebperp);
   eb->SetBranchAddress("isjet", ebisjet, &b_ebisjet);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  dijet tree 
  ////////////////////////////////////////////////////////////////////////
   UPCEvent event;
   TTree *eventTree = new TTree("evt","Event Tree");
   eventTree->SetMaxTreeSize(MAXTREESIZE);
   eventTree->Branch("event",&event.run,eventBranchString.Data());
//   eventTree->Branch("

/*   DiPhoPF dp;
   TTree *dpTree = new TTree("diphoton","Photon Tree");
   dpTree->SetMaxTreeSize(MAXTREESIZE);
   dpTree->Branch("dp",&dp.mass,dpBranchStringPF);
*/
   UPCdiJet dj;
   TTree *djTree = new TTree("dijet","Jet Tree");
   djTree->SetMaxTreeSize(MAXTREESIZE);
   djTree->Branch("dj",&dj.nJet,djBranchString);

   UPCnTrk Track;
   TTree *trkTree = new TTree("Track","Track Tree");
   trkTree->SetMaxTreeSize(MAXTREESIZE);
   trkTree->Branch("trkVar",&Track.nTrack,nTrkString);

   int ntrk;
   float floatntrk;
   const int MAXtrk = 50000; // to accomodate 100 smeared trks, need to be careful with ram
   const int na = 100;
   const int na2 = 1000;
   float trkpt[MAXtrk];
   float trketa[MAXtrk];
   float trkphi[MAXtrk];
   float trkHFplus;
   float trkHFminus;
   float FRP;
   float BRP;

   TTree *newtrkTree = new TTree("fullTrkTree","Track Tree 2");
   newtrkTree->SetMaxTreeSize(MAXTREESIZE);
   newtrkTree->Branch("ntrk",&ntrk,"ntrk/I");
   newtrkTree->Branch("floatntrk",&floatntrk,"floatntrk/F");
   newtrkTree->Branch("pT",trkpt,"pT[ntrk]/F");
   newtrkTree->Branch("Eta",trketa,"Eta[ntrk]/F");
   newtrkTree->Branch("Phi",trkphi,"Phi[ntrk]/F");
   newtrkTree->Branch("trkHFplus",&trkHFplus,"trkHFplus/F");
   newtrkTree->Branch("trkHFminus",&trkHFminus,"trkHFminus/F");
   newtrkTree->Branch("FRP",&FRP,"FRP/F");
   newtrkTree->Branch("BRP",&BRP,"BRP/F");

   const float NE = -2.5;
   const float PE = +2.5;

   float floatnTrkabsEtaover1p5;
   float hemax;

   TTree *trkheTree = new TTree("trkandHE","TrackVsHE"); 
   trkheTree->SetMaxTreeSize(MAXTREESIZE);
   trkheTree->Branch("floatnTrkabsEtaover1p5",&floatnTrkabsEtaover1p5,"floatnTrkabsEtaover1p5/F");
   trkheTree->Branch("hemax",&hemax,"hemax/F");

   Int_t           HBHEn;
   Float_t         HBHEe[1000];   
   Float_t         HBHEet[1000];   
   Float_t         HBHEeta[1000];   
   Float_t         HBHEphi[1000];   
   Float_t         HBHEperp[1000];   
   Bool_t          HBHEisjet[1000];   
   Int_t           HBHEdepth[1000];   
   Float_t         HBmax;
   Float_t         HEmax;
   Int_t           HFn;
   Float_t         HFe[1000];
   Float_t         HFet[1000];
   Float_t         HFeta[1000];
   Float_t         HFphi[1000];
   Float_t         HFperp[1000];
   Bool_t          HFisjet[1000];
   Int_t           HFdepth[1000];
   Float_t         HFtotal;
   Float_t         HFplus;
   Float_t         HFminus;
   Float_t         HFplusmax;
   Float_t         HFminusmax;
   Int_t           EEn;
   Float_t         EEe[1000];
   Float_t         EEet[1000];
   Float_t         EEeta[1000];
   Float_t         EEphi[1000];
   Float_t         EEperp[1000];
   Bool_t          EEisjet[1000];
   Float_t         EEmax;
   Int_t           EBn;
   Float_t         EBe[1000];
   Float_t         EBet[1000];
   Float_t         EBeta[1000];
   Float_t         EBphi[1000];
   Float_t         EBperp[1000];
   Bool_t          EBisjet[1000];
   Float_t         EBmax;

   TTree *calTree = new TTree("Cal","Cal Tree");
   calTree->SetMaxTreeSize(MAXTREESIZE);
   calTree->Branch("HBHEn",&HBHEn,"HBHEn/I");
   calTree->Branch("HBHEe",HBHEe,"HBHEe[HBHEn]/F");
   calTree->Branch("HBHEeta",HBHEeta,"HBHEeta[HBHEn]/F");
   calTree->Branch("HBHEphi",HBHEphi,"HBHEphi[HBHEn]/F");
   calTree->Branch("HBHEperp",HBHEperp,"HBHEperp[HBHEn]/F");
   calTree->Branch("HBHEisjet",HBHEisjet,"HBHEisjet[HBHEn]/F");
   calTree->Branch("HBHEdepth",HBHEdepth,"HBHEdepth[HBHEn]/F"); 
   calTree->Branch("HBmax",&HBmax,"HBmax/F");
   calTree->Branch("HEmax",&HEmax,"HEmax/F");
   calTree->Branch("HFn",&HFn,"HFn/I");
   calTree->Branch("HFe",HFe,"HF[HFn]/F");
   calTree->Branch("HFeta",HFeta,"HFeta[HFn]/F");
   calTree->Branch("HFphi",HFphi,"HFphi[HFn]/F");
   calTree->Branch("HFperp",HFperp,"HFperp[HFn]/F");
   calTree->Branch("HFisjet",HFisjet,"HFisjet[HFn]/F");
   calTree->Branch("HFdepth",HFdepth,"HFdepth[HFn]/F");
   calTree->Branch("HFtotal",&HFtotal,"HFtotal/F");
   calTree->Branch("HFplus",&HFplus,"HFplus/F");
   calTree->Branch("HFminus",&HFminus,"HFminus/F");
   calTree->Branch("HFplusmax",&HFplusmax,"HFplusmax/F");
   calTree->Branch("HFminusmax",&HFminusmax,"HFminusmax/F");
   calTree->Branch("EEn",&EEn,"EEn/I");
   calTree->Branch("EEe",EEe,"EEe[EEn]/F");
   calTree->Branch("EEeta",EEeta,"EEeta[EEn]/F");
   calTree->Branch("EEphi",EEphi,"EEphi[EEn]/F");
   calTree->Branch("EEperp",EEperp,"EEperp[EEn]/F");
   calTree->Branch("EEisjet",EEisjet,"EEisjet[EEn]/F");
   calTree->Branch("EEmax",&EEmax,"EEmax/F");
   calTree->Branch("EBn",&EBn,"EBn/I");
   calTree->Branch("EBe",EBe,"EBe[EBn]/F");
   calTree->Branch("EBeta",EBeta,"EBeta[EBn]/F");
   calTree->Branch("EBphi",EBphi,"EBphi[EBn]/F");
   calTree->Branch("EBperp",EBperp,"EBperp[EBn]/F");
   calTree->Branch("EBisjet",EBisjet,"EBisjet[EBn]/F");
   calTree->Branch("EBmax",&EBmax,"EBmax/F");

   Int_t bin1;
   Int_t bin2;
   Float_t pT1;
   Float_t pT2;
   Float_t eff1;
   Float_t eff2;
   Float_t w1;
   Float_t w2;
   Float_t DjpT;

   TTree *CSTree = new TTree("CS","CS Tree");
   CSTree->SetMaxTreeSize(MAXTREESIZE);
   CSTree->Branch("bin1",&bin1,"bin1/I");
   CSTree->Branch("bin2",&bin2,"bin2/I");
   CSTree->Branch("pT1",&pT1,"pT1/F");
   CSTree->Branch("pT2",&pT2,"pT2/F");
   CSTree->Branch("eff1",&eff1,"eff1/F");
   CSTree->Branch("eff2",&eff2,"eff2/F");
   CSTree->Branch("w1",&w1,"w1/F");
   CSTree->Branch("w2",&w2,"w2/F");
   CSTree->Branch("DjpT",&DjpT,"DjpT/F");

   Int_t aljtnum;
   Float_t aljtpT[100];
   Float_t aljteta[100];
   Float_t aljtphi[100];
   Float_t aljtmass[100];

   TTree *jetTree = new TTree("aljets","jet Tree");
   jetTree->SetMaxTreeSize(MAXTREESIZE);
   jetTree->Branch("aljtnum",&aljtnum,"aljtnum/I");
   jetTree->Branch("aljtpT",aljtpT,"aljtpT[aljtnum]/F");
   jetTree->Branch("aljteta",aljteta,"aljteta[aljtnum]/F");
   jetTree->Branch("aljtphi",aljtphi,"aljtphi[aljtnum]/F");
   jetTree->Branch("aljtmass",aljtmass,"aljtmass[aljtnum]/F");

   Int_t mtwjtnum;
   Float_t mtwjtpT[100];
   Float_t mtwjteta[100];
   Float_t mtwjtphi[100];
   Float_t mtwjtmass[100];

   TTree *mtwjetTree = new TTree("mtwjets","3 or more jet Tree");
   mtwjetTree->SetMaxTreeSize(MAXTREESIZE);
   mtwjetTree->Branch("mtwjtnum",&mtwjtnum,"mtwjtnum/I");
   mtwjetTree->Branch("mtwjtpT",mtwjtpT,"mtwjtpT[mtwjtnum]/F");
   mtwjetTree->Branch("mtwjteta",mtwjteta,"mtwjteta[mtwjtnum]/F");
   mtwjetTree->Branch("mtwjtphi",mtwjtphi,"mtwjtphi[mtwjtnum]/F");
   mtwjetTree->Branch("mtwjtmass",mtwjtmass,"mtwjtmass[mtwjtnum]/F");

   Int_t twjtnum;
   Float_t twjtpT[100];
   Float_t twjteta[100];
   Float_t twjtphi[100];
   Float_t twjtmass[100];

   TTree *twjetTree = new TTree("twjets","2 jet Tree");
   twjetTree->SetMaxTreeSize(MAXTREESIZE);
   twjetTree->Branch("twjtnum",&twjtnum,"twjtnum/I");
   twjetTree->Branch("twjtpT",twjtpT,"twjtpT[twjtnum]/F");
   twjetTree->Branch("twjteta",twjteta,"twjteta[twjtnum]/F");
   twjetTree->Branch("twjtphi",twjtphi,"twjtphi[twjtnum]/F");
   twjetTree->Branch("twjtmass",twjtmass,"twjtmass[twjtnum]/F");

   Int_t thjtnum;
   Float_t thjtpT[100];
   Float_t thjteta[100];
   Float_t thjtphi[100];
   Float_t thjtmass[100];

   TTree *thjetTree = new TTree("thjets","3 jet Tree");
   thjetTree->SetMaxTreeSize(MAXTREESIZE);
   thjetTree->Branch("thjtnum",&thjtnum,"thjtnum/I");
   thjetTree->Branch("thjtpT",thjtpT,"thjtpT[thjtnum]/F");
   thjetTree->Branch("thjteta",thjteta,"thjteta[thjtnum]/F");
   thjetTree->Branch("thjtphi",thjtphi,"thjtphi[thjtnum]/F");
   thjetTree->Branch("thjtmass",thjtmass,"thjtmass[thjtnum]/F");

   Int_t mthjtnum;
   Float_t mthjtpT[100];
   Float_t mthjteta[100];
   Float_t mthjtphi[100];
   Float_t mthjtmass[100];

   TTree *mthjetTree = new TTree("mthjets","4 or more jet Tree");
   mthjetTree->SetMaxTreeSize(MAXTREESIZE);
   mthjetTree->Branch("mthjtnum",&mthjtnum,"mthjtnum/I");
   mthjetTree->Branch("mthjtpT",mthjtpT,"mthjtpT[mthjtnum]/F");
   mthjetTree->Branch("mthjteta",mthjteta,"mthjteta[mthjtnum]/F");
   mthjetTree->Branch("mthjtphi",mthjtphi,"mthjtphi[mthjtnum]/F");
   mthjetTree->Branch("mthjtmass",mthjtmass,"mthjtmass[mthjtnum]/F");

   Int_t nVertex;
   Float_t xVertex[100];
   Float_t yVertex[100];
   Float_t zVertex[100];

   TTree *vtxTree = new TTree("vtx","vertex information");
   vtxTree->SetMaxTreeSize(MAXTREESIZE);
   vtxTree->Branch("nVertex",&nVertex,"nVertex/I");
   vtxTree->Branch("xVertex",xVertex,"xVertex[nVertex]/F");
   vtxTree->Branch("yVertex",yVertex,"yVertex[nVertex]/F");
   vtxTree->Branch("zVertex",zVertex,"zVertex[nVertex]/F");

   Double_t d_phi;

   vector<float> orderpt;
   vector<float> ordere;
   vector<float> jetpt;
   vector<float> jeteta;
   vector<float> jetphi;
   vector<float> jetmass;
   vector<float>::iterator iter_jetpt;

   vector<float> HBe_cand;
   vector<float> HEe_cand;
   vector<float> HFeplus_cand;
   vector<float> HFeminus_cand;
   vector<float> EBe_cand;
   vector<float> EEe_cand;
   vector<float>::iterator iter_HBe_cand;
   vector<float>::iterator iter_HEe_cand;
   vector<float>::iterator iter_HFeplus_cand;
   vector<float>::iterator iter_HFeminus_cand;
   vector<float>::iterator iter_EBe_cand;
   vector<float>::iterator iter_EEe_cand;

   vector<float> FRP_cand;
   vector<float> BRP_cand;
   vector<float>::iterator iter_FRP_cand;
   vector<float>::iterator iter_BRP_cand;

   vector<float> jtcand[40000];
   vector<float> tkcand[40000];
   vector<float>::iterator iter_jtcand;
   vector<float>::iterator iter_tkcand;

   Float_t HB_max;
   Float_t HE_max;
   Float_t HFplus_max;
   Float_t HFminus_max;
   Float_t EB_max;
   Float_t EE_max;

   Float_t v1v2;
   
   TLorentzVector v1, v2;
   TVector2 n;

   Int_t numjt;
   Int_t BAD;
   // event loop start
   //   if(nevt == -1) nevt = pfTree->GetEntries();   

   Float_t rBRP;
   Float_t rFRP;

   Float_t bin = 0.2;
   Int_t numbin = (int)5/bin;

   Int_t numevt = 0;
   if(nevt == -1) nevt = hf->GetEntries(); 
   for (int iev=0; iev<nevt; iev++) {
     if(iev%100000==0)
     { 
       cout << ">>>>> EVENT " << iev << " / " << hf->GetEntries() <<  " ("<<(int)(100.*iev/hf->GetEntries()) << "%)" << endl; 
     }
     HltTree->GetEntry(iev);

   
     if ( (trig != "" ) && (!trigBit) ) {
       continue;
     }
     
     t->GetEntry(iev);
     //hiTree->GetEntry(iev);
     // trigger selection 

     v1v2 = -99;

     rBRP = 0.;
     rFRP = 0.;

     BAD = 0;
     ///// Call the values /////
     pT1 = 0.;
     pT2 = 0.;
     w1 = 0.;
     w2 = 0.;
     eff1 = 0.;
     eff2 = 0.;
     bin1 = 0.;
     bin2 = 0.;
     DjpT = 0.; 
     TLorentzVector totaljtvec;

     HFtotal = 0.;
     HFplus = 0.;
     HFminus = 0.;    
 
     HB_max = 0.;
     HE_max = 0.;
     HFplus_max = 0.;
     HFminus_max = 0.;
     EB_max = 0.;
     EE_max = 0.;

     HBe_cand.clear();
     HEe_cand.clear();
     HFeplus_cand.clear();
     HFeminus_cand.clear();
     EBe_cand.clear();
     EEe_cand.clear();

     HBHEn = hben;
     for(Int_t d = 0; d != hben; d++)
     {
       HBHEe[d] = hbee[d];
       HBHEet[d] = hbeet[d];
       HBHEeta[d] = hbeeta[d];
       HBHEphi[d] = hbephi[d];
       HBHEperp[d] = hbeperp[d];
       HBHEisjet[d] = hbeisjet[d];
       HBHEdepth[d] = hbedepth[d];
       if(hbeeta[d] < 1.5)
       {
         HBe_cand.push_back(hbee[d]);
       }
       if(hbeeta[d] > 1.5 && hbeeta[d] < 3.0)
       {
         HEe_cand.push_back(hbee[d]);
       }
     }
     HFn = hfn;
     for(Int_t e = 0; e != hfn; e++)
     {
       HFe[e] = hfe[e];
       HFet[e] = hfet[e];
       HFeta[e] = hfeta[e];
       HFphi[e] = hfphi[e];
       HFperp[e] = hfperp[e];
       HFisjet[e] = hfisjet[e];
       HFdepth[e] = hfdepth[e];
       if(HFeta[e] > 0)
       {
         HFeplus_cand.push_back(hfe[e]);
       }
       if(HFeta[e] < 0)
       {
         HFeminus_cand.push_back(hfe[e]);
       }
     }
     HFtotal = totalEtHF;
     HFminus = totalEtHFminus;
     HFplus = totalEtHFplus;
     //cout << "evt : " << iev << " HF : " << HFtotal << endl;
     EEn = een;
     for(Int_t f = 0; f != een; f++)
     {
       EEe[f] = eee[f];
       EEet[f] = eeet[f];
       EEeta[f] = eeeta[f];
       EEphi[f] = eephi[f];
       EEperp[f] = eeperp[f];
       EEisjet[f] = eeisjet[f];
       EEe_cand.push_back(eee[f]);
     }
     EBn = ebn;
     for(Int_t g = 0; g != ebn; g++)
     {
       EBe[g] = ebe[g];
       EBet[g] = ebet[g];
       EBeta[g] = ebeta[g];
       EBphi[g] = ebphi[g];
       EBperp[g] = ebperp[g];
       EBisjet[g] = ebisjet[g];
       EBe_cand.push_back(ebe[g]);
     }
     
     iter_HBe_cand = max_element(HBe_cand.begin(), HBe_cand.end());
     iter_HEe_cand = max_element(HEe_cand.begin(), HEe_cand.end());
     iter_HFeplus_cand = max_element(HFeplus_cand.begin(), HFeplus_cand.end());
     iter_HFeminus_cand = max_element(HFeminus_cand.begin(), HFeminus_cand.end());
     iter_EEe_cand = max_element(EEe_cand.begin(), EEe_cand.end());
     iter_EBe_cand = max_element(EBe_cand.begin(), EBe_cand.end()); 
     if(HBe_cand.size() > 0)
     {
       HB_max = *iter_HBe_cand;
     }
     else
     {
       HB_max = 0;
     }
     if(HEe_cand.size() > 0)
     {
       HE_max = *iter_HEe_cand;
     }
     else
     {
       HE_max = 0;
     }
     if(HFeplus_cand.size() > 0)
     {
       HFplus_max = *iter_HFeplus_cand;
     }
     else
     {
       HFplus_max = 0;
     }
     if(HFeminus_cand.size() > 0)
     {
       HFminus_max = *iter_HFeminus_cand;
     }
     else
     {
       HFminus_max = 0;
     }
     if(EEe_cand.size() > 0)
     {
       EE_max = *iter_EEe_cand;
     }
     else
     {
       EE_max = 0;
     }
     if(EBe_cand.size() > 0)
     {
       EB_max = *iter_EBe_cand;
     }
     else
     {
       EB_max = 0;
     }

     HBmax = HB_max;
     HEmax = HE_max;
     hemax = HE_max;
     HFplusmax = HFplus_max;
     HFminusmax = HFminus_max;
     EEmax = EE_max;
     EBmax = EB_max;
     //cout << "HF_max : " << HF_max << " HE_max : " << HE_max << " HB_max : " << HB_max << endl;

     if(HFplus_max >= 3.9 || HFminus_max >= 3.6)
     {
       continue;
     }
//     if(HFplus_max >= 3.9)
//     {
//       continue;
//     }
//     if(HFminus_max >= 3.6)
//     {
//       continue;
//     }    
 
     //cout << "HFplus_max : " << HFplus_max << "  HFminus_max : " << HFminus_max << endl;

     event.clear();
     event.run = run;
     event.lumi = lumi;
     event.event = evt; 
     event.vz = -99;
     event.hfsum = HFtotal; 
     event.hfplus = HFplus; 
     event.hfminus = HFminus;    

/*
     if(Run==262988 && LumiBlock==663 && Event==89487086)
     {
       cout << "262988 : 663 : 89487086  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89469683)
     {
       cout << "262988 : 663 : 89469683  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89482692)
     {
       cout << "262988 : 663 : 89482692  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89472230)
     {
       cout << "262988 : 663 : 89472230  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89553231)
     {
       cout << "262988 : 663 : 89553231  hfsum: " << HFtotal << endl;
     }
*/
     
     dj.clear();

     BRP_cand.clear();
     FRP_cand.clear();

     for(Int_t h = 0; h != 40000; h++)
     {
       jtcand[h].clear();
       tkcand[h].clear();
     }

     Track.clear();
     Track.nTrack = 0;
     ntrk = 0;
     floatntrk = 0;
     trkHFplus = HFplus;
     trkHFminus = HFminus;
     floatnTrkabsEtaover1p5 = 0;
     for(Int_t c = 0; c != nTrk; c++)
     {
       // Yongsun:  Here is where the quality cuts of trak will enter.
       //  pt cutk, pTerr/pt cut, etc etc in near fugure.  ask Daniel and Sasha which track quality cut you have you have to use.  
       if ( fabs(trkEta[c])> 2.4 ) 
	     continue;   // eta cut 
       //
       //
       //
       // From now you should not make any cut below :

       if(highPurity[c] == 1 && trkPtError[c]/trkPt[c] < 0.1 && fabs(trkDz1[c]/trkDzError1[c]) < 3 && fabs(trkDxy1[c]/trkDxyError1[c]) < 3 && trkNHit[c] >= 11 && trkChi2[c]/(float)trkNdof[c]/(float)trkNlayer[c] < 0.15 && trkPt[c] > 0.2)
       {
       trketa[ntrk] = trkEta[c];
       trkphi[ntrk] = trkPhi[c];
       trkpt[ntrk] = trkPt[c];
       //cout << iev << " " << trkEta[c] << " " << trkPhi[c] << " " << trkPt[c] << endl;
       BRP_cand.push_back(abs(trketa[ntrk] - NE));
       FRP_cand.push_back(abs(trketa[ntrk] - PE));
       //cout << BRP_cand.size() << " " << FRP_cand.size() << endl;
       ntrk = ntrk + 1;
       
       floatntrk = (float)ntrk;

       Track.nTrack = Track.nTrack + 1;

       if(fabs(trkEta[c]) > 1.5)
       {
         Track.nTrkabsEtaover1p5 = Track.nTrkabsEtaover1p5 + 1;
         floatnTrkabsEtaover1p5 = floatnTrkabsEtaover1p5 + 1;
       }
       if(fabs(trkEta[c]) < 1.5)
       {
         Track.nTrkabsEtaunder1p5 = Track.nTrkabsEtaunder1p5 + 1;
       }
       if(trkEta[c] >= -2.5 && trkEta[c] < -2.0)
       {
         Track.nTrketam2to2p5 = Track.nTrketam2to2p5 + 1;
       }
       else if(trkEta[c] >= -2 && trkEta[c] < -1.5)
       {
         Track.nTrketam1p5to2 = Track.nTrketam1p5to2 + 1; 
       }
       else if(trkEta[c] >= -1.5 && trkEta[c] < -1.0)
       {
         Track.nTrketam1to1p5 = Track.nTrketam1to1p5 + 1 ;
       }
       else if(trkEta[c] >= -1.0 && trkEta[c] < -0.5)   
       {
         Track.nTrketam0p5to1 = Track.nTrketam0p5to1 + 1;    
       }
       else if(trkEta[c] >= -0.5 && trkEta[c] < 0.0)
       {
         Track.nTrketam0to0p5 = Track.nTrketam0to0p5 + 1;
       }
       else if(trkEta[c] >= 0.0 && trkEta[c] < 0.5)
       {
         Track.nTrketa0to0p5 = Track.nTrketa0to0p5 + 1;
       }
       else if(trkEta[c] >= 0.5 && trkEta[c] < 1.0)
       {
         Track.nTrketa0p5to1 = Track.nTrketa0p5to1 + 1 ;
       }
       else if(trkEta[c] >= 1.0 && trkEta[c] < 1.5)
       {
         Track.nTrketa1to1p5 = Track.nTrketa1to1p5 + 1 ;
       }
       else if(trkEta[c] >= 1.5 && trkEta[c] < 2.0)
       {
         Track.nTrketa1p5to2 = Track.nTrketa1p5to2 + 1 ;
       }
       else if(trkEta[c] >= 2.0 && trkEta[c] < 2.5)
       {
         Track.nTrketa2to2p5 = Track.nTrketa2to2p5 + 1 ;
       }
       }
       else
       {
         continue;
       }
     }
     iter_BRP_cand = min_element(BRP_cand.begin(), BRP_cand.end());
     iter_FRP_cand = min_element(FRP_cand.begin(), FRP_cand.end());
     if(BRP_cand.size() > 0)
     {
       rBRP = *iter_BRP_cand;
     }
     else
     {
       rBRP = -1;
     }
     if(FRP_cand.size() > 0)
     {
       rFRP = *iter_FRP_cand;
     }
     else
     {
       rFRP = -1;
     }

     for(Int_t z = 0; z != numbin; z++)
     {
       if(rBRP >= 0.2*z && rBRP < 0.2*(z + 1))
       {
         BRP = (float)0.2*z;
       }
       if(rFRP >= 0.2*z && rFRP < 0.2*(z + 1))
       {
         FRP = (float)0.2*z;
       }
     }

     nVertex = nVtx;
     for(Int_t p = 0; p != nVtx; p++)
     {
       xVertex[p] = -99;
       yVertex[p] = -99;
       zVertex[p] = -99;

       xVertex[p] = xVtx[p];
       yVertex[p] = yVtx[p];
       zVertex[p] = zVtx[p];
     }  

     numjt = 0;

     orderpt.clear();
     ordere.clear();
     jetpt.clear();
     jeteta.clear();
     jetphi.clear();
     jetmass.clear();

     for(Int_t a = 0; a != nref; a++)      {
       TLorentzVector jt;
       jt.SetPtEtaPhiM( jtpt[a], jteta[a], jtphi[a], jtm[a] );
       // some jet cuts 
       //if(TMath::Abs(jteta[a]) < 2.0 && jtpt[a] > minjPt)
       //
       if(TMath::Abs(jteta[a]) > 1.8) // Should still have eta cut 
    	 continue ;
       
       // From here now, do not introduce any other cuts until the loop is over.
       //orderpt.push_back(jtpt[a]);
       ordere.push_back(jt.E());
       //cout << jt.E() << endl;
       jetpt.push_back(jtpt[a]);
       jeteta.push_back(jteta[a]);
       jetphi.push_back(jtphi[a]);
       jetmass.push_back(jtm[a]);
       totaljtvec += jt;
       numjt += 1;
       //cout << "a : " << a << "  totaljtvector mass : " << totaljtvec.M() << endl;
     }
     //cout << "num of jet : " << numjt << endl;
     if ( numjt == 0 )  { 
       jtpt2[0] = -1 ;
       jteta2[0] = -1 ;
       jtphi2[0] = -1 ;
       jtm2[0] = -1 ;
       jtpt2[1] = -1 ;
       jteta2[1] = -1 ;
       jtphi2[1] = -1 ;
       jtm2[1] = -1 ;

     }
     if ( numjt == 1 ) { 
       jtpt2[0] = jetpt[0];
       jteta2[0] = jeteta[0];
       jtphi2[0] = jetphi[0];
       jtm2[0]   = jetmass[0];

       jtpt2[1] = -1 ;
       jteta2[1] = -1 ;
       jtphi2[1] = -1 ;
       jtm2[1] = -1 ;

     }
     if(numjt > 1) {     // if numjt > 1 
       sort(ordere.begin(), ordere.end(), arrange);
       //cout << numjt << " " << jetpt.size() << " " << jetpt[0] << " " << jetpt[1] << " " << jetpt[2] << " " << jetpt[3] << endl;
       //if(numjt > 0 && jetpt.size() > 0 && numjt == jetpt.size())
       //{
       //  cout << "DAMN!" << endl;
       //}
       for(Int_t b = 0; b != nref; b++) {
       TLorentzVector jt2;
       jt2.SetPtEtaPhiM( jtpt[b], jteta[b], jtphi[b], jtm[b] );
       //cout << jt2.E() << " " << ordere[0] << " " << ordere[1] << endl;
    	 if((float)jt2.E() == (float)ordere[0])
  	   {
	       jtpt2[0] = jtpt[b];
	       jteta2[0] = jteta[b];
	       jtphi2[0] = jtphi[b];
  	     jtm2[0] = jtm[b];
         //cout << jtpt2[0] << endl;
	     }
    	 if((float)jt2.E() == (float)ordere[1])	
       {
  	     jtpt2[1] = jtpt[b];
	       jteta2[1] = jteta[b];
	       jtphi2[1] = jtphi[b];
         jtm2[1] = jtm[b];
         //cout << jtpt2[1] << endl;
	     }
       }
       //cout << ordere[0] << " " << ordere[1] << endl;
       //cout << jtpt2[0] << " " << jtpt2[1] << endl;
       //cout << jtm2[0] << " " << jtm2[1] << endl;
     }
     if(numjt > 0)
     {
       for(Int_t i = 0; i != numjt; i++)
       {
         aljtnum = numjt;
         aljtpT[i] = jetpt[i];
         aljteta[i] = jeteta[i];
         aljtphi[i] = jetphi[i];
         aljtmass[i] = jetmass[i];
         //cout << aljtpT[i] << endl;
         if(numjt == 2)
         {
           twjtnum = numjt;
           twjtpT[i] = jetpt[i];
           twjteta[i] = jeteta[i];
           twjtphi[i] = jetphi[i];
           twjtmass[i] = jetmass[i];
         }
         if(numjt == 3)
         {
           thjtnum = numjt;
           thjtpT[i] = jetpt[i];
           thjteta[i] = jeteta[i];
           thjtphi[i] = jetphi[i];
           thjtmass[i] = jetmass[i];
         }
         if(numjt > 3)
         {
           mthjtnum = numjt;
           mthjtpT[i] = jetpt[i];
           mthjteta[i] = jeteta[i];
           mthjtphi[i] = jetphi[i];
           mthjtmass[i] = jetmass[i];
         }
         if(numjt > 2)
         {
           mtwjtnum = numjt;
           mtwjtpT[i] = jetpt[i];
           mtwjteta[i] = jeteta[i];
           mtwjtphi[i] = jetphi[i];
           mtwjtmass[i] = jetmass[i];    
         }
       }
     }     
      
     d_phi = 0;
     d_phi = TMath::Abs(getDPHI(jtphi2[0], jtphi2[1]));
     if(d_phi > TMath::Pi())
     {
	     d_phi = 2 * 3.141592653589 - d_phi;
     }
     //cout << "d_phi : " << d_phi << endl;
     if(d_phi > 3.141592653589)
     {
	     cout << "Not good" << endl;
     }
      //cout << Run << " " << LumiBlock << " " << Event << endl;
      pt[0] = jtpt2[0];
      eta[0] = jteta2[0];
      phi[0] = jtphi2[0];
      mass[0] = jtm2[0];
      pt[1] = jtpt2[1];
      eta[1] = jteta2[1];
      phi[1] = jtphi2[1];
      mass[1] = jtm2[1];
      //cout << pt[0] << "  " << eta[0] << "  " << pt[1] << "  " << eta[1] << endl;     
      //cout << "pT1 : " << pt[0] << "  pT2 : " << pt[1] << "  d_phi : " << d_phi << "  total vectorsum inv mass : " << totaljtvec.M() << endl;
      TLorentzVector jt1vec, jt2vec, djvec, djrelvec, jt1vec2, jt2vec2;
      jt1vec.SetPtEtaPhiM( pt[0], eta[0], phi[0], mass[0] );
      jt2vec.SetPtEtaPhiM( pt[1], eta[1], phi[1], mass[1] );
      djvec = jt1vec + jt2vec ;
      djrelvec = 0.5*(jt1vec - jt2vec);
     
      if(numjt == 2 && pt[0] > 20 && pt[1] > 15 && totaljtvec.M() > 35 && d_phi > 2)
      {
      //cout << Run << " " << LumiBlock << " " << Event << endl;

      jt1vec2.SetPtEtaPhiE( pt[0], eta[0], phi[0], jt1vec.E() );
      jt2vec2.SetPtEtaPhiE( pt[1], eta[1], phi[1], jt2vec.E() );

      TVector2 p, v1;
      TVector2 q, v2;

      p.Set(jt1vec2[0], jt1vec2[1]);
      q.Set(jt2vec2[0], jt2vec2[1]);

      //cout << jt1vec[0] << " " << jt1vec[1] << endl;

      v1.Set(p.X() + q.X(), p.Y() + q.Y());
      v2.Set(0.5*(p.X() - q.X()), 0.5*(p.Y() - q.Y()));
     
      dj.v1_norm = TMath::Sqrt(v1.X() * v1.X() + v1.Y() * v1.Y());
      dj.v2_norm = TMath::Sqrt(v2.X() * v2.X() + v2.Y() * v2.Y());

      TVector2 v1unit, v2unit;
      v1unit.Set(v1.X() / dj.v1_norm,v1.Y() / dj.v1_norm);
      v2unit.Set(v2.X() / dj.v2_norm,v2.Y() / dj.v2_norm); 

      v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y();

      dj.c12 = v1v2;

      n.Set(v1unit.Y(),-v1unit.X()); 
      dj.n_norm = sqrt(n.X() * n.X() + n.Y() * n.Y());

      dj.s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y());

      dj.a12 = TMath::ATan2(dj.s12, dj.c12);
      if (dj.a12>=0) dj.a12 = dj.a12;
      if (dj.a12<0) dj.a12 = dj.a12 + 2*TMath::Pi();

      //cout << dj.a12 << endl;

      dj.c2phi  = cos(dj.a12) * cos(dj.a12) - sin(dj.a12) * sin(dj.a12);

//      cout << dj.c12 << " " << dj.s12 << " " << dj.a12 << " " << dj.c2phi << endl;

      dj.nJet = numjt;
      dj.mass = djvec.M();
      dj.pt   = djvec.Pt();
      dj.y   = djvec.Rapidity();
      dj.phi   = djvec.Phi();
      dj.eta  = djvec.Eta();
      dj.dphi = TMath::Abs(getDPHI( phi[0], phi[1] ));
      dj.dpt = TMath::Abs(pt[0] - pt[1]);
      dj.deta = TMath::Abs(eta[0] - eta[1]);
      dj.aj = TMath::Abs(pt[0] - pt[1])/(pt[0] + pt[1]);
      dj.pt1 = pt[0];
      dj.rpt1 = rawpt[0];
      dj.eta1 = eta[0];
      dj.phi1 = phi[0];
      dj.e1 = jt1vec.E();
      dj.pt2 = pt[1];
      dj.rpt2 = rawpt[1];
      dj.eta2 = eta[1];
      dj.phi2 = phi[1];
      dj.e2 = jt2vec.E();
    
      pT1 = jtpt2[0];
      pT2 = jtpt2[1]; 
      bin1 = heff->FindBin(pt[0]);
      eff1 = heff->GetBinContent(bin1);
      w1 = 1/eff1;
      bin2 = heff->FindBin(pt[1]);
      eff2 = heff->GetBinContent(bin2);
      w2 = 1/eff2;
      DjpT = djvec.Pt();
   
      eventTree->Fill();
      djTree->Fill();
      trkTree->Fill();   
      newtrkTree->Fill();
      trkheTree->Fill();
      calTree->Fill();
      CSTree->Fill();
      jetTree->Fill();
      mtwjetTree->Fill();
      twjetTree->Fill();
      thjetTree->Fill();
      mthjetTree->Fill();
      vtxTree->Fill();
      }
      else
      {
        continue;
      }
      orderpt.clear(); 
      jetpt.clear();
      jeteta.clear();
      jetphi.clear();
      jetmass.clear(); 
   } //end of event loop
   // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
   //   cout << numevt << endl;
   vtxTree->Write();
   mthjetTree->Write();
   thjetTree->Write();
   twjetTree->Write();
   mtwjetTree->Write();
   jetTree->Write();
   CSTree->Write();
   calTree->Write();
   trkheTree->Write();
   newtrkTree->Write();
   trkTree->Write();
   djTree->Write();
   eventTree->Write();
   
   newfile->Close();
   cout << "THE END" << endl;
} 



TString getDayAndTime() { 
  time_t currentTime;
  struct tm *localTime;
  
  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time
  
  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  //  int Sec    = localTime->tm_sec;
  return Form("%dy_%dm_%dd_%dh_%dm",Year,Month,Day,Hour,Min);
}
