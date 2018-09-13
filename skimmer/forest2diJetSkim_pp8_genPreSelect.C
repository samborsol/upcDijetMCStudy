#include <ctime>
#include <TMath.h>
#include "../interface/cutsAndBin_bgk.h"
#include <TLorentzVector.h>
#include "../interface/commonUtility.h"
#include "Riostream.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <TEventList.h>

using std::stringstream;
using std::vector;
using std::abs;

static const long MAXTREESIZE = 10000000000;

bool arrange(float i, float j) {return (i>j);}
TString getDayAndTime();

void forest2diJetSkim_pp8_genPreSelect(
		TString fname = "/u/user/bekim/758p3_UPC/src/Muon_test/4_UPCTriggers_pp_reco_170427.root",
		TString outputFname = "upcDiJetSkim180827", 
		TString trig = "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1",  
		TString jetCollection = "ak4PFJetAnalyzer", 
		float minjPt = 0,
		int nevt=-1
		) {  

	using namespace std;
	TFile *f1 = new TFile(fname.Data());

	TTree *t = (TTree*)f1->Get(Form("%s/t",jetCollection.Data()));
	TEventList *elist = new TEventList();

        TFile *newfile = new TFile("data/skimmedFiles/genPreSelect_dphi.root","recreate");


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

	t->SetBranchAddress("ngen", &nref, &b_nref);
	t->SetBranchAddress("genpt", jtpt, &b_jtpt);
	t->SetBranchAddress("geneta", jteta, &b_jteta);
	t->SetBranchAddress("geny", jty, &b_jty);
	t->SetBranchAddress("genphi", jtphi, &b_jtphi);
	t->SetBranchAddress("genm", jtm, &b_jtm);

	////////////////////////////////////////////////////////////////////////
	//////////////////  dijet tree 
	////////////////////////////////////////////////////////////////////////
	UPCEvent event;
	TTree *eventTree = new TTree("evt","Event Tree");
	eventTree->SetMaxTreeSize(MAXTREESIZE);
	eventTree->Branch("event",&event.run,eventBranchString.Data());
	UPCdiJet dj;
	TTree *djTree = new TTree("dijet","Jet Tree");
	djTree->SetMaxTreeSize(MAXTREESIZE);
	djTree->Branch("dj",&dj.nJet,djBranchString);

	Int_t bin1;
	Int_t bin2;
	Float_t pT1;
	Float_t pT2;
	Float_t eff1;
	Float_t eff2;
	Float_t w1;
	Float_t w2;
	Float_t DjpT;

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

	Float_t rBRP;
	Float_t rFRP;

	Float_t bin = 0.2;
	Int_t numbin = (int)5/bin;


	Int_t numevt = 0;
	if(nevt == -1) nevt = t->GetEntries(); 
	for (Int_t iev=0; iev<nevt; iev++) {

		t->GetEntry(iev);

		if(nref==2){

			BAD = 0;
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
				ordere.push_back(jt.E());
				jetpt.push_back(jtpt[a]);
				jeteta.push_back(jteta[a]);
				jetphi.push_back(jtphi[a]);
				jetmass.push_back(jtm[a]);
				totaljtvec += jt;
				numjt += 1;
			}
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
				for(Int_t b = 0; b != nref; b++) {
					TLorentzVector jt2;
					jt2.SetPtEtaPhiM( jtpt[b], jteta[b], jtphi[b], jtm[b] );
					if((float)jt2.E() == (float)ordere[0])
					{
						jtpt2[0] = jtpt[b];
						jteta2[0] = jteta[b];
						jtphi2[0] = jtphi[b];
						jtm2[0] = jtm[b];
					}
					if((float)jt2.E() == (float)ordere[1])	
					{
						jtpt2[1] = jtpt[b];
						jteta2[1] = jteta[b];
						jtphi2[1] = jtphi[b];
						jtm2[1] = jtm[b];
					}
				}
			}

			if(numjt == 2){

				d_phi = 0;
				d_phi = TMath::Abs(getDPHI(jtphi2[0], jtphi2[1]));
				if(d_phi > TMath::Pi())
				{
					d_phi = 2 * 3.141592653589 - d_phi;
				}
				if(d_phi > 3.141592653589)
				{
					cout << "Not good" << endl;
				}

				pt[0] = jtpt2[0];
				eta[0] = jteta2[0];
				phi[0] = jtphi2[0];
				mass[0] = jtm2[0];
				pt[1] = jtpt2[1];
				eta[1] = jteta2[1];
				phi[1] = jtphi2[1];
				mass[1] = jtm2[1];
				TLorentzVector jt1vec, jt2vec, djvec, djrelvec, jt1vec2, jt2vec2;
				jt1vec.SetPtEtaPhiM( pt[0], eta[0], phi[0], mass[0] );
				jt2vec.SetPtEtaPhiM( pt[1], eta[1], phi[1], mass[1] );
				djvec = jt1vec + jt2vec ;
				djrelvec = 0.5*(jt1vec - jt2vec);

				if(numjt == 2 && pt[0] > 20 && pt[1] > 15 && totaljtvec.M() > 35 && abs(eta[0])<1.8 && abs(eta[1])<1.8 )
				{


					jt1vec2.SetPtEtaPhiE( pt[0], eta[0], phi[0], jt1vec.E() );
					jt2vec2.SetPtEtaPhiE( pt[1], eta[1], phi[1], jt2vec.E() );

					TVector2 p, v1;
					TVector2 q, v2;

					p.Set(jt1vec2[0], jt1vec2[1]);
					q.Set(jt2vec2[0], jt2vec2[1]);


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

					dj.c2phi  = cos(dj.a12) * cos(dj.a12) - sin(dj.a12) * sin(dj.a12);

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

					djTree->Fill();
					elist->Enter(iev);
				}
			}
		}
	}
	TTree *HltTree = (TTree*)f1->Get("hltanalysis/HltTree");
	TTree *hiTree  = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");
	TTree *trackTree  = (TTree*)f1->Get("ppTrack/trackTree");
	//TTree *hbhe = (TTree*)f1->Get("rechitanalyzer/hbhe");
	//TTree *hf = (TTree*)f1->Get("rechitanalyzer/hf");
	//TTree *ee = (TTree*)f1->Get("rechitanalyzer/ee");
	//TTree *eb = (TTree*)f1->Get("rechitanalyzer/eb");


	t->SetEventList(elist);
	TTree *s = t->CopyTree("");
	TTree *s2 = HltTree->CopyTree("");
	TTree *s3 = hiTree->CopyTree("");
	TTree *s4 = trackTree->CopyTree("");

	djTree->Write();
	s->Write();
	s2->Write();
	s3->Write();
	s4->Write();

	newfile->Close();
} 

