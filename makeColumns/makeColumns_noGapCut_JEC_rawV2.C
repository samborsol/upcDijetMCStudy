#include "../interface/cutsAndBin_bgk.h"
#include "../interface/commonUtility.h"
#include "TProfile.h"
#include "TProfile2D.h"

//No gap cuts applied

void makeColumns_noGapCut_JEC_rawV2(){
	//This script is for plotting the V2 vs PT. The event mixing effect is subtracted. 
	//Both BRP and FRP are summed up. GEN-level RAPGAP is also histogramed.

	TString treeFolder = "/home/samboren/Workspace/upcAnalysis/mcClosure/data/skimmedFiles/";
	TFile *fMCTrig = new TFile(treeFolder+"TrkCutsppreco_closureMC_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0.root");
	TFile *fDATA = new TFile(treeFolder+"TrkCutsppreco_upcDiJetSkim180827_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0_2018y_8m_27d_23h_34m.root");
	treeFolder = "/home/samboren/Workspace/upcAnalysis/mcClosure/data/skimmedFiles/";
	TFile *fGEN = new TFile(treeFolder+"genPreSelect.root");

	TTree *dijetTree_DATA   = (TTree*)fDATA->Get("dijet");
	TTree *trkTree_DATA     = (TTree*)fDATA->Get("fullTrkTree");
	TTree *calTree_DATA = (TTree*)fDATA->Get("Cal");
	TTree *evtTree_DATA = (TTree*)fDATA->Get("evt");
	TTree *vtxTree_DATA = (TTree*)fDATA->Get("vtx");
	dijetTree_DATA->AddFriend(trkTree_DATA);
	dijetTree_DATA->AddFriend(calTree_DATA);
	dijetTree_DATA->AddFriend(evtTree_DATA);
	dijetTree_DATA->AddFriend(vtxTree_DATA);

	UPCdiJet djObj;
	Float_t HFplusmax;
	Float_t HFminusmax;
	TLorentzVector jet1a;
	TLorentzVector jet2a;
	TLorentzVector jet1b;
	TLorentzVector jet2b;
	TLorentzVector jet1c;
	TLorentzVector jet2c;
	Float_t vz1;

	TLorentzVector jetplus;
	TLorentzVector jetminus;

	TLorentzVector jetplusGen;
	TLorentzVector jetminusGen;

	Float_t zee = 0.5;
	Float_t anglePerp;
	Float_t anglePerpGen;
	Float_t rawAnglePerp;
	Float_t redAnglePerp;

	Float_t QT;
	Float_t rawQT;
	Float_t redQT;

	Float_t d_phi;

	dijetTree_DATA->SetBranchAddress("dj",&djObj);
	calTree_DATA->SetBranchAddress("HFplusmax",&HFplusmax);
	calTree_DATA->SetBranchAddress("HFminusmax",&HFminusmax);

	Float_t dataVtx;
	vtxTree_DATA->SetBranchAddress("zVertex",&dataVtx);

	UPCEvent event;
	evtTree_DATA->SetBranchAddress("event",&event);	

	Float_t nEvents;

	TLorentzVector v01, v02, vv1, vv2, vqt, vpt, m, w, g, h, p, q;
	TVector2 n;

	//Define variables
	Double_t pi = TMath::Pi();
	float v1_norm, v2_norm, n_norm, m_norm;
	float c12, s12, v1v2, a12;
	Float_t pt1, pt2, eta1, eta2, phi1, phi2, e1, e2;
	Float_t px1, py1, pz1, ee1, px2, py2, pz2, ee2;
	float ptvqt, ptvpt;
	float sign=0.0;

	Float_t vtx1,vtx2,frp1,brp1,frp2,brp2,phiAngle;
	Float_t redPhiAngle;
	gROOT->Reset();

	TFile *fmd = new TFile("../data/dataColumns/columnsNoGap_JEC.root","RECREATE");

	//The pure data, nominal cuts
	nEvents=dijetTree_DATA->GetEntries();
	TNtuple *dataNtuple = new TNtuple("dataNtuple","data from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:vtx1:phiAngle:rpt1:rpt2:rawPhiAngle:QT:rawQT");
	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_DATA->GetEntry(i);
		if( djObj.mass < 35 || djObj.pt1<20 || djObj.pt2<15 || djObj.nJet!=2 || abs(djObj.eta1)>1.8 || abs(djObj.eta2)>1.8 || djObj.dphi<2) continue;

		vz1 = event.vz;

		if(djObj.e1>djObj.e2){
			jet1a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);

			jet1b.SetPtEtaPhiE(djObj.rpt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2b.SetPtEtaPhiE(djObj.rpt2,djObj.eta2,djObj.phi2,djObj.e2);

		}else{
			jet2a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);

			jet2b.SetPtEtaPhiE(djObj.rpt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1b.SetPtEtaPhiE(djObj.rpt2,djObj.eta2,djObj.phi2,djObj.e2);
		}

		jetplus=jet1a+jet2a;
		jetminus=(1-zee)*jet1a - zee*jet2a;
		TVector2 p, v1;
		TVector2 q, v2;
		p.Set(jet1a[0],jet1a[1]);
		q.Set(jet2a[0],jet2a[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		TVector2 v1unit, v2unit;
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);			
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		anglePerp=a12;
		QT=jetplus.Pt();

		jetplus=jet1b+jet2b;
		jetminus=(1-zee)*jet1b - zee*jet2b;
		p.Set(jet1b[0],jet1b[1]);
		q.Set(jet2b[0],jet2b[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		rawAnglePerp=a12;
		rawQT=jetplus.Pt();

		dataNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2a.E(),jet2a.Px(),jet2a.Py(),jet2a.Pz(),dataVtx,anglePerp,jet1b.Pt(),jet2b.Pt(),rawAnglePerp,QT,rawQT);
	}


	TTree *dijetTree_GEN = (TTree*)fGEN->Get("dijet");
	UPCdiJet djObjGen;
	dijetTree_GEN->SetBranchAddress("dj",&djObjGen);
	nEvents = dijetTree_GEN->GetEntries();
	//nEvents = 100;

	TNtuple *MCNtuple = new TNtuple("MCNtuple","gen mc from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:vtx1:phiAngle:QT");
	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_GEN->GetEntry(i);
		if(  djObjGen.mass < 35 || djObjGen.pt1<20 || djObjGen.pt2<15 || djObjGen.nJet!=2 || abs(djObjGen.eta1)>1.8 || abs(djObjGen.eta2)>1.8 || djObjGen.dphi<2) continue;

		vz1 = event.vz;
		if(djObjGen.e1>djObjGen.e2){
			jet1a.SetPtEtaPhiE(djObjGen.pt1,djObjGen.eta1,djObjGen.phi1,djObjGen.e1);
			jet2a.SetPtEtaPhiE(djObjGen.pt2,djObjGen.eta2,djObjGen.phi2,djObjGen.e2);
		}else{
			jet2a.SetPtEtaPhiE(djObjGen.pt1,djObjGen.eta1,djObjGen.phi1,djObjGen.e1);
			jet1a.SetPtEtaPhiE(djObjGen.pt2,djObjGen.eta2,djObjGen.phi2,djObjGen.e2);
		}

		jetplusGen=jet1a+jet2a;
		jetminusGen=(1-zee)*jet1a - zee*jet2a;
		TVector2 p, v1;
		TVector2 q, v2;
		p.Set(jet1a[0],jet1a[1]);
		q.Set(jet2a[0],jet2a[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		TVector2 v1unit, v2unit;
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);			
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		anglePerpGen=a12;
		MCNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2a.E(),jet2a.Px(),jet2a.Py(),jet2a.Pz(),vz1,anglePerpGen);
	}

	//Now doing the reco MC shit
	TTree *dijetTree_MC = (TTree*)fMCTrig->Get("dijet");
	dijetTree_MC->SetBranchAddress("dj",&djObj);
        TTree *vtxTree_MC = (TTree*)fMCTrig->Get("vtx");
	dijetTree_MC->AddFriend(vtxTree_MC);
	vtxTree_MC->SetBranchAddress("zVertex",&vz1);	

	nEvents = dijetTree_MC->GetEntries();

	TNtuple *MCRecoNtuple = new TNtuple("MCRecoNtuple","reco mc from tree","e1:px1:py1:pz1:e2:px2:py2:pz2:phiAngle:rpt1:rpt2:rawPhiAngle:rawQT:vtx1");
	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_MC->GetEntry(i);
		if( djObj.mass < 35 || djObj.pt1<20 || djObj.pt2<15|| djObj.nJet!=2 || abs(djObj.eta1)>1.8 || abs(djObj.eta2)>1.8 || djObj.dphi<2) continue;

		if(djObj.e1>djObj.e2){
			jet1a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);

			jet1b.SetPtEtaPhiE(djObj.rpt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2b.SetPtEtaPhiE(djObj.rpt2,djObj.eta2,djObj.phi2,djObj.e2);

			jet1c.SetPtEtaPhiE(1.089*djObj.rpt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet2c.SetPtEtaPhiE(1.078*djObj.rpt2,djObj.eta2,djObj.phi2,djObj.e2);

		}else{
			jet2a.SetPtEtaPhiE(djObj.pt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1a.SetPtEtaPhiE(djObj.pt2,djObj.eta2,djObj.phi2,djObj.e2);

			jet2b.SetPtEtaPhiE(djObj.rpt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1b.SetPtEtaPhiE(djObj.rpt2,djObj.eta2,djObj.phi2,djObj.e2);

			jet2c.SetPtEtaPhiE(1.078*djObj.rpt1,djObj.eta1,djObj.phi1,djObj.e1);
			jet1c.SetPtEtaPhiE(1.089*djObj.rpt2,djObj.eta2,djObj.phi2,djObj.e2);

		}


		jetplus=jet1a+jet2a;
		jetminus=(1-zee)*jet1a - zee*jet2a;
		TVector2 p, v1;
		TVector2 q, v2;
		p.Set(jet1a[0],jet1a[1]);
		q.Set(jet2a[0],jet2a[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		TVector2 v1unit, v2unit;
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);			
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		anglePerpGen=a12;
		QT=jetplus.Pt();

		jetplus=jet1b+jet2b;
		jetminus=(1-zee)*jet1b - zee*jet2b;
		p.Set(jet1b[0],jet1b[1]);
		q.Set(jet2b[0],jet2b[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		rawAnglePerp=a12;
		rawQT=jetplus.Pt();

		jetplus=jet1c+jet2c;
		jetminus=(1-zee)*jet1c - zee*jet2c;
		p.Set(jet1c[0],jet1c[1]);
		q.Set(jet2c[0],jet2c[1]);
		//Computing Qt and Pt as 2-Vectors
		v1.Set(p.X() + q.X(),p.Y()+q.Y());
		v2.Set(0.5*(p.X()-q.X()),0.5*(p.Y()-q.Y()));
		//computing the norm of Qt and Pt
		v1_norm = sqrt (  v1.X() * v1.X() + v1.Y() * v1.Y()    );
		v2_norm = sqrt ( v2.X() * v2.X() + v2.Y() * v2.Y()   );
		//Making unit vectors of Qt and Pt, resulting in Qt-hat and Pt-hat, unit vectors
		v1unit.Set(v1.X() / v1_norm,v1.Y() / v1_norm);
		v2unit.Set(v2.X() / v2_norm,v2.Y() / v2_norm);
		//Computing the dot product of Qt-hat and Pt-hat
		v1v2 = v1unit.X() * v2unit.X() + v1unit.Y() * v2unit.Y()   ;
		//The dot product is the cosine of the angle
		c12 = v1v2  ;
		//Define a perpendicular angle to Qt-hat, in order to compute the sine of the angle
		n.Set(v1unit.Y(),-v1unit.X()) ;
		n_norm = sqrt ( n.X()*n.X() + n.Y()*n.Y()  );
		//Sine of the angle
		s12 = (n.X()*v2unit.X() + n.Y()*v2unit.Y()  ) ;
		//Computing the angle by using arctan2 function, and considering the sign, so will be from 0,2*pi
		a12 = atan2(s12, c12);
		if (a12>=0) a12 = a12;
		if (a12<0) a12 = a12 + 2*pi;
		//Computing the cos(2phi) using trigonometry expression
		c12  = cos(a12) * cos(a12) - sin(a12) * sin(a12);
		redQT = jetplus.Pt();
		redAnglePerp=a12;

		MCRecoNtuple->Fill(jet1a.E(),jet1a.Px(),jet1a.Py(),jet1a.Pz(),jet2a.E(),jet2a.Px(),jet2a.Py(),jet2a.Pz(),anglePerpGen,jet1b.Pt(),jet2b.Pt(),rawAnglePerp,rawQT,vz1);

	}

	fmd->Write();
}
