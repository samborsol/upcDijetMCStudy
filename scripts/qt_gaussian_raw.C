#include "../interface/cutsAndBin_bgk.h"
#include "../interface/commonUtility.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "../interface/tdrstyle.C"
#include "../interface/CMS_lumi.C"


void SetTPaveTxt( TPaveText * txtemplate, int txtsize ) {
	txtemplate->SetFillColor(0);
	txtemplate->SetBorderSize(0);
	txtemplate->SetTextFont(43);
	txtemplate->SetTextAlign(12);
	txtemplate->SetTextSize(txtsize);
}

void SetLegend( TLegend * legtemplate, int legsize ) {
	legtemplate->SetFillColor(0);
	legtemplate->SetBorderSize(0);
	legtemplate->SetTextFont(43);
	legtemplate->SetTextSize(legsize);
}

void qt_gaussian_raw(){

	/*
	   Float_t min = 35;
	   Float_t max = 40;
	   Float_t butt = 3.0;
	   Float_t minG = 35;
	   Float_t maxG = 40;
	   Float_t buttG = 3.0;
	 */

	Float_t min = 20;
	Float_t max = 25;
	Float_t butt = 0;

	Float_t minG = min;
	Float_t maxG = max;
	Float_t buttG = butt;

	setTDRStyle();

	TFile *f1 = new TFile("../data/dataColumns/columnsNoGap_JEC.root");
	TTreeReader reader("T", f1);

	TH1F *h1 = new TH1F("h1","h1",40,0,40);
	TH1F *h2 = new TH1F("h2","h2",40,0,40);
	TH1F *h3 = new TH1F("h3","h3",40,0,40);
	TH1F *h4 = new TH1F("h4","h3",40,0,40);

	Float_t eng1, eng2, px1, px2, py1, py2, pz1, pz2, phiAngle, rawPhiAngle, QT, rawQT, redQT, redPhiAngle;
	Float_t rpt1, rpt2;
	Int_t nentries;
	TVector3 v1, v2, v1r, v2r;

	TTree *t1 = (TTree*)f1->Get("MCNtuple");
	t1->SetBranchAddress("phiAngle",&phiAngle);
	t1->SetBranchAddress("QT",&QT);
	t1->SetBranchAddress("e1",&eng1);
	t1->SetBranchAddress("e2",&eng2);
	t1->SetBranchAddress("px1",&px1);
	t1->SetBranchAddress("px2",&px2);
	t1->SetBranchAddress("py1",&py1);
	t1->SetBranchAddress("py2",&py2);
	t1->SetBranchAddress("pz1",&pz1);
	t1->SetBranchAddress("pz2",&pz2);
	nentries=t1->GetEntries();
	for(Int_t i=0; i<nentries; i++){
		t1->GetEntry(i);
		v1.SetXYZ(px1,py1,pz1);
		v2.SetXYZ(px2,py2,pz2);
		v1 = v1 + v2;
		h1->Fill(v1.Pt());
	}

	t1 = (TTree*)f1->Get("MCRecoNtuple");
	t1->SetBranchAddress("phiAngle",&phiAngle);
	t1->SetBranchAddress("rawPhiAngle",&rawPhiAngle);
	t1->SetBranchAddress("redPhiAngle",&redPhiAngle);
	t1->SetBranchAddress("redQT",&redQT);	
	t1->SetBranchAddress("rawQT",&rawQT);
	t1->SetBranchAddress("rpt1",&rpt1);
	t1->SetBranchAddress("rpt2",&rpt2);
	t1->SetBranchAddress("e1",&eng1);
	t1->SetBranchAddress("e2",&eng2);
	t1->SetBranchAddress("px1",&px1);
	t1->SetBranchAddress("px2",&px2);
	t1->SetBranchAddress("py1",&py1);
	t1->SetBranchAddress("py2",&py2);
	t1->SetBranchAddress("pz1",&pz1);
	t1->SetBranchAddress("pz2",&pz2);
	nentries=t1->GetEntries();
	for(Int_t i=0; i<nentries; i++){
		t1->GetEntry(i);

		h2->Fill(rawQT);


		v1.SetXYZ(px1,py1,pz1);
		v2.SetXYZ(px2,py2,pz2);
		v1 = v1+v2;
		h3->Fill(v1.Pt());
	}

	TCanvas *c1 = new TCanvas("c1","c1",800,800);

	gStyle->SetStatX(0.8);
	gStyle->SetStatY(0.5);

	h1->SetLineWidth(2);
	h1->SetMarkerColor(kGreen);
	h1->SetLineColor(kGreen);
//	h1->Draw();
	h2->GetXaxis()->SetTitle("Uncorr. QT (GeV/c)");
	h2->GetYaxis()->SetTitle("Events per bin");

	
	h2->SetMarkerColor(kBlue);
	h2->SetLineColor(kBlue);
	h2->Draw("");

	h3->SetMarkerColor(kBlue);
	h3->SetLineColor(kBlue);
	h3->SetLineWidth(2);
//	h3->Draw("same");
	h3->SetLineStyle(2);
	h3->SetMarkerStyle(24);
	TF1* ft1 = new TF1("ft1","gaus",0,35);
	h2->Fit(ft1);

	TLegend *l1 = new TLegend(0.60,0.75,0.95,0.90);
//	l1->AddEntry(h1,"Gen","p");
//	l1->AddEntry(gaus,"Gaussian Fit to Reco, JEC","p");
	l1->AddEntry(h2,"Reco, Uncorr.","l");
        l1->AddEntry(ft1,"Gaussian Fit to Reco, Uncorr.","l");
	l1->Draw();
	c1->SaveAs("../plots/png/qt_gaussian_raw.png");

}
