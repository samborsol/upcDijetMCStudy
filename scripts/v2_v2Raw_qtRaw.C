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

void v2_v2Raw_qtRaw(){

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

	TProfile *p1 = new TProfile("p1","p1",10,0,40);
	TProfile *p2 = new TProfile("p2","p2",10,0,40);
	TProfile *p3 = new TProfile("p3","p3",10,0,40);
	TProfile *p4 = new TProfile("p4","p4",10,0,40);

	TH1F *h1 = new TH1F("h1","h1",100,0,100);
	TH1F *h2 = new TH1F("h2","h2",100,0,100);
	TH1F *h3 = new TH1F("h3","h3",100,0,100);
	TH1F *h4 = new TH1F("h4","h3",100,0,100);

	TH1F *h1b = new TH1F("h1b","h1b",100,-3,3);
	TH1F *h2b = new TH1F("h2b","h2b",100,-3,3);
	TH1F *h3b = new TH1F("h3b","h3b",100,-3,3);
	TH1F *h4b = new TH1F("h4b","h3b",100,-3,3);

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
		h1b->Fill(v1.Eta());
		h1b->Fill(v2.Eta());
		cout<<"hello"<<endl;
		v1=v1+v2;
		p1->Fill(v1.Pt(),cos(2*phiAngle));
		cout<<v1.Pt()<<" "<<cos(2*phiAngle)<<endl;
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
		t1->GetEntry(i);
		v1.SetXYZ(px1,py1,pz1);
		v2.SetXYZ(px2,py2,pz2);


		p2->Fill(rawQT,cos(2*rawPhiAngle));
		p3->Fill(redQT,cos(2*redPhiAngle));

		h2b->Fill(v1.Eta());
		h2b->Fill(v2.Eta());
		v1=v1+v2;
		p4->Fill(v1.Pt(),cos(2*phiAngle));

	}

	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	p1->SetMaximum(1.0);
	p1->SetMinimum(0);
	p1->SetMarkerColor(kGreen);
	p1->SetLineColor(kGreen);
	p1->Draw();
	p2->GetXaxis()->SetTitle("QT (GeV/c)");
	p2->GetYaxis()->SetTitle("V2, <cos(2*phi)>");
	p2->SetMarkerColor(kRed);
	p2->SetLineColor(kRed);
	p2->Draw("same");
	p3->SetMarkerColor(kBlue);
	p3->SetLineColor(kBlue);
	p3->Draw("same");
	p4->SetLineColor(kBlue);
	p4->SetMarkerColor(kBlue);
	p3->SetLineStyle(2);
	p3->SetMarkerStyle(24);
	p4->Draw("same");
	TLegend *l1 = new TLegend(0.60,0.15,0.95,0.25);
	l1->AddEntry(p1,"Gen","p");
	l1->AddEntry(p2,"Reco, Uncorr.","p");
	l1->AddEntry(p3,"Reco, Red. JEC","p");
	l1->AddEntry(p4,"Reco, JEC","p");
	l1->Draw();
	c1->SaveAs("../plots/png/v2_v2Raw_qtRaw.png");

}
