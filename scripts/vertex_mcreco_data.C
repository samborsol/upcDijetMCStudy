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

void vertex_mcreco_data(){

	setTDRStyle();

	TFile *f1 = new TFile("../data/dataColumns/columnsNoGap_JEC.root");
	TTreeReader reader("T", f1);

	TH1F *h1 = new TH1F("h1","h1",40,-25,25);
	TH1F *h2 = new TH1F("h2","h2",40,-25,25);
	TH1F *h3 = new TH1F("h3","h3",40,-25,25);
	TH1F *h4 = new TH1F("h4","h3",40,-25,25);

	Float_t vtx1, eng1, eng2, px1, px2, py1, py2, pz1, pz2, phiAngle, rawPhiAngle, QT, rawQT, redQT, redPhiAngle;
	Float_t rpt1, rpt2;
	Int_t nentries;
	TVector3 v1, v2, v1r, v2r;


	TTree *t1 = (TTree*)f1->Get("MCRecoNtuple");
	t1->SetBranchAddress("phiAngle",&phiAngle);
	t1->SetBranchAddress("rawPhiAngle",&rawPhiAngle);
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
	t1->SetBranchAddress("vtx1",&vtx1);
	nentries=t1->GetEntries();
	for(Int_t i=0; i<nentries; i++){
		t1->GetEntry(i);
		h1->Fill(vtx1);
	}

        t1 = (TTree*)f1->Get("dataNtuple");
        t1->SetBranchAddress("phiAngle",&phiAngle);
        t1->SetBranchAddress("rawPhiAngle",&rawPhiAngle);
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
	t1->SetBranchAddress("vtx1",&vtx1);
        nentries=t1->GetEntries();
        for(Int_t i=0; i<nentries; i++){
                t1->GetEntry(i);
		h2->Fill(vtx1);
        }


	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->SetLogy();
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);
	h1->SetMarkerColor(kRed);
	h1->SetLineColor(kRed);
	h2->Draw();
	h2->GetYaxis()->SetTitle("Events per bin");
	h2->GetXaxis()->SetTitle("Vz (cm) Dist.");
	h2->SetMarkerColor(kBlack);
	h2->SetLineColor(kBlack);
	h1->Draw("same");
	TLegend *l1 = new TLegend(0.20,0.9,0.4,0.8);
	l1->AddEntry(h1,"Reco MC","p");
	l1->AddEntry(h2,"DATA","p");
	l1->Draw();
	c1->SaveAs("../plots/png/vertex_mcreco_data.png");

}
