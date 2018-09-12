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


void dphi_stat(){

	TString treeFolder = "/home/samboren/Workspace/upcAnalysis/mcClosure/data/skimmedFiles/";
	TFile *fMCTrig = new TFile(treeFolder+"TrkCutsppreco_dphiMC_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0.root");

	//Now doing the reco MC shit
	TTree *dijetTree_MC = (TTree*)fMCTrig->Get("dijet");
	UPCdiJet djObj;
	dijetTree_MC->SetBranchAddress("dj",&djObj);
	TTree *jetTree_MC = (TTree*)fMCTrig->Get("t");	
	dijetTree_MC->AddFriend(jetTree_MC);

	Int_t ngen;
	dijetTree_MC->SetBranchAddress("ngen",&ngen);	
	
	TH1F *h1 = new TH1F("h1","h1",100,0,3.14159);
        TH1F *h2 = new TH1F("h2","h2",100,0,3.14159);

	Int_t nEvents = dijetTree_MC->GetEntries();
	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_MC->GetEntry(i);


		if( ngen!=2 || djObj.nJet!=2 || abs(djObj.eta1)>1.8 || abs(djObj.eta2)>1.8 || djObj.pt1<20 || djObj.pt2<15 || djObj.mass<35) continue	
		h1->Fill(djObj.dphi);
	}


        setTDRStyle();

        TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->SetLogy();
        gStyle->SetStatX(0.9);
        gStyle->SetStatY(0.6);
        h1->SetLineWidth(2);
        h1->GetXaxis()->SetTitle("#Delta #phi (azimuthal)");
        h1->GetYaxis()->SetTitle("Events per bin");
	h1->Draw();
        c1->SaveAs("../plots/png/mult_fake_stat.png");

}
