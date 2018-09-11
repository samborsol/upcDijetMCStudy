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


void ang_resolution(){

	TString treeFolder = "/home/samboren/Workspace/upcAnalysis/mcClosure/data/skimmedFiles/";
	TFile *fMCTrig = new TFile(treeFolder+"TrkCutsppreco_closureMC_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0.root");

	//Now doing the reco MC shit
	TTree *dijetTree_MC = (TTree*)fMCTrig->Get("dijet");
	UPCdiJet djObj;
	dijetTree_MC->SetBranchAddress("dj",&djObj);
	TTree *jetTree_MC = (TTree*)fMCTrig->Get("t");	
	dijetTree_MC->AddFriend(jetTree_MC);

	Float_t jetpt[2];
	Float_t genpt[2];

	TH1F *h1 = new TH1F("h1","h1",100,.8,1.2);

	jetTree_MC->SetBranchAddress("jtphi",&jetpt);
	jetTree_MC->SetBranchAddress("genphi",&genpt);

	Int_t nEvents = dijetTree_MC->GetEntries();
	for(Long64_t i=0; i<nEvents; i++){
		dijetTree_MC->GetEntry(i);


		if( djObj.nJet!=2 || abs(djObj.eta1)>1.8 || abs(djObj.eta2)>1.8 || djObj.dphi<2 || djObj.pt1<20 || djObj.pt2<15 || djObj.mass<35) continue;

		if( abs(jetpt[0]-genpt[0])<abs(jetpt[1]-genpt[0]) ){

			h1->Fill(jetpt[0]/genpt[0]);				
			h1->Fill(jetpt[1]/genpt[1]);

		}
		else{
			h1->Fill( jetpt[0]/genpt[1] );
			h1->Fill( jetpt[1]/genpt[0] );

		}
	}

        setTDRStyle();

        TCanvas *c1 = new TCanvas("c1","c1",800,800);

        gStyle->SetStatX(0.9);
        gStyle->SetStatY(0.6);

        h1->SetLineWidth(2);
        h1->GetXaxis()->SetTitle("Jet Resolution in Azimuthal Angle, #phi");
        h1->GetYaxis()->SetTitle("Events per bin");
	h1->Draw();
        TF1* ft1 = new TF1("ft1","gaus",0,3);
        h1->Fit(ft1);
        TLegend *l1 = new TLegend(0.60,0.75,0.95,0.90);
        l1->AddEntry(h1,"Reco-#phi/Gen-#phi","l");
        l1->Draw();
        c1->SaveAs("../plots/png/ang_resolution.png");

}
