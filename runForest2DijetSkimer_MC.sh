
#void forest2diJetSkim(
#          TString fname = "/home/kbg777/CMSwork/4_UPCTriggers_Beomgon.root",
#	  TString outputFname = "upcDiJetSkim",
#	  TString trig = "",  # if trig=="" then the trigger cut is not functioning
#	  TString jetCollection = "ak5PFJetAnalyzer", or  "akPu5PFJetAnalyzer",
#	  float minjPt = 30,
#	  int nevt=-1  ) 

# pp di-jet data 
#Input file is at MIT server 
inputFile="data/hiforest/genPreSelect.root" 
jetCollection="ak4PFJetAnalyzer"
#trig="HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1"
trig=""
for  minjPt in 0
do
    root -l -b -q 'skimmer/forest2diJetSkim_pp6.C+("'$inputFile'","closureMC","'$trig'","'$jetCollection'",'$minjPt',-1)'
done
# pp PYTHIA MC :
# input file is at XXX

# PbPb UPC triggered event 

# PbPb high pt jet triggered event 

# STARLIGHT DPMJET event


