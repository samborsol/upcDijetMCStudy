#Source the input hiforest
inputFile="data/hiforest/rg_eGlu_qQbar_fullMe_ipro13_v17.root" 
#Assign the jet collection from the hiforest
jetCollection="ak4PFJetAnalyzer"
#pick out a HLT path, or leave it blank for no trigger
#trig="HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1"
trig=""
#if you want to set some minimum jet pt values, for loop over the list of minima. I don't want to do that, so I have it set to 0. 
for  minjPt in 0
do
    root -l -b -q 'skimmer/forest2diJetSkim_pp6_genPreSelect.C+("'$inputFile'","closureMC","'$trig'","'$jetCollection'",'$minjPt',-1)'
done

