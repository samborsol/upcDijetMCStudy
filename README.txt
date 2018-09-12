Hello!

First, you will want to make a plots folders

>>mkdir plots
>>mkdir plots/png
>>mkdir plots/pdf

Now, you will want to make the data folders

>>mkdir data
>>mkdir data/skimmedFiles
>>mkdir data/hiforest
>>mkdir data/dataColumns

If you have the hiforest for the MC, put it in the folder "data/hiforest". Sasha has prepared a very good bash script for creating the hiforest in CMSSW. I have added this to the repository, under "runMC.sh", but he may have a more up to date version.

If you have the data-skimmed-file, put it in "data/skimmedFiles"

Run the "genPreSelect" forest skimmer. This will give you a skimmed-hiforest that only has gen jets that pass these cuts:

1. Only two gen jets are in the event.
2. The leading jet has genpt > 20 GeV.
3. The subleading jet has genpt > 15 GeV.
4. Both genjets have abs(geneta) < 1.8
5. The phi difference between the genjets is > 2. 
6. The invariant mass of the gen dijet is > 35 GeV.

Now, with this hiforest that only has what we want in it, perform Beomgon's dijet skim. 

Run the script "" to create a basic ntuple that is easy to make plots for. The script selects events that pass our coherent dijet selection, and contains MC truth, MC reco, and data.

Descriptions of plotting scripts:

1.
2.
3.
4.
5.
6.
