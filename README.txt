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

If you have the hiforest for the MC, put it in the folder "data/hiforest"

If you have the data-skimmed-file, put it in "data/skimmedFiles"

Run the two forest skimmers. These will give you a hiforest that only has gen jets that pass these cuts:

1.) Only two gen jets
2.) The leading jet is genpt>20
3.) The subleading jet is genpt>15
4.) Both genjets have abs(geneta)<1.8
5.) The phi difference between the genjets is greater than 2. 

