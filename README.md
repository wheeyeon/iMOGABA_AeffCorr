# iMOGABA_AeffCorr
A python script to correct for the re-determined KVN antenna efficiencies without re-running the iMOGABA pipeline.
This code uses the Python libraries astropy and numpy to calculate the observed source elevation per station per visibility.
This, along with the ANTAB DPFU+GC and newly-determined DPFU+GC are used to re-calibrate each visibility.
The observation date, frequency, and source coordinates are read from the UVFITS header. These values must be accurate. 
For station coordinates, we use the values presented in the EAVN status report, assuming that ANT1, 2, 3 are KTN, KUS, and KYS respectively. The code will need some minor edits if this is not the case. Future updates will read the station information directly from the AN table.
