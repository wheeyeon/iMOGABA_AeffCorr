# iMOGABA_AeffCorr
A python script to correct for the re-determined KVN antenna efficiencies without re-running the iMOGABA pipeline.
This code uses the Python libraries astropy and numpy to calculate the observed source elevation per station per visibility.
This, along with the ANTAB DPFU+GC and newly-determined DPFU+GC are used to re-calibrate each visibility.
The observation date, frequency, and source coordinates are read from the UVFITS header. These values must be accurate. 
For station coordinates, we use the values presented in the EAVN status report, assuming that ANT1, 2, 3 are KTN, KUS, and KYS respectively. The code will need some minor edits if this is not the case. Future updates will read the station information directly from the AN table.
The KVNSD_Tbin_imogaba.csv contains the start and end date of each iMOGABA epoch, as well as an additional scaling factor for KYS K-band observations to correct for additional flux loss most likely caused by phase instability for certain epochs. It is assumed that both the .py file and this .csv file are in the same directory.
Ultimately, the code is run by passing the path of the input and output UVFITS file to the function "scale_data()". One may edit the code manually for each UVFITS file, or alternatively create a list of input and output file paths (ex : 
total_file_list = os.popen("ls %s*.KVN27"%data_dir).read().split("\n")[:-1]) and iterate over this list.
The typical time required to correct an individual iMOGABA UVFITS file is approximately 1 ~ 2 minutes. However, it was noted that this code is very slow when run on Windows, with times exceeding a few hours per source. We recomment that this code be run on a Linux environment. 
