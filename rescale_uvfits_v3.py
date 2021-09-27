# Rework of separte codes for K/Q-band observations
# For use with K/Q-band LCP observations after 2015-09 (a.k.a. since KVN-2015B)
# Will Update to use station data from AN table directly in the future
# Assumes ANT1 = KTN, ANT2 = KUS, ANT3 = KYS

from astropy.io import fits, ascii
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5       # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude     # Angles
import numpy as np
from tqdm import tqdm

###############################################
# Edit here
file_name_mas = '0235+164'
infile_uvfits = '{0}.KVN15'.format(file_name_mas,)
otfile_uvfits = '{0}_AeffCorr3.KVN15'.format(file_name_mas,)
###############################################

###############################################
#     Util functions related to telescope     #
###############################################
def get_angle(val_deg, val_min, val_sec):
	output = Angle((val_deg, val_min, val_sec),unit=u.deg)
	return output

def tele_loc_lib(telename):
    teleloclib={'GBT'      : EarthLocation(lat = get_angle(38,25,59), lon = -79.839722*u.deg, height = 807.43*u.m),
                'nobeyama' : EarthLocation.from_geocentric(x=-3871025.4987, y=3428107.3984, z=3724038.7361, unit=u.m),
                'takahagi' : EarthLocation.from_geocentric(x=-3961881.8250, y=3243372.4800, z=3790687.4490, unit=u.m),
                'tianma'   : EarthLocation.from_geocentric(x=-2826708.6380, y=4679237.0440, z=3274667.5330, unit=u.m),
                'nanshan'  : EarthLocation.from_geocentric(x=  228310.1700, y=4631922.7550, z=4367064.0740, unit=u.m),
                'miznao20' : EarthLocation.from_geocentric(x=-3857244.9731, y=3108782.9179, z=4003899.1695, unit=u.m),
                'iriki'    : EarthLocation.from_geocentric(x=-3521719.8813, y=4132174.6817, z=3336994.1132, unit=u.m),
                'ogasa20'  : EarthLocation.from_geocentric(x=-4491068.3826, y=3481545.2394, z=2887399.8018, unit=u.m),
                'ishigaki' : EarthLocation.from_geocentric(x=-3263995.2619, y=4808056.3902, z=2619948.6347, unit=u.m),
                'kvnyonsei': EarthLocation.from_geocentric(x=-3042281.0183, y=4045902.6730, z=3867374.3296, unit=u.m),
                'kvnulsan' : EarthLocation.from_geocentric(x=-3287268.6453, y=4023450.1367, z=3687379.9886, unit=u.m),
                'kvntamna' : EarthLocation.from_geocentric(x=-3171731.6665, y=4292678.5393, z=3481038.7880, unit=u.m),
                'mopra'    : EarthLocation.from_geocentric(x=-4682769.05850, y=2802619.04217, z=-3291759.33837, unit=u.m),
                'yebes'    : EarthLocation.from_geocentric(x=4848761.7579, y=-261484.0570, z=4123085.1343, unit=u.m),}
    return teleloclib[telename]

def kys_anomoly(obs_date):
    # KYS scaling factor due to K-band phase instability
    # obs_date in JD
    # Notes : IM50 SF may be different per source
    #       : IM52 observation failed
    #       : IM61 not in image repository. Setting to 1 atm.
    kvn_epoch_dat   = ascii.read('KVNSD_Tbin_imogaba.csv')
    kvn_epoch_mask  = kvn_epoch_dat[Time(kvn_epoch_dat['LowerTbin']).to_value('jd') <= obs_date]
    kvn_epoch_final = kvn_epoch_mask[Time(kvn_epoch_mask['UpperTbin']).to_value('jd') >= obs_date]
    if len(kvn_epoch_final) == 0:
        return ('Unknown',1)
    else:
        return (kvn_epoch_final['Epoch'][0], kvn_epoch_final['KYS_SF'][0])

def get_gc_info(obs_date, obs_freq):
    # For returning old and new gc values to use
    # Will eventually update to query external table
    if obs_freq < 30:
        ##########################################################################
        antab_gc3 = (0.06716256, 1.04235321,-0.00080929, 0.00000376)# KYS-ANTAB
        if Time(obs_date) < Time('2016-03-15'):
            antab_gc2 = (0.07686854, 1.01825705, 0.00019177,-0.00000438)# KUS-ANTAB (2013 ver)
        else:
            antab_gc2 = (0.07625360, 1.01825705, 0.00019177,-0.00000438)# KUS-ANTAB (2016 ver)
        antab_gc1 = (0.07237465, 1.01635444, 0.00042356,-0.00000673)# KTN-ANTAB
        ##########################################################################
        wycjupgc3 = (0.07508119, 0.99688671, 0.00014536,-0.00000171)# KYS
        if Time(obs_date) < Time('2020-08-31'):
            wycjupgc1 = (0.08155963, 0.98931915, 0.00038022,-0.00000338)# KTN
        else:
            wycjupgc1 = (0.07694384, 0.98931915, 0.00038022,-0.00000338)# KTN
        if Time(obs_date) < Time('2018-08-31'):
            wycjupgc2 = (0.07196967, 1.00000000,-0.00052698, 0.00000442)# KUS
        else:
            wycjupgc2 = (0.08249817, 1.00000000,-0.00052698, 0.00000442)# KUS
        ##########################################################################
    else:
        ##########################################################################
        antab_gc3 = (0.07263276, 1.01009575, 0.00128409,-0.00001551)# KYS-ANTAB
        if Time(obs_date) < Time('2016-03-15'):
            antab_gc2 = (0.07240779, 1.06107499, 0.00050293,-0.00001313)# KUS-ANTAB (2013 ver)
        else:
            antab_gc2 = (0.07299647, 1.06107499, 0.00050293,-0.00001313)# KUS-ANTAB (2016 ver)
        antab_gc1 = (0.07374193, 1.00319584, 0.00234228,-0.00002642)# KTN-ANTAB
        ##########################################################################
        if Time(obs_date) < Time('2016-02-09'): #<--date may be off
            wycjupgc3 = (0.07689246, 0.94946968, 0.00208465,-0.00002150)# KYS
        elif Time(obs_date) < Time('2016-08-31'):
            wycjupgc3 = (0.07185862, 0.94946968, 0.00208465,-0.00002150)# KYS
        elif Time(obs_date) < Time('2017-05-19'): #<--date may be off by couple of days
            wycjupgc3 = (0.06279863, 0.94946968, 0.00208465,-0.00002150)# KYS
        else:
            wycjupgc3 = (0.08013851, 0.99688671, 0.00014536,-0.00002150)# KYS
        if Time(obs_date) < Time('2016-11-20'): #<--date may be off
            wycjupgc1 = (0.07375735, 0.97748671, 0.00100397,-0.00001119)# KTN
        else:
            wycjupgc1 = (0.08054063, 0.97748671, 0.00100397,-0.00001119)# KTN, Aeff scatter "very" large. May need to check individual epochs.
        if Time(obs_date) < Time('2016-01-15'): #<--date may be off by couple of days
            wycjupgc2 = (0.07569678, 1.00000000,-0.00078043,-0.00000001)# KUS
        else:
            wycjupgc2 = (0.08374899, 1.00000000,-0.00078043,-0.00000001)# KUS
        ##########################################################################
    return (antab_gc1, antab_gc2, antab_gc3, wycjupgc1, wycjupgc2, wycjupgc3)


def get_tar_altaz(tar_coord, datentime, format = 'jd', obs_tele = 'mopra'):
    # Code for getting Altitude/Azimuth of target at telescope
    # tar_coord : Target coordinates as astropy ICRS object
    # datentime : Date and time
    # obs_tele  : Name of telescope to use
    aptime    = Time(datentime,format=format)
    apframe   = AltAz(obstime = aptime, location = tele_loc_lib(obs_tele))
    tar_coord = SkyCoord(tar_coord)
    tar_altaz = tar_coord.transform_to(apframe)
    return tar_altaz


###############################################
#        Main UVFITS Editting Functions       #
###############################################

def gain_curve_poly(ele,gc_param_dat):
    # ele : elevation in degrees.
    dpfu,a0,a1,a2 = gc_param_dat
    return  dpfu*(a0*ele**2 + a1*ele + a2)

def scale_data(file_name, save_name):
    # For scaling data in bl1 
    # UVFITS BL convension
    # BL_12 = 2**(7+ant1)+ant2
    input_data = fits.open(file_name)
    data_table = input_data['PRIMARY'].data
    data_headr = input_data['PRIMARY'].header
    obs_date   = data_headr.get('DATE-OBS')
    obs_source = data_headr.get('OBJECT')
    obs_cfreq  = data_headr.get('CRVAL4')*1.e-9 # Hz to GHz
    obsepochif = kys_anomoly(data_table[0]['DATE'])
    print('Found source {0}, observed on {1} ({2}) at {3:.1f} GHz'.format(obs_source, obs_date, obsepochif[0], obs_cfreq))
    if Time(obs_date) < Time('2015-08-31'):
        print('Observation date not supported by this code.')
        print('Terminating...')
        return
    if obs_cfreq > 50: # No W/D band correction yet
        print('Observation frequency not supported by this code.')
        print('Terminating...')
        return
    #
    if data_headr.get('CTYPE6') == 'RA' and data_headr.get('CTYPE7') == 'DEC':
        # Coord in Deg.
        obs_RA     = data_headr.get('CRVAL6')
        obs_DEC    = data_headr.get('CRVAL7')
        obs_skycrd = ICRS(ra=obs_RA*u.deg, dec=obs_DEC*u.deg)
    else:
        print('Found incompatible header! Check manually!')
        print('Terminating...')
        return
    #
    # To get appropriate GC info 
    antab_gc1, antab_gc2, antab_gc3, wycjupgc1, wycjupgc2, wycjupgc3 = get_gc_info(obs_date, obs_cfreq)
    #
    print('Correcting Data')
    for i in tqdm(range(len(data_table))):
        bl_i    = data_table[i]['BASELINE']
        ant_2_i = int(bl_i % 256)
        ant_1_i = int(bl_i // 256)
        bl_st_l = [ant_1_i, ant_2_i]
        #print(bl_st_l)
        bl_i_da = data_table[i]['DATA']
        scale_fact = 1
        if 1 in bl_st_l:
            ant1_source_elev = (get_tar_altaz(obs_skycrd, data_table[i]['DATE'], obs_tele = 'kvntamna').alt).value
            sf_ant1 = gain_curve_poly(ant1_source_elev, antab_gc1) / gain_curve_poly(ant1_source_elev, wycjupgc1)
            scale_fact = np.sqrt(sf_ant1) * scale_fact
        if 2 in bl_st_l:
            ant2_source_elev = (get_tar_altaz(obs_skycrd, data_table[i]['DATE'], obs_tele = 'kvnulsan').alt).value
            sf_ant2 = gain_curve_poly(ant2_source_elev, antab_gc2) / gain_curve_poly(ant2_source_elev, wycjupgc2)
            scale_fact = np.sqrt(sf_ant2) * scale_fact
        if 3 in bl_st_l:
            ant3_source_elev = (get_tar_altaz(obs_skycrd, data_table[i]['DATE'], obs_tele = 'kvnyonsei').alt).value
            sf_ant3 = gain_curve_poly(ant3_source_elev, antab_gc3) / gain_curve_poly(ant3_source_elev, wycjupgc3)
            scale_fact = np.sqrt(sf_ant3) * scale_fact * obsepochif[1]
        #else:
        #    scale_fact = 1
        data_shape = np.shape(bl_i_da)
        for dat_i1 in range(data_shape[0]):
            for dat_i2 in range(data_shape[1]):
                for dat_i3 in range(data_shape[2]):
                    for dat_i4 in range(data_shape[3]):
                        for dat_i5 in range(data_shape[4]):
                                data_table[i]['DATA'][dat_i1][dat_i2][dat_i3][dat_i4][dat_i5][0] = scale_fact*bl_i_da[dat_i1][dat_i2][dat_i3][dat_i4][dat_i5][0]
                                data_table[i]['DATA'][dat_i1][dat_i2][dat_i3][dat_i4][dat_i5][1] = scale_fact*bl_i_da[dat_i1][dat_i2][dat_i3][dat_i4][dat_i5][1]
                                data_table[i]['DATA'][dat_i1][dat_i2][dat_i3][dat_i4][dat_i5][2] = (bl_i_da[dat_i1][dat_i2][dat_i3][dat_i4][dat_i5][2])/(scale_fact**2)
    #
    print('Saving to new file...')
    input_data.writeto(save_name)
    print('Done!')

scale_data(infile_uvfits, otfile_uvfits)
