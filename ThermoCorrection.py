# -*- coding: utf-8 -*-
"""
======================
ThermoCorrection.py
======================
version: 1.3
Author: Bene Hiebl
Date: 05.11.2021

An atmospheric and emissivity correction for georeferenced thermal infrared images from close range sensing
based on Caselles et al. (1996), Kodimalar et al. (2020), Wiecek (2011) and Minkina et al. (2016).

Input requirements:
•	Georeferenced thermo raster images (.tif) + timeinfo
•	DSM
•	FVC
•	Red Band (optional)
•	Camera location
•	Vertical gradients air temperature and humidity
•	Emissivity for soil/vegetation
•	Radiance air, downwelling radiance (optional)

define Input in Correction_Parameters.txt

"""


import os
import math as m
import rasterio
import numpy as np
import pandas as pd
from scipy.integrate import quad
import warnings




##############################################################################################
# FUNCTIONS
##############################################################################################

# reads information of paramterfile and returns dictionary
def ReadParameterfile(wd, file_name):
    "Reads Parameterfile from txt and give information to dictionary, returns dictionary"

    pfile = open(wd + "/" + file_name, "r")

    # create dictionary from parameters
    try:
        pdict = {}
        for line in pfile:

            # delete # comments
            if line.startswith("#") == False:
                # delete certain spaces that are not part of the file
                line = line.strip(" ")
                # delete \n
                line = line.strip("\n")
                # tabulator as seperator
                sep_line = line.split("\t")
                sep_line = list(sep_line)
                # values to dictionary
                pdict[sep_line[0]] = sep_line[2]

        pfile.close()
    except:
        print(
            "Something's wrong with the Parameter file! Probably you missed a /t (Tab)...")

    # convert numbers to float
    for i in pdict:
        try:
            pdict[i] = float(pdict[i])
        except:
            pdict[i] = str(pdict[i])

    return pdict

# create output directory


def OutDirectory(out_dir):
    "creates output directory"

    try:
        os.mkdir(out_dir)
    except OSError:
        print(f"Creation of the directory {out_dir} failed")
    else:
        print(f"Successfully created the directory {out_dir} ")
    return out_dir

# creates filelist


def FileList(in_dir):
    "creates string list of files in specific folder (format: foldername, no /)"

    image_dir = os.scandir(in_dir)
    imagenames = []
    for file in image_dir:
        name = file.name
        imagenames.append(name)
    return imagenames

# create fractional vegetation cover


def FvFromVI(VI_array, soil_limit, vegetation_limit):
    "Calculates fractional Vegetation cover from Vegetation index, returns Fv array"
    fv_cover = VI_array.copy()
    fv_cover[fv_cover < soil_limit] = 0
    fv_cover = np.where(VI_array > soil_limit, ((
        VI_array-soil_limit)/(vegetation_limit-soil_limit))**2, fv_cover)
    fv_cover = np.where(VI_array > vegetation_limit, 1, fv_cover)

    return fv_cover

# calculates land surface emissivity from Fractional Cover with thresholds for soil and vegetation


def LseFromFv(fv_cover, veg_lse, soil_lse, cavity):
    "Calculates land surface emissivity from fractional vegetation cover (Caselles 1996),returns lse as array"

    # de_value = fv_cover.copy()
    # de_value[de_value < 1] = cavity
    # de_value = np.where(fv_cover == 0, 0, de_value)
    lse = np.where((fv_cover < 1) & (fv_cover >= 0), veg_lse*fv_cover +
                   soil_lse*(1-fv_cover) + (4 * cavity * fv_cover*(1-fv_cover)), -99999)
    lse = np.where(fv_cover >= 1, veg_lse + cavity, lse)

    return lse

# calculate max humidity for given temperature


def MaxHumidity(T):
    "Calculates max possible humidity at temperature T from empirical formula, returns humidity"

    H = ((611.2*m.exp((22.46*T)/(272.62+T)))/(461.51*T))
    return H

# get raster coordinates from pixel in array of certain geotransformation


def RastCoordinates(i, j, transform, dsm_array):
    "computes raster coordinates from array coordinates, returns list=(x,y,z)"

    x = j*transform[0] + transform[2]
    y = i*transform[4] + transform[5]
    z = dsm_array[i, j]
    loc = (x, y, z)
    return loc

# profile_dsm["transform"]
# RastCoordinates(490,390,profile_dsm["transform"],dsm)


def Boltz():
    o = (5454781984210512994952000000 * (m.pi)**5) / \
        29438455734650141042413712126365436049
    return o


def SatVapourPressure(T_atm):
    "Calculates saturation vapour pressure in Pa from atmospheric temperature in K"

    if T_atm >= 273.15:
        ps = 288.68 * (1.098 + ((T_atm - 273.15)/100))**8.02
    else:
        ps = 4.689 * ((1.486 + ((T_atm - 273.15)/100))**12.3)

    return ps


def SatVapourPressureKG(T_atm):
    "Calculates saturation vapour pressure g/kg-3 from atmospheric temperature in K"

    if T_atm >= 273.15:
        ps = 6.112 * m.exp(17.269*((T_atm-273.17)/(T_atm-35.86)))
    else:
        ps = 6.112 * m.exp(21.874*((T_atm-273.17)/(T_atm-7.66)))

    return ps

def planck(lamda,T):
    return (1.19104*10**8) / (lamda**5 * (np.exp(14387.7 / (lamda * T)) - 1))

def planck_1(lamda,L):
    return 14387.7 / (lamda * np.log((1.19104*10**8/(lamda**5*L))+1))

# interpolated passman larmore tables for 7.5 to 13 microns for H2O and CO2 against newton height and distance d
index_passman = [7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14]
cols_passman = [0.0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
passman_h2o = [[1, 0.947, 0.874, 0.762, 0.582, 0.258, 0.066, 0.000],                  
               [1, 0.990, 0.975, 0.951, 0.904, 0.777, 0.603, 0.365],
               [None,None,None,None,None,None,None],
               [1, 0.997, 0.992, 0.984, 0.968, 0.921, 0.848, 0.719],
               [None,None,None,None,None,None,None],
               [1, 0.998, 0.994, 0.988, 0.975, 0.940, 0.883, 0.780],
               [None,None,None,None,None,None,None],
               [1, 0.998, 0.994, 0.988, 0.975, 0.940, 0.883, 0.779],
               [None,None,None,None,None,None,None],
               [1, 0.997, 0.993, 0.987, 0.974, 0.937, 0.878, 0.770],
               [None,None,None,None,None,None,None],
               [1, 0.997, 0.992, 0.984, 0.967, 0.921, 0.846, 0.718],
               [None,None,None,None,None,None,None],
               [None,None,None,None,None,None,None]]
passman_h2o = pd.DataFrame(passman_h2o, index=index_passman, columns=cols_passman)
passman_h2o = passman_h2o.interpolate(method="linear")

passman_co2 = [[1, 1, 1, 1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1, 1, 1, 1],
               [None,None,None,None,None,None,None],
               [1, 1, 1, 1, 1, 1, 1, 1],
               [None,None,None,None,None,None,None],
               [1, 1, 1, 0.999, 0.997, 0.994, 0.989, 0.978],
               [None,None,None,None,None,None,None],
               [1, 1, 0.999, 0.999, 0.997, 0.993, 0.986, 0.973],
               [None,None,None,None,None,None,None],
               [1, 1, 1, 0.999, 0.999, 0.997, 0.993, 0.986],
               [None,None,None,None,None,None,None],
               [1, 0.991, 0.977, 0.955, 0.912, 0.794, 0.630, 0.397],
               [None,None,None,None,None,None,None],
               [None,None,None,None,None,None,None]]
passman_co2 = pd.DataFrame(passman_co2, index=index_passman, columns=cols_passman)
passman_co2 = passman_co2.interpolate(method="linear")




############################################################################################################################################################################
# Start MAIN SCRIPT
############################################################################################################################################################################


# Read Parameter file and give to dictionary
print("Reading Parameter file.")
# Set working directory
wd = os.path.dirname(os.path.realpath("ThermoCorrection.py"))
os.chdir(wd)

pdict = ReadParameterfile(wd, "Correction_Parameters.txt")

out_dir = pdict["OUTPUT_LOC"]
dsm_name = pdict["DSM_INPUT"]
thermo_loc = pdict["THERMO_FOLDER"]
time_info = pdict["TIME_INFO"]
vi_name = pdict["VI_INPUT"]
fv_name = pdict["FV_INPUT"]
red_name = pdict["RED_INPUT"]

cam_location = (float(pdict["CAMLOCATION_X"]), float(
    pdict["CAMLOCATION_Y"]), float(pdict["CAMLOCATION_Z"]))

vertical_gradient = pdict["VERTICAL_GRAD"]
em_const = pdict["EM_CONST"]
if em_const == "NA":
    em_const = 0.005
else:
    em_const = float(em_const)

soil_lse = (float(pdict["LSE_VEGETATION"]), float(pdict["LSE_SOIL"]))
snow_lse = float(pdict["LSE_SNOW"])
lamda_min = float(pdict["LAMBDA_MIN"])
lamda_max = float(pdict["LAMBDA_MAX"])
    
cloud_fraction = pdict["CLOUD_COVER"]
if cloud_fraction == "NA":
    cloud_fraction = 0.0
else:
    cloud_fraction = float(cloud_fraction)
    


# create output directory
out_dir = OutDirectory(out_dir)

##########################################################################################################################
# READ INPUT FILES
##########################################################################################################################

# create filelist of input data and meteo data
files_thermo_all = FileList(thermo_loc)
files_thermo = []
for name in files_thermo_all:
    if name.endswith(".tif"):
        files_thermo.append(name)

# read fractional cover file and store as array
with rasterio.open(vi_name, "r") as img:
    vi = img.read(1)
    profile_vi = img.profile

# read fractional cover file and store as array
with rasterio.open(fv_name, "r") as img:
    fv_cover = img.read(1)
    profile_fv = img.profile

# read DSM file and store as array
with rasterio.open(dsm_name, "r") as img:
    dsm = img.read(1)
    profile_dsm = img.profile


# read meteo vertical profile
gradient = pd.read_csv(vertical_gradient, sep=",")
# read time info for thermo raster
timeinfo = pd.read_csv(time_info, sep=",")


# indexing dates in gradient data
gradient_timeindex = []
for i in range(0, len(gradient["time"])):
    time_val = 24*60*int(gradient["time"][i][:2]) + 60*int(
        gradient["time"][i][11:13]) + int(gradient["time"][i][14:16])
    gradient_timeindex.append(time_val)
gradient["timeindex"] = gradient_timeindex
timeinfo = timeinfo.set_index(timeinfo["Name"])

# indexing dates in tmeinfo thermo data
timeinfo_index = []
for i in range(0, len(timeinfo["Datetime"])):
    time_val = 24*60*int(timeinfo["Datetime"][i][:2]) + 60*int(
        timeinfo["Datetime"][i][11:13]) + int(timeinfo["Datetime"][i][14:16])
    timeinfo_index.append(time_val)
timeinfo["timeindex"] = timeinfo_index
timeinfo = timeinfo.set_index(timeinfo["Name"])



#################################################################################################################################
# EMISSIVITY
#################################################################################################################################
print("Start calculating EMISSIVITY...")

# Calculate emissivity from Fv cover (Caselles 1996, kodimalar 2020)
print("   Calculating LSE")

# LSE for vegetation and soil
lse_veg = LseFromFv(fv_cover, soil_lse[0], soil_lse[1], em_const)

# LSE for snow
if red_name != "NA":
    with rasterio.open(red_name, "r") as img:
        red = img.read(1)
    lse_veg = np.where(red >= 250, snow_lse, lse_veg)
    #lse_veg = np.where((vi < 0) & (lse_veg != -99999), 0.94, lse_veg)
else:
    print("   no red channel available.")

# LSE for rock
#lse_rock = np.where(landclass == rock, 0.93, lse_snow)

# LSE whole#
# lse = lse_rock.copy
lse = lse_veg.copy()


##############################################################################################
# Atmospheric correction
#############################################################################################

# Path length
print("Calculating PATH LENGTH")

# Creat x,y,z arrays for path length calculation
transform = profile_dsm["transform"]

y_ar = np.zeros(np.shape(dsm))
for i in range(np.shape(dsm)[0]):
    y_ar[i, :] = i*transform[4] + transform[5]

x_ar = np.zeros(np.shape(dsm))
for j in range(np.shape(dsm)[1]):
    x_ar[:, j] = j*transform[0] + transform[2]

# Calculate ground distance and height difference between XYZ and cam location
ground_dist = np.where(
    dsm != -99999, np.sqrt((cam_location[0]-x_ar)**2 + (cam_location[1]-y_ar)**2), -99999)
height_diff = np.where(
    dsm != -99999, np.sqrt((dsm - cam_location[2])**2), -99999)

# Calculate path distance between XYZ and cam location
path_ar = np.where(dsm != -99999, (np.sqrt(ground_dist **
                   2 + height_diff**2)) / 1000, -99999)


###################################################################################################
# Atmospheric Correction
print("Starting Atmospheric Correction")

# create dictionary for recalcultaion of T from Plancks law as errors occur during calculation
# T1 = np.arange(220,330,0.01)
# T2 = []
# for l in np.arange(220,330,0.01):
#     T = l
#     AT = quad(planck,7.5,14,args=(T))[0]        
#     L = AT/6.5
#     A = quad(planck_1,7.5,14,args=(L))[0]
#     T2.append(A/6.5)
# T_dict = dict(zip(np.round(T1,decimals=2),T1-T2))


# iterating over all thermo images
for nr in range(0, len(files_thermo)):
    print("   Preparing " + files_thermo[nr])

    # read time info and extract time of thermo img
    for a in range(0, len(timeinfo["timeindex"])):
        if timeinfo["Name"][a] == files_thermo[nr][:-11]:
            thermo_timeindex = timeinfo["timeindex"][a]

    # get correct coefficient and intersection for given timestep for T and RH
    # get index of forwarding value
    idx2 = pd.Index(gradient["timeindex"])
    gradientindex = idx2.get_loc(thermo_timeindex, method="ffill")
    # get
    T_gradient = gradient["TL_coefficient"][gradientindex] + (thermo_timeindex - gradient["timeindex"][gradientindex]) * (
        (gradient["TL_coefficient"][gradientindex+1]-gradient["TL_coefficient"][gradientindex]) / (gradient["timeindex"][gradientindex+1]-gradient["timeindex"][gradientindex]))
    T_intercept = gradient["TL_intercept"][gradientindex] + (thermo_timeindex - gradient["timeindex"][gradientindex]) * (
        (gradient["TL_intercept"][gradientindex+1]-gradient["TL_intercept"][gradientindex]) / (gradient["timeindex"][gradientindex+1]-gradient["timeindex"][gradientindex]))

    RH_gradient = gradient["RH_coefficient"][gradientindex] + (thermo_timeindex - gradient["timeindex"][gradientindex]) * (
        (gradient["RH_coefficient"][gradientindex+1]-gradient["RH_coefficient"][gradientindex]) / (gradient["timeindex"][gradientindex+1]-gradient["timeindex"][gradientindex]))
    RH_intercept = gradient["RH_intercept"][gradientindex] + (thermo_timeindex - gradient["timeindex"][gradientindex]) * (
        (gradient["RH_intercept"][gradientindex+1]-gradient["RH_intercept"][gradientindex]) / (gradient["timeindex"][gradientindex+1]-gradient["timeindex"][gradientindex]))

    T_cam = T_gradient * cam_location[2] + T_intercept
    RH_cam = RH_gradient * cam_location[2] + RH_intercept

    # correct if not Kelvin
    if T_cam < 200:
        T_cam = T_cam + 273.15

    # at cam values for RH
    # correct if not in 0. form
    if RH_cam > 2:
        RH_cam = RH_cam / 100
    H_air = RH_cam * MaxHumidity(T_cam)

    #######################################################################################################################
    # Mean path temperature and relative humidity
    print("      Calculating Path variables")
    
    #read thermo file
    with rasterio.open(thermo_loc+"/"+files_thermo[nr]) as img:
        thermo1 = img.read(1)
        profile_thermo = img.profile

    # check if in Kelvin or degree
    thermo2 = np.where(thermo1 == -99999, np.nan, thermo1)
    mean_th = np.nanmean(thermo2)
    del thermo2
    
    if mean_th < 100:
        thermo = np.nan_to_num(np.where(thermo1 != -99999, thermo1 + 273.15, -99999),nan=-99999)
    else:
        thermo = np.nan_to_num(thermo1,nan=-99999)
    thermo = np.where(thermo<=230, -99999, thermo)
        
        
    # Calculating mean Pathtemperature
    T_pix = np.where(dsm != -99999, (T_gradient * dsm + T_intercept) + 273.15, -99999)
    pathtemp = np.where(dsm != -99999, ((T_pix + T_cam) / 2) - 273.15, -99999)

    # Calculating mean Relative Humidity along path
    # from measured T at cam location
    #RH_pix = np.where(dsm != -99999, H_air / ((611.2*np.exp((22.46*T_pix)/(272.62+T_pix)))/(461.51*T_pix)), -99999)
    # from gradient calculated from stations
    RH_pix = np.where(dsm != -99999, RH_gradient * dsm + RH_intercept, -99999)
    pathRH = np.where(dsm != -99999, (RH_pix + RH_cam) / 2, -99999)

    #########################################################################################################################################
    # Calculating Atmospheric transmission as suggested by Minkina (2016) via Passman Larmore tables ()
    print("         Atmospheric Transmittance")
    # Calculating tao from passman larmore tables (experimentally) wavelength 7.5-13ym
    # as suggested by Minkina 2016
    
    # calculate cylinder height newton method    
    h_ar = np.where(thermo != -99999, (0.00016667 * pathtemp**3 + 0.01 * pathtemp**2 + 0.38333 * pathtemp + 5) * (pathRH/100) * path_ar, -99999)
    # calculate saturation vapour pressure
    p_sat = np.where((T_pix >= 273.15), 288.68 * (1.098 + ((T_pix - 273.15)/100))**8.02, 4.689 * ((1.486 + ((T_pix - 273.15)/100))**12.3))
    p_sat_kg = np.where((pathtemp + 273.15) >= 273.15, 6.112 * np.exp(17.269*(((pathtemp + 273.15)-273.17)/((pathtemp + 273.15)-35.86))), 6.112 * np.exp(21.874*(((pathtemp + 273.15)-273.17)/((pathtemp + 273.15)-7.66))))
    # emissivity air for spectral broadband 7.5 to 14
    e_h2o = np.where(thermo != -99999, 0.54 * (pathRH/100 * p_sat_kg)**(1/7), -99999)
    
    
    

    # iterating over all wavelengths
    
    tao_atm_list = []
    
    for lamda in index_passman: #np.arange(lamda_min,lamda_max+0.5,0.5):
        print("            "+str(lamda))
        
        p_h2o = np.ones(np.shape(thermo))
        p_co2 = np.ones(np.shape(thermo))
        
        for nr1 in range(0, len(cols_passman)-1):
            
            
            ## linearly interpolate values in passman larmore tables
            # for H2O            
            p_h2o = np.where((h_ar > cols_passman[nr1]) & (h_ar <= cols_passman[nr1+1]) & (thermo != -99999), ((h_ar - cols_passman[nr1]) / (cols_passman[nr1+1] - cols_passman[nr1])) * (
                passman_h2o[cols_passman[nr1+1]].loc[lamda] - passman_h2o[cols_passman[nr1]].loc[lamda]) + passman_h2o[cols_passman[nr1]].loc[lamda], p_h2o)
            # for CO2        
            p_co2 = np.where((path_ar > cols_passman[nr1]) & (path_ar <= cols_passman[nr1+1]) & (thermo != -99999), ((path_ar - cols_passman[nr1]) / (cols_passman[nr1+1] - cols_passman[nr1])) * (
                passman_co2[cols_passman[nr1+1]].loc[lamda] - passman_co2[cols_passman[nr1]].loc[lamda]) + passman_co2[cols_passman[nr1]].loc[lamda], p_co2)
        # calculate tao    
        tao_lamda = p_h2o * p_co2
        tao_atm_list.append(tao_lamda)
    
    tao = np.nanmean(tao_atm_list,axis=0)
    del tao_atm_list    
        
    
    ########################################################################################################
    ### ATMOSPHERIC CORRECTION
    ########################################################################################################
    print("       Start correcting...this may take a while.")
        
    T_surf = np.empty(np.shape(thermo)) * np.nan
    
    for i in range(0,np.shape(thermo)[0]):
        for j in range(0,np.shape(thermo)[1]):
            if (thermo[i,j] <= 230) | (lse[i,j] == -99999):
                T_surf[i,j] = -99999
                
            else:
                warnings.filterwarnings('ignore')
    
                ##### calculate downwelling radiance (Konzelmann et al. 1994), clear sky emissivity (Klok and Oerlemans 2002)
                #e_cs = 0.23 + 0.433 * ((RH_pix * SatVapourPressure(T_pix)) / T_pix)**(1/8)
                # total emissivity sky with cloud cover (e_cl = 0.976 given by Gruell et al. 1997)
                #e_sky = (0.23 + 0.433 * ((RH_pix * SatVapourPressure(T_pix)) / T_pix)**(1/8)) * (1 - cloud_fraction**2) + 0.976 * cloud_fraction**2    
                T = T_pix[i,j]
                #L_down_int = quad(planck,7.5,14,args=(T))    
                L_down = ((0.23 + 0.433 * (((RH_pix[i,j]/100) * p_sat[i,j]) / T)**(1/8)) * (1 - cloud_fraction**2) + 0.976 * cloud_fraction**2) * quad(planck,lamda_min,lamda_max,args=(T))[0]   #Boltz() * T_pix**4
        
        
                ### path radiance = LWIR atmosphärische eigenstrahlung corresponding to pathdistance/km Brutsaert 1975 and Siqueira Katul 2009
                T = pathtemp[i,j] + 273.15
                #L_path_int = quad(planck,7.5,14,args=(T))
                L_path = e_h2o[i,j] * (1-tao[i,j]) * quad(planck,lamda_min,lamda_max,args=(T))[0]   #* Boltz() * (pathtemp + 273.15)**4
                
        
                ### Corection based on jimenez 2009 and Richter 2002 and atmospheric transfer equation
                # Calculate longwave radiation at camera for wavelength lamda # at cam radiance by plancks law as mean for wavelengths 7.5 to 14
                T = thermo[i,j]
                L_cam = quad(planck,lamda_min,lamda_max,args=(T))[0]
                # Correcting Surface radiance for each pixel
                L = ((L_cam - L_path - tao[i,j] * (1 - lse[i,j]) * (L_down/m.pi)) / (tao[i,j] * lse[i,j])) / (lamda_max-lamda_min)
                
                # calculate T by reverse plancks law
                T_surf[i,j] = quad(planck_1,lamda_min,lamda_max,args=(L))[0] / (lamda_max-lamda_min)
                
                
    
   

    # Surface temperature in different units reverse of planck function
    if mean_th < 100:
        T_surf = np.where(T_surf != -99999, T_surf - 273.15, T_surf)
    
   
    #### save corrected surface temperature
    print("      Correcting "+files_thermo[nr])
    
    profile_thermo["dtype"] = T_surf.dtype
    with rasterio.open(f"{out_dir}/{files_thermo[nr][:-4]}_corrected.tif", "w", **profile_thermo) as dst:
        dst.write(T_surf, 1)


print("   done.")
    
#### Wiecek 2011 variation for tao calculation    
    # T_ref = 288.15
    # RH_ref = 50
    # d_ref = 1000
    # p_atm_ref = 0.8160

    # ps_alpha = SatVapourPressure(T_ref)
    # alpha = -((Boltz() * T_ref)/((RH_ref*100) * ps_alpha * d_ref))*m.log(p_atm_ref)

    # ### Calculate atmospheric transmission coefficient (tao)
    # print("      calculate Tao")

    # ### calculate saturation vapour pressure Minkina/Wiecek
    # p_sat = np.where(pathtemp >= 273.15, 6.112 * np.exp(17.269*((pathtemp-273.17)/(pathtemp-35.86))), 6.112 * np.exp(21.874*((pathtemp-273.17)/(pathtemp-7.66))))
    # p_sat = np.where(pathtemp == -99999, -99999, p_sat)

    # ### Calculate transmission coefficient
    # tao_atm = np.where(pathtemp != -99999, np.exp(-(alpha) * path_ar * ((p_sat*(pathRH*100))/(Boltz()*pathtemp))), -99999)




print("Save additional info")

# save lse
profile_dsm["dtype"] = lse.dtype
with rasterio.open(f"{out_dir}/LSE.tif", "w", **profile_dsm) as dst:
    dst.write(lse, 1)

# save tao
profile_thermo["dtype"] = tao_lamda.dtype
with rasterio.open(f"{out_dir}/Tao.tif", "w", **profile_thermo) as dst:
    dst.write(tao_lamda, 1)

