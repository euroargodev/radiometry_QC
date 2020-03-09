##############################################################################
# This program was designed in order to perform the RT-QC for the chlorophyll
# in the framework of Bioargo following the ADMT 14. on the OAO platform
#
# Based on Pabim White Book, Xing 2012 (NPQ)
#
# C. Schmechtig  February 2014
# -adjusted in August 2014 in order to perform Adjusted CHLA no longer on request but automatically
# -adjusted in October 2014 in order to calculate BBP and to perform the RT_QC on it 
# -calculate the quenching in 90% in the MLD
# -Adjusted in February 2015 in order to perform the RT_QC on radiometry with emanuele Organelli's work  
# -Change the CHLA_ADJUSTED_QC in NPQ zone to 8 (according to coriolis recommendations)
# -Adjusted in december 2015 to perform RTQC on CDOM, small corrections on BBP QC, and on radiometry QC
# -Adjusted in March 2016 to be able to finish the QC processing even if the CTD is not working properly and restore the spike test for BBP 
# -Change in May 2017 to adapt to RNetCDF library 
#
# The values of DEPTH_LIMIT and DELTA_DEPTH
# 
##############################################################################


##### loading RNetCDF library 

library(ncdf4)

##### Require oceanographic routines to perform MLD 

require(oce)

##### Require some Filter function from Xiao Gang  
source("/home/bioargo/PROGRAM_TEST/RT_QC/RunningFilter.R")

##### Require the calculation of betasw from zhang et al 
source("betasw124_ZHH2009.R")

##### Require Emanuele's Routine for RT_QC radiometry
source("RT_QC_radiometry_function_oao_2.R")

##### Some tests on the solar angle for the Radiometry
source("julian_to_date.R")

#### get the mission name 

uf=commandArgs()

mission  <- uf[2]

mis=substr(mission,start=7,stop=10)

##### Set up some threshold values on MLD 

MLD_LIMIT=0.03

DEPTH_LIMIT=200

DELTA_DEPTH=50

#### Global Range Values 

MIN_RANGE=-0.1

MAX_RANGE=50

#### FOR CDOM

MIN_RANGE_CDOM=-0.5

MAX_RANGE_CDOM=375

# FOR BBP

MIN_RANGE_BBP532=-0.000025

MIN_RANGE_BBP700=-0.000005

MAX_RANGE_BBP=0.1

#### Threshold for potential density 

POTDENS_THR=0.03

#### mission with 2 BB for Eaims and UK bio-Argo

mission_type=substr(mission,1,3)

if ( mission_type == "met" | mission_type == "bas" | mission_type == "imr")
{
mission_2BB=TRUE
} else {
mission_2BB=FALSE
}

#### Factory Calibration Coefficient !!!!! Pick Up in the File !!!!

file_calib=paste("/home/bioargo/PROGRAM_TEST/PROFIL_NC/PREPROCESSING/TECHNIQUE/calibration_CHLA_",mission,".txt",sep="")

calib=read.table(file=file_calib,header=FALSE)

DARK_CHLA=calib$V1     

SCALE_CHLA=calib$V2

if(!mission_2BB)
{

#### Factory Calibration Coefficient !!!!! Pick Up in the File !!!!

file_calib_CDOM=paste("/home/bioargo/PROGRAM_TEST/PROFIL_NC/PREPROCESSING/TECHNIQUE/calibration_CDOM_",mission,".txt",sep="")

calib_CDOM=read.table(file=file_calib_CDOM,header=FALSE)

DARK_CDOM=calib_CDOM$V1     

SCALE_CDOM=calib_CDOM$V2

}

#### Creating the list of files to RT_QC 

repNCDF=paste("/var/www/oao/BD_FLOAT/NETCDF/",mission,"/",sep="") # Directory

#LIST_nc=list.files(repNCDF,pattern=mission) # all Files

liste_to_do_file=paste(repNCDF,"liste_to_do",sep="") 

################################
# testing the size of the file"

fileinfo=file.info(liste_to_do_file)

if(fileinfo$size != 0){

liste_to_do=read.table(liste_to_do_file)

} else {

stop("RT_QC already done")

}
###############################

LIST_nc=liste_to_do$V1

####  Init Hist_CHLA Value

last_file=paste(repNCDF,"lastNC_",mis,sep="")

last_calib=paste(repNCDF,"lastDARK_",mis,sep="")

if(!mission_2BB)
{

if(file.exists(c(last_calib))) {

lastDARK=read.table(last_calib)
HIST_CHLA=lastDARK$V1
FOND_CDOM_HIST=lastDARK$V2
FIRST_CDOM_FLAG=lastDARK$V3
SHIFT_CDOM=lastDARK$V4

} else {

HIST_CHLA=DARK_CHLA
FOND_CDOM_HIST=0
FIRST_CDOM_FLAG=FALSE
SHIFT_CDOM=0
}

} else {

if(file.exists(c(last_calib))) {

lastDARK=read.table(last_calib)
HIST_CHLA=lastDARK$V1
print(HIST_CHLA)

} else {

HIST_CHLA=DARK_CHLA

}
}

####################################################################
## Beginning of the loop on every profile
####################################################################

for (IDnc in LIST_nc) {

# Get the name of the file without extension

file=substr(IDnc,40,(nchar(IDnc)-3))

# Fichier de sortie du profil PRES,CHLA,CHLA_ADJUSTED et Flags

path_out_peigne=paste("/var/www/oao/BD_FLOAT/NETCDF/PEIGNE/",file,"_prof.txt",sep="")

# Fichier Netcdf

file_in=paste(IDnc,sep="")

filenc=nc_open(file_in,readunlim=FALSE,write=TRUE)

##################################################
#### 1. Get Data from the Netcdf File 
##################################################

TEMP=ncvar_get(filenc,"TEMP")

PSAL=ncvar_get(filenc,"PSAL")

PRES=ncvar_get(filenc,"PRES")

CHLA=ncvar_get(filenc,"CHLA")

CHL_RAW=ncvar_get(filenc,"CHL_RAW")

CHLA_QC=ncvar_get(filenc,"CHLA_QC")

CHLA_ADJUSTED=ncvar_get(filenc,"CHLA_ADJUSTED")

CHLA_ADJUSTED_QC=ncvar_get(filenc,"CHLA_ADJUSTED_QC")

BACKSCATTERING=ncvar_get(filenc,"BACKSCATTERING")

BACKSCATTERING_QC=ncvar_get(filenc,"BACKSCATTERING_QC")

CDOM_RAW=ncvar_get(filenc,"CDOM_RAW")

CDOM_ADJUSTED=ncvar_get(filenc,"CDOM_ADJUSTED")

CDOM_ADJUSTED_QC=ncvar_get(filenc,"CDOM_ADJUSTED_QC")

CDOM=ncvar_get(filenc,"CDOM")

CDOM_QC=ncvar_get(filenc,"CDOM_QC")

# in case of 2BB mission the BB at 532 is stored in the CDOM
if(mission_2BB)
{
BBP532=CDOM
BBP532_QC=CDOM_QC
BETASW532=CDOM
}

## Initializing BBP to BACKSCATTERING and also BETASW,to have the right dimension 

BBP700=BACKSCATTERING

#BBP700_ADJUSTED=BACKSCATTERING

BBP700_QC=BACKSCATTERING_QC

#BBP700_ADJUSTED_QC=CHLA_QC

BETASW700=BACKSCATTERING

#### Radiometry

JULD=ncvar_get(filenc,"JULD") # to test solar angle

LATITUDE=ncvar_get(filenc,"LATITUDE") # to test solar angle

LONGITUDE=ncvar_get(filenc,"LONGITUDE") # to test solar angle

IRR_380=ncvar_get(filenc,"DOWNWELLING_IRRADIANCE_380") #irradiance at 380 nm

DOWNWELLING_IRRADIANCE_380_QC=ncvar_get(filenc,"DOWNWELLING_IRRADIANCE_380_QC")

IRR_412=ncvar_get(filenc,"DOWNWELLING_IRRADIANCE_412") #irradiance at 412 nm

DOWNWELLING_IRRADIANCE_412_QC=ncvar_get(filenc,"DOWNWELLING_IRRADIANCE_412_QC")

IRR_490=ncvar_get(filenc,"DOWNWELLING_IRRADIANCE_490") #irradiance at 490 nm

DOWNWELLING_IRRADIANCE_490_QC=ncvar_get(filenc,"DOWNWELLING_IRRADIANCE_490_QC")

PAR=ncvar_get(filenc,"PAR") # Photosynthetically Available Radiation (PAR)

PAR_QC=ncvar_get(filenc,"PAR_QC") # Photosynthetically Available Radiation (PAR)

ind_norad=(IRR_380>=99999.)

ind_rad=which(IRR_380!=99999.)

test_radiometry=solar_angle_test(JULD,LATITUDE,LONGITUDE)

##################################################
#### 1.b Initialising the QC for ACRI stat 
##################################################

# Global Range
CHLA_QC_GR=rep(0,length(PRES))

# Negative Spike Test
CHLA_QC_ST=rep(0,length(PRES))

# Mixed layer
CHLA_QC_ML=rep(0,length(PRES))

# Quenching Test
CHLA_QC_QT=rep(0,length(PRES))

# Density Inversion
CHLA_QC_DI=rep(0,length(PRES))

# STRONG Adjustment
CHLA_QC_SA=rep(0,length(PRES))

# Global range BBP700
BBP700_QC_GR=rep(0,length(PRES))

# Global range BBP532
BBP532_QC_GR=rep(0,length(PRES))

# Negative Spike Test BBP700 
BBP700_QC_ST=rep(0,length(PRES))

# Negative Spike Test BBP532
BBP532_QC_ST=rep(0,length(PRES))

# Type radiometry 380 
Type_380=rep(0,length(PRES))

# Type radiometry 412 
Type_412=rep(0,length(PRES))

# Type radiometry 490 
Type_490=rep(0,length(PRES))

# Type radiometry PAR 
Type_PAR=rep(0,length(PRES))

##################################################################
# 2. Calculate MLD
# We only take care of Pressure, associated with CTD Measurements
###################################################################

ind_psal=which(PSAL!=99999.)

PSAL[which(PSAL>=99999.)]=NA

TEMP[which(TEMP>=99999.)]=NA

PRES_CTD=PRES[ind_psal]

TEMP_CTD=TEMP[ind_psal]

PSAL_CTD=PSAL[ind_psal]

# calcul de la temperature potentielle 

THETA=swTheta(PSAL_CTD,TEMP_CTD,PRES_CTD)

# calcul de l anomalie de densite potentielle
 
POTDENS=swSigmaTheta(PSAL_CTD,TEMP_CTD,PRES_CTD)

#######################################################
# 2.1        TESTING DENSITY INVERSION (Test. 14) 
#######################################################

FLAG_BAD_POTDENS=rep(FALSE,length(PRES_CTD))

if (length(PRES_CTD)>=5) {

for(i in 2 : length(PRES_CTD)) {
        
	if(POTDENS[i]<=(POTDENS[i-1]-0.03)) FLAG_BAD_POTDENS[i-1]=TRUE 
}

# without forgetting the last value

if(POTDENS[length(PRES_CTD)]<=(POTDENS[length(PRES_CTD)-1]-0.03)) FLAG_BAD_POTDENS[length(PRES_CTD)]=TRUE

}

CHLA_QC_DI[ind_psal[FLAG_BAD_POTDENS]]=1

#########################################################
# 2.2  	Calculating the MLD 
#########################################################

POTDENS[which(FLAG_BAD_POTDENS==TRUE)]=NA

ind_pot=which(FLAG_BAD_POTDENS==FALSE)

PRES_POT=PRES_CTD[ind_pot]

POTDENS_POT=POTDENS[ind_pot]

abs_pres_10=abs(PRES_POT-10)

# on calcule la densite potentielle a 10 m

POTDENS_10=max(POTDENS_POT[which(abs_pres_10==min(abs_pres_10))])

# index pour trouver la profondeur de la MLD
# on initialise au max de profondeur au cas ou on ne trouverait pas de MLD

MLD=max(PRES_POT)

if (length(PRES_CTD)>=5) {

MLD_CALC=(POTDENS_POT-POTDENS_10)

for(i in 1 : length(PRES_POT)) {
 
        if (MLD_CALC[i]>MLD_LIMIT){ 

	MLD=PRES_POT[i]
	break

	}

}
}

###################################
# evaluate the quenching in 90% MLD
###################################
MLD_NPQ=0.9*MLD

###################################################################
# 3. Value at Depth
###################################################################

# TEMPORARY
# define the index for which chlorophyll is correct or not !!!TEMPORARY!!!

ind_chl=which(CHL_RAW!=99999.)

CHL_RAW_CHL=CHL_RAW[ind_chl]

ind_nochl=which(CHL_RAW>=99999.)

CHL_RAW[ind_nochl]=NA

CHLA[ind_nochl]=NA

CHLA_ADJUSTED[ind_nochl]=NA

CHLA_QC[ind_chl]=1

CHLA_QC[ind_nochl]=" "

CHLA_ADJUSTED_QC[ind_chl]=1

CHLA_ADJUSTED_QC[ind_nochl]=" "

PRES_CHL=PRES[ind_chl]

CDOM_RAW[ind_nochl]=NA

CDOM_CHL=CDOM[ind_chl]

CDOM[ind_nochl]=NA

CDOM_ADJUSTED[ind_nochl]=NA

CDOM_QC[ind_chl]=1

CDOM_QC[ind_nochl]=" "

CDOM_ADJUSTED_QC[ind_chl]=1

CDOM_ADJUSTED_QC[ind_nochl]=" "

# END OF TEMPORARY

FLAG_MLD=FALSE # Not mixed
	
FLAG_NEWDARK=FALSE # No New Dark

FLAG_SP=FALSE # No Small pressure Profile

FLAG_STRONG=FALSE # 

NEW_DARK=HIST_CHLA 

###################################################################################
# if the calibration changed from the factory calibration the CHLA_QC is set to 3 
# and we try to estimate a new value at depth, the 3 values is for the whole profile
###################################################################################

if(HIST_CHLA!=DARK_CHLA) CHLA_QC[ind_chl]=3

###################################################################################

if (length(PRES_CHL)>=5) {

if (max(PRES_CHL)<(MLD+DEPTH_LIMIT+DELTA_DEPTH)){

FLAG_SP=TRUE
FLAG_MLD=TRUE
CHLA_ADJUSTED=SCALE_CHLA*(CHL_RAW-HIST_CHLA)
CHLA_ADJUSTED_QC[ind_chl]=2 # No deep value estimated because of shallow profile so the adjusted value might be good 
CHLA_QC_ML[ind_chl]=1

} else { # estimate a new deep value

ind_fond=which((PRES_CHL>=(max(PRES_CHL)-DELTA_DEPTH)))

NEW_DARK=round(median(CHL_RAW_CHL[ind_fond]),0)

if(NEW_DARK!=HIST_CHLA)FLAG_NEWDARK=TRUE

if(abs(NEW_DARK-DARK_CHLA)>0.2*DARK_CHLA){

FLAG_STRONG=TRUE
CHLA_QC_SA[ind_chl]=1
CHLA_QC[ind_chl]=3 # new dark value so the CHLOROPHYLL might be bad
CHLA_ADJUSTED=SCALE_CHLA*(CHL_RAW-HIST_CHLA)
CHLA_ADJUSTED_QC[ind_chl]=3 # Too strong variation for the chlorophyll adjustment

} else {

CHLA_QC[ind_chl]=3 # new dark value so the CHLOROPHYLL might be bad
CHLA_ADJUSTED=SCALE_CHLA*(CHL_RAW-NEW_DARK)
CHLA_ADJUSTED_QC[ind_chl]=1 # new dark value so the CHLOROPHYLL ADJUSTED might be good
HIST_CHLA=NEW_DARK

} # fin difference sup 20%

} # fin else shallow profile or MLD

}

###################################################################################
# FOR CDOM we should stay at the same value at depth if the water mass doesn't change
###################################################################################
if(!mission_2BB)
{

if (length(PRES_CHL)>=5) {

if (max(PRES_CHL)<(MLD+DEPTH_LIMIT+DELTA_DEPTH)){

CDOM_ADJUSTED=SCALE_CDOM*(CDOM_RAW-DARK_CDOM)+SHIFT_CDOM
CDOM_ADJUSTED_QC[ind_chl]=2 # No deep value estimated because of shallow profile so the adjusted value might be good 

} else { # estimate a new deep value

ind_fond=which((PRES_CHL>=(max(PRES_CHL)-DELTA_DEPTH)))

FOND_CDOM=median(CDOM_CHL[ind_fond])

if(!FIRST_CDOM_FLAG){
FIRST_CDOM_FLAG=TRUE
FOND_CDOM_HIST=FOND_CDOM
}
# on calcule le shift par rapport a la valeur de fond
if(FOND_CDOM!=FOND_CDOM_HIST){

SHIFT_CDOM=FOND_CDOM_HIST-FOND_CDOM
CDOM_QC[ind_chl]=3
}

CDOM_ADJUSTED=SCALE_CDOM*(CDOM_RAW-DARK_CDOM)+SHIFT_CDOM

CDOM_ADJUSTED_QC[ind_chl]=1 # value adjusted might be good 

} # fin else shallow profile or MLD for CDOM

} # fin length PRES_CHL

new_calibration_cdom=paste("Shift CDOM= ",SHIFT_CDOM," FOND_CDOM_HIST= ",FOND_CDOM_HIST,sep="")

##### And setting up the new value 
CDOM_CHL=CDOM[ind_chl]
CDOM_FILTERED_CHL=CDOM[ind_chl]
CDOM_ADJUSTED_CHL=CDOM_ADJUSTED[ind_chl]
CDOM_ADJUSTED_FILTERED_CHL=CDOM_ADJUSTED[ind_chl]

}

new_calibration=paste("Factory Calibration, Dark= ",DARK_CHLA,", Actual calibration, Dark= ",HIST_CHLA,sep="")

CHLA_ADJUSTED_CHL=CHLA_ADJUSTED[ind_chl]
CHLA_ADJUSTED_FILTERED_CHL=CHLA_ADJUSTED[ind_chl]

###################################################
### 4. NPQ correction       #######################
###################################################

SURF_CHLA_NPQ=0
PRES_NPQ=0
CHLA_NPQ=0
SURF_NORM=0

###### Before anything we must filter the data

#no sliding median if too few values
if(length(ind_chl)>=5)
{ 

MED_CHL=RunningFilter(2,CHLA_ADJUSTED_CHL,na.fill=T, ends.fill=T, Method="Median")

} else {

MED_CHL=CHLA_ADJUSTED_CHL

}

CHLA_ADJUSTED_FILTERED_CHL=MED_CHL

RES_CHLA_ADJUSTED_CHL=CHLA_ADJUSTED_CHL-MED_CHL

Q10=rep(2*quantile(RES_CHLA_ADJUSTED_CHL,0.10),length(PRES_CHL))

Q90=rep(2*quantile(RES_CHLA_ADJUSTED_CHL,0.90),length(PRES_CHL))

CHLA_ADJUSTED_FILTERED_CHL[which(RES_CHLA_ADJUSTED_CHL>Q90)]=NA

if(!FLAG_MLD){

ind_MLD=which(PRES_CHL<=MLD_NPQ)
CHLA_NPQ=max(CHLA_ADJUSTED_FILTERED_CHL[ind_MLD],na.rm=TRUE)
#print(CHLA_NPQ)

# index du max de chloro
i_NPQ=which.max(CHLA_ADJUSTED_FILTERED_CHL[ind_MLD])

PRES_CHL_NPQ=PRES_CHL[ind_MLD]
# Pression du max de chloro
PRES_NPQ=PRES_CHL_NPQ[i_NPQ]

# on corrige la chloro depuis la valeur du max a la surface
CHLA_ADJUSTED[which(PRES<=PRES_NPQ)]=CHLA_NPQ

CHLA_ADJUSTED[ind_nochl]=NA

# on affecte les QC 

if( !FLAG_STRONG & !FLAG_SP ) CHLA_ADJUSTED_QC[which(PRES<=PRES_NPQ)]=8

CHLA_QC[which(PRES<=PRES_NPQ)]=3

CHLA_QC[ind_nochl]=" "

CHLA_ADJUSTED_QC[ind_nochl]=" " 

CHLA_QC_QT[which(PRES<=PRES_NPQ)]=1

CHLA_QC_QT[ind_nochl]=0

}

CHLA_ADJUSTED_NONPQ=CHLA_ADJUSTED[ind_chl]

###################################################
### 5. Test 6, Global Range #######################
###################################################

ind_out_of_range=which((CHLA>MAX_RANGE) | (CHLA<MIN_RANGE))

CHLA_QC[ind_out_of_range]=4

ind_out_of_range_adjusted=which((CHLA_ADJUSTED>MAX_RANGE) | (CHLA_ADJUSTED<MIN_RANGE))

CHLA_ADJUSTED_QC[ind_out_of_range_adjusted]=4

CHLA_QC_GR[ind_out_of_range_adjusted]=1

###################################################
### 5. Test 6, Global Range CDOM   ################
###################################################
if(!mission_2BB)
{

ind_out_of_range_CDOM=which((CDOM>MAX_RANGE_CDOM) | (CDOM<MIN_RANGE_CDOM))

CDOM_QC[ind_out_of_range_CDOM]=4

ind_out_of_range_CDOM_adjusted=which((CDOM_ADJUSTED>MAX_RANGE_CDOM) | (CDOM_ADJUSTED<MIN_RANGE_CDOM))

CDOM_ADJUSTED_QC[ind_out_of_range_CDOM_adjusted]=4

}
###################################################
### 6. Test 9, Spike Test   #######################
###################################################

spike=rep(FALSE,length(PRES_CHL))

if(!FLAG_MLD){
spike[which((RES_CHLA_ADJUSTED_CHL<Q10) & (PRES_CHL > MLD))]=TRUE
} else {
spike[which(RES_CHLA_ADJUSTED_CHL<Q10)]=TRUE
}

ind_spike=ind_chl[spike]
CHLA_QC[ind_spike]=4
CHLA_ADJUSTED_QC[ind_spike]=4
CHLA_QC_ST[ind_spike]=1

###################################################
### 6. Test 9, Spike Test   CDOM ##################
###################################################
if(!mission_2BB)
{

spike1_cdom=rep(FALSE,length(PRES_CHL))

spike2_cdom=rep(FALSE,length(PRES_CHL))

# Test 6.1

CDOM_RAW_CHL=CDOM_RAW[ind_chl]

CDOM_Q25=quantile(CDOM_RAW_CHL,0.25)

CDOM_Q75=quantile(CDOM_RAW_CHL,0.75)

IQR=(CDOM_Q75-CDOM_Q25)

spike1_cdom[which((CDOM_RAW_CHL>CDOM_Q75+2*IQR)|(CDOM_RAW_CHL<CDOM_Q25-2*IQR))]=TRUE

ind_spike1_cdom=ind_chl[spike1_cdom]

# Test 6.2

#no sliding median if too few values
if(length(ind_chl)>=51)
{ 

MED_CDOM_smooth=RunningFilter(25,CDOM_RAW_CHL,na.fill=T, ends.fill=T, Method="Average")

} else {

MED_CDOM_smooth=CDOM_RAW_CHL

}

spike2_cdom[which(abs(CDOM_RAW_CHL-MED_CDOM_smooth)>4)]=TRUE

ind_spike_cdom=ind_chl[spike1_cdom | spike2_cdom]

CDOM_QC[ind_spike_cdom]=4

CDOM_ADJUSTED_QC[ind_spike_cdom]=4

}
##################################################################
### 7. Additional Test : Pressure monotonically increasing
##################################################################

#FLAG_BAD_PRES=rep(FALSE,length(PRES_CHL))

#for(i in 2 : length(PRES_CHL)) {

#        if((PRES_CHL[i]<=PRES_CHL[i-1])) FLAG_BAD_PRES[i-1]=TRUE 
#}
# without forgetting the last value
#if(PRES_CHL[length(PRES_CHL)]<=PRES_CHL[length(PRES_CHL)-1]) FLAG_BAD_PRES[length(PRES_CHL)]=TRUE

#ind_bad_pres=ind_chl[FLAG_BAD_PRES]
#CHLA_QC[ind_bad_pres]=4
#CHLA_ADJUSTED_QC[ind_bad_pres]=4


###############################################################################
#######################################################################
# 7.a Work on BBP
#######################################################################

# Interpolate TEMP and PSAL through the whole profile
# and also extrapolate when needed rule=2
if (length(PRES_CTD)>=5)
{
TEMP_PRES=approx(PRES_CTD,TEMP_CTD,PRES,method="linear",rule=2)

TEMP_CHL=TEMP_PRES$y[ind_chl]

PSAL_PRES=approx(PRES_CTD,PSAL_CTD,PRES,method="linear",rule=2)

PSAL_CHL=PSAL_PRES$y[ind_chl]

# initialisation 
# TEST sur le nom de mission au cas ou il y aurait deux BBP
BETASW124_CHL=PRES_CHL

# BB a 700

lambda=700

# Must set up a loop to calculate Betasw

for(i in 1 : length(PRES_CHL)) {

BETASW124_CHL[i]=calc_betasw124(lambda,PSAL_CHL[i],TEMP_CHL[i])
}

BETASW700[ind_chl]=BETASW124_CHL

BETASW700[ind_nochl]=NA

BACKSCATTERING[ind_nochl]=NA

BBP700=2*pi*1.076*(BACKSCATTERING-BETASW700)

BBP700[ind_nochl]=NA

BBP700_CHL=BBP700[ind_chl]

##############################
# missions a deux BB
##############################

if(mission_2BB)
{

lambda=532

# Must set up a loop to calculate Betasw

for(i in 1 : length(PRES_CHL)) {

BETASW124_CHL[i]=calc_betasw124(lambda,PSAL_CHL[i],TEMP_CHL[i])
}

BETASW532[ind_chl]=BETASW124_CHL

BETASW532[ind_nochl]=NA

CDOM[ind_nochl]=NA

BBP532=2*pi*1.076*(CDOM-BETASW532)

BBP532[ind_nochl]=NA

BBP532_CHL=BBP532[ind_chl]

#######################################################################
# PERFORM RT_QC for BBP532
#######################################################################
# a priori data are considered as correct

BBP532_QC[ind_chl]=1
BBP532_QC[ind_nochl]=" "

#### Global Range test ####

ind_out_of_range_BBP532=which((BBP532>MAX_RANGE_BBP) | (BBP532<MIN_RANGE_BBP532))

BBP532_QC[ind_out_of_range_BBP532]=3

BBP532_QC_GR[ind_out_of_range_BBP532]=1

#### BBP bad Offset and spike test ####

###### Before anything we must filter the data (sliding Median on 5 values)
if(length(ind_chl)>=5)
{

BBP532_MED_CHL=RunningFilter(2,BBP532_CHL,na.fill=T, ends.fill=T, Method="Median")

} else {

BBP532_MED_CHL=BBP532_CHL

}

# Spike test
RES_BBP532_CHL=BBP532_CHL-BBP532_MED_CHL

BBP532_Q10=rep(2*quantile(RES_BBP532_CHL,0.10),length(PRES_CHL))

spike_BBP532=rep(FALSE,length(PRES_CHL))

spike_BBP532[which(RES_BBP532_CHL<BBP532_Q10)]=TRUE

ind_spike_BBP532=ind_chl[spike_BBP532]

BBP532_QC[ind_spike_BBP532]=4

BBP532_QC_ST[ind_spike_BBP532]=1

# Bad Offset Test

BAD_OFFSET_532=-20*min(BBP532_MED_CHL)

BO_BBP532_H=rep(FALSE,length(PRES_CHL))
BO_BBP532_M=rep(FALSE,length(PRES_CHL))

ind_thresh_BBP532=which(BBP532_MED_CHL<MIN_RANGE_BBP532)

if(length(ind_thresh_BBP532)>0 ){
if(BBP532_CHL[ind_thresh_BBP532]>=BAD_OFFSET_532)
{
BO_BBP532_M[ind_thresh_BBP532]=TRUE
} else {
print("ok")
BO_BBP532_H[ind_thresh_BBP532]=TRUE
}
}

ind_BO_BBP532_M=ind_chl[BO_BBP532_M]

ind_BO_BBP532_H=ind_chl[BO_BBP532_H]

if(length(ind_BO_BBP532_H)>0){
if(BBP532_QC[ind_BO_BBP532_H]!=4)BBP532_QC[ind_BO_BBP532_H]=3
}

if(length(ind_BO_BBP532_M)>0){
if(BBP532_QC[ind_BO_BBP532_M]!=4)BBP532_QC[ind_BO_BBP532_M]=2
}

}# End mission 2BB

#######################################################################
# 7.b PERFORM RT_QC for BBP
#######################################################################
# a priori data are considered as correct

BBP700_QC[ind_chl]=1
BBP700_QC[ind_nochl]=" "

#BBP700_ADJUSTED_QC[ind_chl]=1
#BBP700_ADJUSTED_QC[ind_nochl]=" "

#### Global Range test ####

ind_out_of_range_BBP=which((BBP700>MAX_RANGE_BBP) | (BBP700<MIN_RANGE_BBP700))

BBP700_QC[ind_out_of_range_BBP]=3

BBP700_QC_GR[ind_out_of_range_BBP]=1

#### BBP bad Offset and spike test ####

###### Before anything we must filter the data (sliding Median on 5 values)
if(length(ind_chl)>=5)
{

BBP700_MED_CHL=RunningFilter(2,BBP700_CHL,na.fill=T, ends.fill=T, Method="Median")

} else {

BBP700_MED_CHL=BBP700_CHL

}

# Restore spike test
# Lets calculate the residual between the signal and the sliding median
RES_BBP700_CHL=BBP700_CHL-BBP700_MED_CHL

BBP700_Q10=rep(2*quantile(RES_BBP700_CHL,0.10),length(PRES_CHL))

spike_BBP700=rep(FALSE,length(PRES_CHL))

spike_BBP700[which(RES_BBP700_CHL<BBP700_Q10)]=TRUE

ind_spike_BBP700=ind_chl[spike_BBP700]

BBP700_QC[ind_spike_BBP700]=4

BBP700_QC_ST[ind_spike_BBP700]=1

## Bad Offset test

BAD_OFFSET_700=-20*min(BBP700_MED_CHL)

BO_BBP700_H=rep(FALSE,length(PRES_CHL))
BO_BBP700_M=rep(FALSE,length(PRES_CHL))

ind_thresh_BBP700=which(BBP700_MED_CHL<MIN_RANGE_BBP700)

if(length(ind_thresh_BBP700)>0 ){
if(BBP700_CHL[ind_thresh_BBP700]>=BAD_OFFSET_700)
{
BO_BBP700_M[ind_thresh_BBP700]=TRUE
} else {
print("ok")
BO_BBP700_H[ind_thresh_BBP700]=TRUE
}
}

ind_BO_BBP700_M=ind_chl[BO_BBP700_M]

ind_BO_BBP700_H=ind_chl[BO_BBP700_H]

if(length(ind_BO_BBP700_H)>0){
if(BBP700_QC[ind_BO_BBP700_H]!=4)BBP700_QC[ind_BO_BBP700_H]=3
}

if(length(ind_BO_BBP700_M)>0){
if(BBP700_QC[ind_BO_BBP700_M]!=4)BBP700_QC[ind_BO_BBP700_M]=2
}

} else {

BBP700_QC[ind_chl]=4

if(mission_2BB)

{
BBP532_QC[ind_chl]=4
}

}# end on length of the CTD Values

#######################################################################
# 7.c Define new variables to store the BBP
#######################################################################
# first we test if it exists or not, if not then we create !!


if (is.null(filenc$var$BBP700$name))
{
# get the dimension of the array
dim_nprof=filenc$dim[["N_PROF"]]

dim_nlev=filenc$dim[["N_LEVELS"]]

# define Missvalue
mv=99.99

#define the new variables
var_bbp700=ncvar_def( "BBP700","m-1", list(dim_nlev,dim_nprof), mv,prec="float")
var_bbp700_adjusted=ncvar_def( "BBP700_ADJUSTED","m-1", list(dim_nlev,dim_nprof),mv,prec="float")

var_bbp700_adjusted_qc=ncvar_def( "BBP700_ADJUSTED_QC" ,units=" ", dim=list(dim_nlev,dim_nprof),missval=" ",prec="char")
var_bbp700_qc=ncvar_def( "BBP700_QC" , units=" ", dim=list(dim_nlev,dim_nprof),missval=" ",prec="char")


filenc=ncvar_add( filenc, var_bbp700)
filenc=ncvar_add( filenc, var_bbp700_adjusted)
filenc=ncvar_add( filenc, var_bbp700_adjusted_qc)
filenc=ncvar_add( filenc, var_bbp700_qc)

ncatt_put( filenc, "BBP700","_FillValue",99.99,prec="float")
#ncatt_put( filenc, "BBP700_ADJUSTED","_FillValue",99.99)

if(mission_2BB)
{
var_bbp532=ncvar_def( "BBP532","m-1", list(dim_nlev,dim_nprof), mv)
var_bbp532_adjusted=ncvar_def( "BBP532_ADJUSTED","m-1", list(dim_nlev,dim_nprof),mv)

var_bbp532_adjusted_qc=ncvar_def( "BBP532_ADJUSTED_QC" ,units=" ", dim=list(dim_nlev,dim_nprof),missval=" ")
var_bbp532_qc=ncvar_def( "BBP532_QC" , units=" ", dim=list(dim_nlev,dim_nprof),missval=" ")


filenc=ncvar_add( filenc, var_bbp532)
filenc=ncvar_add( filenc, var_bbp532_adjusted)
filenc=ncvar_add( filenc, var_bbp532_adjusted_qc)
filenc=ncvar_add( filenc, var_bbp532_qc)

#ncatt_put( filenc, "BBP532","_FillValue",99.99,prec="float")
#ncatt_put( filenc, "BBP532_ADJUSTED","_FillValue",99.99)

}# end mission 2BB

}# end definition new variables

ncatt_put( filenc, "BBP700","_FillValue",99.99,prec="float")
ncvar_put( filenc, "BBP700_QC", paste(BBP700_QC,collapse=""))
#ncvar_put( filenc, "BBP700_ADJUSTED_QC", paste(BBP700_ADJUSTED_QC,collapse=""))
ncvar_put( filenc, "BBP700", BBP700)
#ncvar_put( filenc, "BBP700_ADJUSTED", BBP700_ADJUSTED)

if(mission_2BB)
{
ncvar_put( filenc, "BBP532_QC", paste(BBP532_QC,collapse=""))

ncvar_put( filenc, "BBP532", BBP532)

}

###################################################
### 7.d RT_QC on radiometry #######################
###################################################
if(test_radiometry){

if(length(ind_rad) <= 5){ 
DOWNWELLING_IRRADIANCE_380_QC[ind_rad]="3"
DOWNWELLING_IRRADIANCE_412_QC[ind_rad]="3"
DOWNWELLING_IRRADIANCE_490_QC[ind_rad]="3"
PAR_QC[ind_rad]="3"
Type_380=rep(3,length(PRES))
Type_412=rep(3,length(PRES))
Type_490=rep(3,length(PRES))
Type_PAR=rep(3,length(PRES))

} else {

RADIOMETRY_QC <- RT_QC_radiometry(PRES,IRR_380,IRR_412,IRR_490,PAR)

# affect QC and remoce NA
DOWNWELLING_IRRADIANCE_380_QC=RADIOMETRY_QC[[1]]
# affect QC and remoce NA
DOWNWELLING_IRRADIANCE_412_QC=RADIOMETRY_QC[[2]]
# affect QC and remoce NA
DOWNWELLING_IRRADIANCE_490_QC=RADIOMETRY_QC[[3]]
# affect QC and remoce NA
PAR_QC=RADIOMETRY_QC[[4]]

# profile type
Type_380=RADIOMETRY_QC[[5]]
# profile type
Type_412=RADIOMETRY_QC[[6]]
# profile type
Type_490=RADIOMETRY_QC[[7]]
# profile type
Type_PAR=RADIOMETRY_QC[[8]]
}

} else {

DOWNWELLING_IRRADIANCE_380_QC[ind_rad]="3"
DOWNWELLING_IRRADIANCE_412_QC[ind_rad]="3"
DOWNWELLING_IRRADIANCE_490_QC[ind_rad]="3"
PAR_QC[ind_rad]="3"
Type_380=rep(3,length(PRES))
Type_412=rep(3,length(PRES))
Type_490=rep(3,length(PRES))
Type_PAR=rep(3,length(PRES))
print("pas routine") 
}

DOWNWELLING_IRRADIANCE_380_QC[ind_norad]=" "

DOWNWELLING_IRRADIANCE_412_QC[ind_norad]=" "

DOWNWELLING_IRRADIANCE_490_QC[ind_norad]=" "

PAR_QC[ind_norad]=" "


###################################################
### 8. Writing to File      #######################
###################################################

CHLA_ADJUSTED[ind_nochl]=99.99

ncatt_put( filenc, "CHLA_ADJUSTED","_FillValue",99.99,prec="float") 

ncvar_put( filenc, "CHLA_QC", paste(CHLA_QC,collapse=""))

ncvar_put( filenc, "CHLA_ADJUSTED_QC", paste(CHLA_ADJUSTED_QC,collapse=""))

ncvar_put( filenc, "CHLA_ADJUSTED", CHLA_ADJUSTED)

# QC for CDOM
if(!mission_2BB)
{
ncvar_put( filenc, "CDOM_QC", paste(CDOM_QC,collapse=""))

ncvar_put( filenc, "CDOM_ADJUSTED_QC", paste(CDOM_ADJUSTED_QC,collapse=""))

ncvar_put( filenc, "CDOM_ADJUSTED", CDOM_ADJUSTED)

}
# QC for Radiometry

ncvar_put( filenc, "DOWNWELLING_IRRADIANCE_380_QC",paste(DOWNWELLING_IRRADIANCE_380_QC,collapse=""))

ncvar_put( filenc, "DOWNWELLING_IRRADIANCE_412_QC",paste(DOWNWELLING_IRRADIANCE_412_QC,collapse=""))

ncvar_put( filenc, "DOWNWELLING_IRRADIANCE_490_QC",paste(DOWNWELLING_IRRADIANCE_490_QC,collapse=""))

ncvar_put( filenc, "PAR_QC",paste(PAR_QC,collapse=""))

# QC for Radiometry Profiles  

ncvar_put( filenc, "PROFILE_DOWNWELLING_IRRADIANCE_380_QC",paste(Type_380[1],collapse=""))

ncvar_put( filenc, "PROFILE_DOWNWELLING_IRRADIANCE_412_QC",paste(Type_412[1],collapse=""))

ncvar_put( filenc, "PROFILE_DOWNWELLING_IRRADIANCE_490_QC",paste(Type_490[1],collapse=""))

ncvar_put( filenc, "PROFILE_PAR_QC",paste(Type_PAR[1],collapse=""))


# Equation for CHLA ADJUSTMENT
if(!mission_2BB)
{
ncvar_put( filenc, "SCIENTIFIC_CALIB_EQUATION",start=c(1,1,1,1),count=c(256,2,1,1),c(new_calibration,new_calibration_cdom))
} else {
ncvar_put( filenc, "SCIENTIFIC_CALIB_EQUATION",count=c(256,1,1,1),new_calibration)
}
nc_close(filenc)

##### Writing the Peigne for ACRI_ST#################

profil=data.frame(PRES,CHLA,CHLA_QC,CHLA_ADJUSTED,CHLA_ADJUSTED_QC,CHLA_QC_GR,CHLA_QC_ST,CHLA_QC_ML,CHLA_QC_SA,CHLA_QC_DI,CHLA_QC_QT,BBP700_QC_GR,BBP700_QC_ST,BBP532_QC_GR,BBP532_QC_ST,Type_380,Type_412,Type_490,Type_PAR)

write.table(profil,row.names = T,file=path_out_peigne,quote=TRUE)

###################################################
### 9. Plotting ChLA, CHLA_ADJUSTED
### fichier de sortie pdf
###################################################
if(length(CHLA[ind_chl])!=0){

MAXCHLA=max(CHLA[ind_chl])

CHLA_ADJUSTED_BAD=CHLA_ADJUSTED[which(CHLA_ADJUSTED_QC==4)]
PRES_BAD=PRES[which(CHLA_ADJUSTED_QC==4)]

path_out_pdf=paste("/var/www/oao/BD_FLOAT/NETCDF/GRAPHE/",file,"_RT.pdf",sep="")

pdf(file=path_out_pdf)

matplot(CHLA[ind_chl],PRES[ind_chl],col=8,type="l",ylab="Depth [m]",xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,MAXCHLA+0.5),ylim=rev(c(0, max(PRES))))

matplot(CHLA_ADJUSTED[ind_chl],PRES[ind_chl],col=1,type="l",ylab="Depth [m]",xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,MAXCHLA+0.5),ylim=rev(c(0, max(PRES))),add=TRUE)

matplot(CHLA_ADJUSTED_BAD,PRES_BAD,type="p",pch=1, col=2,cex=2,add=TRUE)

legend("bottomright",c("Chl-a","Chl-a adjusted","Chl-a adjusted QC=4"),pch=c(20,20,20),col=c(8,1,2))
dev.off()

}
###################################################
### 9. Plotting CDOM, CDOM_ADJUSTED
### fichier de sortie pdf
###################################################
if(!mission_2BB)
{
if(length(CDOM[ind_chl])!=0){

MAXCDOM=max(CDOM[ind_chl])

CDOM_BAD=CDOM[which(CDOM_QC==4)]
PRES_BAD=PRES[which(CDOM_QC==4)]

path_out_pdf=paste("/var/www/oao/BD_FLOAT/NETCDF/GRAPHE/",file,"_CDOM_RT.pdf",sep="")

pdf(file=path_out_pdf)

matplot(CDOM[ind_chl],PRES[ind_chl],col=8,type="l",ylab="Depth [m]",xlab=expression("CDOM [ppb]"),xlim=c(-0.2,MAXCDOM+0.5),ylim=rev(c(0, max(PRES))))

matplot(CDOM_ADJUSTED[ind_chl],PRES[ind_chl],col=1,type="l",ylab="Depth [m]",xlab=expression("CDOM [ppb]"),xlim=c(-0.2,MAXCDOM+0.5),ylim=rev(c(0, max(PRES))),add=TRUE)

matplot(CDOM_BAD,PRES_BAD,type="p",pch=1, col=2,cex=2,add=TRUE)

legend("bottomright",c("CDOM","CDOM adjusted","CDOM QC=4"),pch=c(20,20,20),col=c(8,1,2))
dev.off()

}

}
###################################################
### 9.b Plotting BBP700, BBP700_QC
### fichier de sortie pdf
###################################################

if(length(BBP700[ind_chl])!=0){

MAXBBP700=max(BBP700[ind_chl])

BBP700_BAD=BBP700[which(BBP700_QC==4)]
PRES_BAD=PRES[which(BBP700_QC==4)]

path_out_pdf=paste("/var/www/oao/BD_FLOAT/NETCDF/GRAPHE/",file,"_BBP_RT.pdf",sep="")

pdf(file=path_out_pdf)

matplot(BBP700[ind_chl],PRES[ind_chl],col=8,type="l",ylab="Depth [m]",xlab=expression("BBP700 ["*m ^ -1 * "]"),xlim=c(0.,MAXBBP700),ylim=rev(c(0, max(PRES))))

matplot(BBP700_BAD,PRES_BAD,type="p",pch=1, col=2,cex=2,add=TRUE)

legend("bottomright",c("BBP700","BBP700 QC=4"),pch=c(20,20),col=c(8,2))

dev.off()

}

if(mission_2BB)
{

if(length(BBP532[ind_chl])!=0){

MAXBBP532=max(BBP532[ind_chl])

BBP532_BAD=BBP532[which(BBP532_QC==4)]
PRES_BAD=PRES[which(BBP532_QC==4)]

path_out_pdf=paste("/var/www/oao/BD_FLOAT/NETCDF/GRAPHE/",file,"_BBP532_RT.pdf",sep="")

pdf(file=path_out_pdf)

matplot(BBP532[ind_chl],PRES[ind_chl],col=8,type="l",ylab="Depth [m]",xlab=expression("BBP532 ["*m ^ -1 * "]"),xlim=c(0.,MAXBBP532),ylim=rev(c(0, max(PRES))))

matplot(BBP532_BAD,PRES_BAD,type="p",pch=1, col=2,cex=2,add=TRUE)

legend("bottomright",c("BBP532","BBP532 QC=4"),pch=c(20,20),col=c(8,2))

dev.off()
}# end mission 2BB

}

} 
# fin de boucle sur les fichiers 

#####################################################
# 10. Writing the name of the last controlled File  
# and last value of the HIST_CHLA
####################################################

write(IDnc,last_file)
if(!mission_2BB)
{

write(c(HIST_CHLA,FOND_CDOM_HIST,FIRST_CDOM_FLAG,SHIFT_CDOM),last_calib)

} else {

write(HIST_CHLA,last_calib)

}





