require(nortest)# Lilliefors test of normality

RT_QC_radiometry <- function(PRES,IRR_380,IRR_412,IRR_490,PAR) {

############################################################################################################################################################################################
# This program was designed in order to perform the RT-QC for radiometry profiles acquired by Bio-Argo floats
#
# Modified from the NASA Protocol vol.3 2003 
#
# Flag 1 to 4 are assigned according to the ARGO phylosophy during 4 steps 
#
# 1: DARK IDENTIFICATION   2a:MAJOR CLOUDS IDENTIFICATION  2b:SWINGING PROFILE IDENTIFICATION 3:TYPE PROFILE IDENTIFICATION 4a:TYPE2_FLAG ASSIGNEMENT 4b:TYPE1_FLAG ASSIGNEMENT
#
# E. Organelli  January 2015
# 
#############################################################################################################################################################################################

    TAB=cbind(PRES,IRR_380,IRR_412, IRR_490, PAR) #create a matrix using all the parameters 
    TAB=as.data.frame(TAB) #convert TAB in dataframe
    #head(TAB)
    #str(TAB)
  
    ##################################################
    #### 3. Create a data.frame without NA 
    ##################################################
    TAB_complete=TAB[complete.cases(TAB$IRR_380),]#to remove rows with NA. Same rows for any irradiance and PAR measurements.
    #str(TAB_complete)
    TAB_NA=TAB[!complete.cases(TAB$IRR_380),] #dataframe with only NA values. Useful when assignin flags.
  
    ##################################################
    #### 4. STEP 1: DARK IDENTIFICATION 
    ##################################################
  
    PARAM_NAMES = c("IRR_380", "IRR_412", "IRR_490", "PAR")
  
    #### 4.A Lilliefors Test for detecting lowest good irradiance measurement
    # alfa selected at 0.01
  
    Lilliefors_All = NULL
    for (param_name in PARAM_NAMES) {
        Lilliefors_param = rep(NA, length(TAB_complete$IRR_380)-4)
        for (l in 1 :(length(TAB_complete$IRR_380)-4)) {
            Lilliefors_param[l] = lillie.test(TAB_complete[[param_name]][ l:length(TAB_complete[[param_name]]) ])[2] #select p-value
        }
        Lilliefors_All[[param_name]] = unlist(Lilliefors_param)
    }

    limAll = NULL
    for (param_name in PARAM_NAMES) {
        i_param = which(abs(Lilliefors_All[[param_name]]) > 0.01)
        limAll[[param_name]] = i_param[1] -1
    }

    lim380=i_380[1]-1
    lim412=i_412[1]-1
    lim490=i_490[1]-1
    limPAR=i_PAR[1]-1

#### 4.B Check for negative values
if(!is.na(lim380)) {

OK380=TRUE
neg_380=which(TAB_complete$IRR_380[1:lim380]<=0)

{
  if(length(neg_380) !=0)
    lim380_bis=length(TAB_complete$IRR_380[1:neg_380[1]])-1
  else lim380_bis=lim380
}

} else {
	# 380
	OK380=FALSE 
        FLAG_380_QC=rep("3", length(TAB_complete$IRR_380)) #ASSIGNEMENT CORRESPONDING to STEP 3
        TAB_complete$FLAG_380_QC=FLAG_380_QC
        TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
        type380="3" #catherine
        TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
        #str(newdata_380)
}

if(!is.na(lim412)) {

OK412=TRUE
neg_412=which(TAB_complete$IRR_412[1:lim412]<=0)

{
  if(length(neg_412) !=0)
    lim412_bis=length(TAB_complete$IRR_412[1:neg_412[1]])-1
  else lim412_bis=lim412
}

} else {	
	# 412
	OK412=FALSE
        FLAG_412_QC=rep("3", length(TAB_complete$IRR_412)) #ASSIGNEMENT CORRESPONDING to STEP 3
        TAB_complete$FLAG_412_QC=FLAG_412_QC
        TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
        type412="3" #catherine
        TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
        #str(newdata_412)
}

if(!is.na(lim490)) {

OK490=TRUE
neg_490=which(TAB_complete$IRR_490[1:lim490]<=0)

{
  if(length(neg_490) !=0)
    lim490_bis=length(TAB_complete$IRR_490[1:neg_490[1]])-1
  else lim490_bis=lim490
}

} else {	
	# 490
	OK490=FALSE 
        FLAG_490_QC=rep("3", length(TAB_complete$IRR_490)) #ASSIGNEMENT CORRESPONDING to STEP 3
        TAB_complete$FLAG_490_QC=FLAG_490_QC
        TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
        type490="3" #catherine
        TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
        #str(newdata_490)
}

if(!is.na(limPAR)) {

OKPAR=TRUE
neg_PAR=which(TAB_complete$PAR[1:limPAR]<=0)

{
  if(length(neg_PAR) !=0)
    limPAR_bis=length(TAB_complete$PAR[1:neg_PAR[1]])-1
  else limPAR_bis=limPAR
} #logical test for negative PAR values 


} else {
	# PAR
	OKPAR=FALSE
        FLAG_PAR_QC=rep("3", length(TAB_complete$PAR)) #ASSIGNEMENT CORRESPONDING to STEP 3
        TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC
        TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
        typePAR="3" #catherine
        TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
        #str(newdata_PAR)
}

#### 4.C SD calculation for the removed part (Calculation just useful for statistics etc not for QC)
# dark_380=TAB_complete$IRR_380[(lim380_bis+1):length(TAB_complete$IRR_380)]
# sd_dark_380=sd(dark_380, na.rm=T)
# n_dark_380=length(dark_380)
# ave_dark_380=mean(dark_380, na.rm=T)
# CV_dark_380=(sd_dark_380/ave_dark_380)*100
# 
# dark_412=TAB_complete$IRR_412[(lim412_bis+1):length(TAB_complete$IRR_412)]
# sd_dark_412=sd(dark_412, na.rm=T)
# n_dark_412=length(dark_412)
# ave_dark_412=mean(dark_412, na.rm=T)
# CV_dark_412=(sd_dark_412/ave_dark_412)*100
# 
# dark_490=TAB_complete$IRR_490[(lim490_bis+1):length(TAB_complete$IRR_490)]
# sd_dark_490=sd(dark_490, na.rm=T)
# n_dark_490=length(dark_490)
# ave_dark_490=mean(dark_490, na.rm=T)
# CV_dark_490=(sd_dark_490/ave_dark_490)*100
# 
# dark_PAR=TAB_complete$PAR[(limPAR_bis+1):length(TAB_complete$PAR)]
# sd_dark_PAR=sd(dark_PAR, na.rm=T)
# n_dark_PAR=length(dark_PAR)
# ave_dark_PAR=mean(dark_PAR, na.rm=T)
# CV_dark_PAR=(sd_dark_PAR/ave_dark_PAR)*100

#####################################################################################
#### 5. STEP 2, STEP 3, STEP 4 and FLAG assignment for all radiometric channel
######################################################################################

#### 5.A _ CHANNEL: 380 nm ####
  # y is Irradiance at 380 nm (in log)
  # x is PRES (linear)
  # STEP 2a: MAJOR CLOUDS IDENTIFICATION
if(OK380){
if (sum(!is.na(TAB_complete$IRR_380[1:lim380_bis]))>5 ) {
#print("380 bis")
#print(sum(!is.na(TAB_complete$IRR_380[1:lim380_bis])))
pol_380_1=lm(log(TAB_complete$IRR_380[1:lim380_bis]) ~ TAB_complete$PRES[1:lim380_bis] + I(TAB_complete$PRES[1:lim380_bis]^2) + I(TAB_complete$PRES[1:lim380_bis]^3) + I(TAB_complete$PRES[1:lim380_bis]^4)) #not orthogonal fit
summary(pol_380_1)
mean_380=mean(pol_380_1$residuals, na.rm=T)
sd_380=sd(pol_380_1$residuals, na.rm=T)
lim_sd2_380=2*sd_380 # 2 standard deviation #to detect FLAG 3 in STEP 2
lim_sd3_380=3*sd_380 # 3 standard deviation #to detect FLAG 4 in STEP 2
flag3_380=pol_380_1$residuals<(mean_380-lim_sd2_380) | pol_380_1$residuals>(mean_380+lim_sd2_380) # FLAG 3 points in STEP 2
pts.flag3_380=which(flag3_380==T) # FLAG 3 points in STEP 2
flag4_380=pol_380_1$residuals<(mean_380-lim_sd3_380) | pol_380_1$residuals>(mean_380+lim_sd3_380) # FLAG 4 points in STEP 2
pts.flag4_380=which(flag4_380==T) # FLAG 4 points in STEP 2
rsquared_380=summary(pol_380_1)$r.squared
# intercept_380=pol_380_1$coefficients[1]
# x1_380=pol_380_1$coefficients[2]
# x2_380=pol_380_1$coefficients[3]
# x3_380=pol_380_1$coefficients[4]
# x4_380=pol_380_1$coefficients[5]
#END of STEP 2a #REMEMBER! IN THE CODE THE CORRESPONDING STEP 2a FLAGS ARE ASSIGNED TOGETHER THE FLAGS ASSIGNED DURING FOLLOWING STEPS  

no.cloudy_380=which(flag3_380==F) #good records that will be used in STEP 3 
i_depth=which(TAB_complete$PRES>15) #for identifying swinging profiles in STEP 2
flag3plus_380= pol_380_1$residuals[i_depth]>(mean_380+lim_sd2_380)  #for checking sunny points during ovecast profiles below the surface 
pts.flag3plus_380=which(flag3plus_380==T) #for checking sunny points during ovecast profiles below the surface 
flag3minus_380=pol_380_1$residuals[i_depth]<(mean_380-lim_sd2_380) #for checking sunny points during ovecast profiles  below the surface
pts.flag3minus_380=which(flag3minus_380==T) #for checking sunny points during ovecast profiles below the surface

{
  if(length(pts.flag3plus_380)>length(pts.flag3minus_380) & rsquared_380<0.995) { #THIS is STEP 2b:SWINGING PROFILE IDENTIFICATION
    FLAG_380_QC=rep("3", length(TAB_complete$IRR_380)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_380_QC=FLAG_380_QC
    TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
    type380="3" #catherine    
    TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
    #str(newdata_380) 
    # newdata_380$FLAG_380_QC is the right string of characters to be included in NETCDF file
  } 
  #END of STEP 2b
  
  else { #BEGIN OF STEP 3: TYPE PROFILE IDENTIFICATION
    #print("380 no cloudy")
    pol_380_2=lm(log(TAB_complete$IRR_380[no.cloudy_380]) ~ TAB_complete$PRES[no.cloudy_380] + I(TAB_complete$PRES[no.cloudy_380]^2) + I(TAB_complete$PRES[no.cloudy_380]^3) + I(TAB_complete$PRES[no.cloudy_380]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
    #summary(pol_380_2)
    rsquared_380_2=summary(pol_380_2)$r.squared
#     intercept_380_2=pol_380_2$coefficients[1] #intercept
#     x1_380_2=pol_380_2$coefficients[2] #coefficient
#     x2_380_2=pol_380_2$coefficients[3] #coefficient
#     x3_380_2=pol_380_2$coefficients[4] #coefficient
#     x4_380_2=pol_380_2$coefficients[5] #coefficient
         
{
  if(rsquared_380_2<=0.990) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
    FLAG_380_QC=rep("3", length(TAB_complete$IRR_380)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4	
    TAB_complete$FLAG_380_QC=FLAG_380_QC  
    TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
    type380="3" #catherine 
    TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
    #str(newdata_380) 
    # newdata_380$FLAG_380_QC is the right string of characters to be included in NETCDF file
  }
  
  else { 
    if(rsquared_380_2>0.990 & rsquared_380_2<=0.997) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
      FLAG_380_QC=rep("3", length(TAB_complete$IRR_380)) #stop analysis #dark values are still assigned to FLAG 3 from STEP 1
      FLAG_380_QC[pts.flag3_380]="3" #FLAG 3 pts assigned in STEP 2a #no-useful command but it is kept just to remind the STEP of assignement
      FLAG_380_QC[pts.flag4_380]="3" #FLAG 4 pts assigned in STEP 2a
      TAB_complete$FLAG_380_QC=FLAG_380_QC   
      TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
      type380="3" #catherine 
      TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
      newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
      #str(newdata_380)  
      # newdata_380$FLAG_380_QC is the right string of characters to be included in NETCDF file
    }
    
    else{
      if(rsquared_380_2>0.997 & rsquared_380_2<=0.999) { #begin of STEP 4a: TYPE2_FLAG ASSIGNEMENT
        pol_380_3=lm(log(TAB_complete$IRR_380[no.cloudy_380]) ~ TAB_complete$PRES[no.cloudy_380] + I(TAB_complete$PRES[no.cloudy_380]^2) + I(TAB_complete$PRES[no.cloudy_380]^3) + I(TAB_complete$PRES[no.cloudy_380]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_380_3)
        mean_380_3=mean(pol_380_3$residuals, na.rm=T)
        sd_380_3=sd(pol_380_3$residuals, na.rm=T)
        lim_sd2_380_3=2*sd_380_3 # 2 standard deviation #to detect FLAG 3 in STEP 4a
        lim_sd3_380_3=3*sd_380_3 # 3 standard deviation #to detect FLAG 4 in STEP 4a
        flag3_380_ADJ=pol_380_3$residuals<(mean_380_3-lim_sd2_380_3) | pol_380_3$residuals>(mean_380_3+lim_sd2_380_3) #FLAG 3 points in STEP 4a
        pts.flag3_380_ADJ=which(flag3_380_ADJ==T) #FLAG 3 points in STEP 4a
        flag4_380_ADJ=pol_380_3$residuals<(mean_380_3-lim_sd3_380_3) | pol_380_3$residuals>(mean_380_3+lim_sd3_380_3) #FLAG 4 points in STEP 4a
        pts.flag4_380_ADJ=which(flag4_380_ADJ==T) #FLAG 4 points in STEP 4a
        rsquared_380_3=summary(pol_380_3)$r.squared
       #intercept_380_3=pol_380_3$coefficients[1] 
       #x1_380_3=pol_380_3$coefficients[2]
       #x2_380_3=pol_380_3$coefficients[3]
       #x3_380_3=pol_380_3$coefficients[4]
       #x4_380_3=pol_380_3$coefficients[5]
        FLAG_380_QC=rep("2", length(TAB_complete$IRR_380)) # It is a profile without FLAGS 1 # assignment corresponding to STEP 3
        FLAG_380_QC[(lim380_bis+1):length(FLAG_380_QC)]="3" #dark values assigned in STEP 1
        FLAG_380_QC[pts.flag3_380]="3" #FLAG 3 assigned in STEP 2a
        FLAG_380_QC[pts.flag4_380]="3" #FLAG 4 assigned in STEP 2a
        FLAG_380_QC[no.cloudy_380[pts.flag3_380_ADJ]]="3" #FLAG 3 assigned in STEP 4a
        FLAG_380_QC[no.cloudy_380[pts.flag4_380_ADJ]]="3" #FLAG 4 assigned in STEP 4a
        TAB_complete$FLAG_380_QC=FLAG_380_QC
        TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
	type380="2" #catherine 
        TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
        #str(newdata_380)
# newdata_380$FLAG_380_QC is the right string of characters to be included in NETCDF file
      }
      
      else{ #begin of STEP 4b : TYPE1_FLAG ASSIGNEMENT
           pol_380_3=lm(log(TAB_complete$IRR_380[no.cloudy_380]) ~ TAB_complete$PRES[no.cloudy_380] + I(TAB_complete$PRES[no.cloudy_380]^2) + I(TAB_complete$PRES[no.cloudy_380]^3) + I(TAB_complete$PRES[no.cloudy_380]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
           #summary(pol_380_3)
           mean_380_3=mean(pol_380_3$residuals, na.rm=T)
           sd_380_3=sd(pol_380_3$residuals, na.rm=T)
           lim_sd1_380_3=1*sd_380_3 # 1 standard deviation #to detect FLAG 2
           lim_sd2_380_3=2*sd_380_3 # 2 standard deviation #to detect FLAG 3
           lim_sd3_380_3=3*sd_380_3 # 3 standard deviation #to detect FLAG 4
           flag2_380_ADJ=pol_380_3$residuals<(mean_380_3-lim_sd1_380_3) | pol_380_3$residuals>(mean_380_3+lim_sd1_380_3) #FLAG 2 points in STEP 4b
           pts.flag2_380_ADJ=which(flag2_380_ADJ==T) #FLAG 2 points in STEP 4b
           flag3_380_ADJ=pol_380_3$residuals<(mean_380_3-lim_sd2_380_3) | pol_380_3$residuals>(mean_380_3+lim_sd2_380_3) #FLAG 3 points in STEP 4b
           pts.flag3_380_ADJ=which(flag3_380_ADJ==T) #FLAG 3 points in STEP 4b
           flag4_380_ADJ=pol_380_3$residuals<(mean_380_3-lim_sd3_380_3) | pol_380_3$residuals>(mean_380_3+lim_sd3_380_3) #FLAG 4 points in STEP 4b
           pts.flag4_380_ADJ=which(flag4_380_ADJ==T) #FLAG 4 points in STEP 4b
           rsquared_380_3=summary(pol_380_3)$r.squared
           #intercept_380_3=pol_380_3$coefficients[1]
           #x1_380_3=pol_380_3$coefficients[2]
           #x2_380_3=pol_380_3$coefficients[3]
           #x3_380_3=pol_380_3$coefficients[4]
           #x4_380_3=pol_380_3$coefficients[5]
           FLAG_380_QC=rep("1", length(TAB_complete$IRR_380)) #ASSIGNEMENT CORRESPONDING to STEP 3
           FLAG_380_QC[(lim380_bis+1):length(FLAG_380_QC)]="3" #dark flags assigned in STEP 1
           FLAG_380_QC[pts.flag3_380]="3" #FLAG 3 assigned in STEP 2a
           FLAG_380_QC[pts.flag4_380]="3" #FLAG 4 assigned in STEP 2a
           FLAG_380_QC[no.cloudy_380[pts.flag2_380_ADJ]]="2" #FLAG 2 assigned in STEP 4b
           FLAG_380_QC[no.cloudy_380[pts.flag3_380_ADJ]]="3" #FLAG 3 assigned in STEP 4b
           FLAG_380_QC[no.cloudy_380[pts.flag4_380_ADJ]]="3" #FLAG 4 assigned in STEP 4b
           TAB_complete$FLAG_380_QC=FLAG_380_QC   
           TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
	   type380="1" #catherine 
           TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
           newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
           #str(newdata_380)
           # newdata_380$FLAG_380_QC is the right string of characters to be included in NETCDF file
      }}}}}
} 
} else {
	FLAG_380_QC=rep("3", length(TAB_complete$IRR_380)) #ASSIGNEMENT CORRESPONDING to STEP 3
	TAB_complete$FLAG_380_QC=FLAG_380_QC   
	TAB_NA$FLAG_380_QC=rep("NA", length(TAB_NA$PRES))
	type380="3" #catherine 
        TAB_380=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_380 <- TAB_380[order(as.numeric(row.names(TAB_380))),]  #to sort for the initial row order, NA included
        #str(newdata_380)
}
}
# END of RT-QC for the CHANNEL 380 nm
  
  
#### 5.B _ CHANNEL: 412 nm ####
# y is Irradiance at 412 nm (in log)
# x is PRES (linear)
# STEP 2a: MAJOR CLOUDS IDENTIFICATION
if(OK412){
if (sum(!is.na(TAB_complete$IRR_412[1:lim412_bis]))>5 )
{
pol_412_1=lm(log(TAB_complete$IRR_412[1:lim412_bis]) ~ TAB_complete$PRES[1:lim412_bis] + I(TAB_complete$PRES[1:lim412_bis]^2) + I(TAB_complete$PRES[1:lim412_bis]^3) + I(TAB_complete$PRES[1:lim412_bis]^4)) #not orthogonal fit
summary(pol_412_1)
mean_412=mean(pol_412_1$residuals, na.rm=T)
sd_412=sd(pol_412_1$residuals, na.rm=T)
lim_sd2_412=2*sd_412 # 2 standard deviation #to detect FLAG 3 in STEP 2
lim_sd3_412=3*sd_412 # 3 standard deviation #to detect FLAG 4 in STEP 2
flag3_412=pol_412_1$residuals<(mean_412-lim_sd2_412) | pol_412_1$residuals>(mean_412+lim_sd2_412) # FLAG 3 points in STEP 2
pts.flag3_412=which(flag3_412==T) # FLAG 3 points in STEP 2
flag4_412=pol_412_1$residuals<(mean_412-lim_sd3_412) | pol_412_1$residuals>(mean_412+lim_sd3_412) # FLAG 4 points in STEP 2
pts.flag4_412=which(flag4_412==T) # FLAG 4 points in STEP 2
rsquared_412=summary(pol_412_1)$r.squared
# intercept_412=pol_412_1$coefficients[1]
# x1_412=pol_412_1$coefficients[2]
# x2_412=pol_412_1$coefficients[3]
# x3_412=pol_412_1$coefficients[4]
# x4_412=pol_412_1$coefficients[5]
#END of STEP 2a #REMEMBER! IN THE CODE THE CORRESPONDING STEP 2a FLAGS ARE ASSIGNED TOGETHER THE FLAGS ASSIGNED DURING FOLLOWING STEPS  

no.cloudy_412=which(flag3_412==F) #good records that will be used in STEP 3 
i_depth=which(TAB_complete$PRES>15) #for identifying swinging profiles in STEP 2
flag3plus_412= pol_412_1$residuals[i_depth]>(mean_412+lim_sd2_412)  #for checking sunny points during ovecast profiles below the surface 
pts.flag3plus_412=which(flag3plus_412==T) #for checking sunny points during ovecast profiles below the surface 
flag3minus_412=pol_412_1$residuals[i_depth]<(mean_412-lim_sd2_412) #for checking sunny points during ovecast profiles  below the surface
pts.flag3minus_412=which(flag3minus_412==T) #for checking sunny points during ovecast profiles below the surface

{
  if(length(pts.flag3plus_412)>length(pts.flag3minus_412) & rsquared_412<0.995) { #THIS is STEP 2b:SWINGING PROFILE IDENTIFICATION
    FLAG_412_QC=rep("3", length(TAB_complete$IRR_412)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_412_QC=FLAG_412_QC   
    TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
    type412="3" #catherine  
    TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
    #str(newdata_412) 
    # newdata_412$FLAG_412_QC is the right string of characters to be included in NETCDF file
  } 
  #END of STEP 2b
  
  else { #BEGIN OF STEP 3: TYPE PROFILE IDENTIFICATION
    pol_412_2=lm(log(TAB_complete$IRR_412[no.cloudy_412]) ~ TAB_complete$PRES[no.cloudy_412] + I(TAB_complete$PRES[no.cloudy_412]^2) + I(TAB_complete$PRES[no.cloudy_412]^3) + I(TAB_complete$PRES[no.cloudy_412]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
    summary(pol_412_2)
    rsquared_412_2=summary(pol_412_2)$r.squared
    #     intercept_412_2=pol_412_2$coefficients[1] #intercept
    #     x1_412_2=pol_412_2$coefficients[2] #coefficient
    #     x2_412_2=pol_412_2$coefficients[3] #coefficient
    #     x3_412_2=pol_412_2$coefficients[4] #coefficient
    #     x4_412_2=pol_412_2$coefficients[5] #coefficient
    
{
  if(rsquared_412_2<=0.990) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
    FLAG_412_QC=rep("3", length(TAB_complete$IRR_412)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_412_QC=FLAG_412_QC   
    TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
    type412="3" #catherine  
    TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
    #str(newdata_412) 
    # newdata_412$FLAG_412_QC is the right string of characters to be included in NETCDF file
  }
  
  else { 
    if(rsquared_412_2>0.990 & rsquared_412_2<=0.997) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
      FLAG_412_QC=rep("3", length(TAB_complete$IRR_412)) #stop analysis #dark values are still assigned to FLAG 3 from STEP 1
      FLAG_412_QC[pts.flag3_412]="3" #FLAG 3 pts assigned in STEP 2a #no-useful command but it is kept just to remind the STEP of assignement
      FLAG_412_QC[pts.flag4_412]="3" #FLAG 4 pts assigned in STEP 2a
      TAB_complete$FLAG_412_QC=FLAG_412_QC   
      TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
      type412="3" #catherine  
      TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
      newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
      #str(newdata_412)  
      # newdata_412$FLAG_412_QC is the right string of characters to be included in NETCDF file
    }
    
    else{
      if(rsquared_412_2>0.997 & rsquared_412_2<=0.998) { #begin of STEP 4a: TYPE2_FLAG ASSIGNEMENT
        pol_412_3=lm(log(TAB_complete$IRR_412[no.cloudy_412]) ~ TAB_complete$PRES[no.cloudy_412] + I(TAB_complete$PRES[no.cloudy_412]^2) + I(TAB_complete$PRES[no.cloudy_412]^3) + I(TAB_complete$PRES[no.cloudy_412]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_412_3)
        mean_412_3=mean(pol_412_3$residuals, na.rm=T)
        sd_412_3=sd(pol_412_3$residuals, na.rm=T)
        lim_sd2_412_3=2*sd_412_3 # 2 standard deviation #to detect FLAG 3 in STEP 4a
        lim_sd3_412_3=3*sd_412_3 # 3 standard deviation #to detect FLAG 4 in STEP 4a
        flag3_412_ADJ=pol_412_3$residuals<(mean_412_3-lim_sd2_412_3) | pol_412_3$residuals>(mean_412_3+lim_sd2_412_3) #FLAG 3 points in STEP 4a
        pts.flag3_412_ADJ=which(flag3_412_ADJ==T) #FLAG 3 points in STEP 4a
        flag4_412_ADJ=pol_412_3$residuals<(mean_412_3-lim_sd3_412_3) | pol_412_3$residuals>(mean_412_3+lim_sd3_412_3) #FLAG 4 points in STEP 4a
        pts.flag4_412_ADJ=which(flag4_412_ADJ==T) #FLAG 4 points in STEP 4a
        rsquared_412_3=summary(pol_412_3)$r.squared
        #intercept_412_3=pol_412_3$coefficients[1] 
        #x1_412_3=pol_412_3$coefficients[2]
        #x2_412_3=pol_412_3$coefficients[3]
        #x3_412_3=pol_412_3$coefficients[4]
        #x4_412_3=pol_412_3$coefficients[5]
        FLAG_412_QC=rep("2", length(TAB_complete$IRR_412)) # It is a profile without FLAGS 1 # assignment corresponding to STEP 3
        FLAG_412_QC[(lim412_bis+1):length(FLAG_412_QC)]="3" #dark values assigned in STEP 1
        FLAG_412_QC[pts.flag3_412]="3" #FLAG 3 assigned in STEP 2a
        FLAG_412_QC[pts.flag4_412]="3" #FLAG 4 assigned in STEP 2a
        FLAG_412_QC[no.cloudy_412[pts.flag3_412_ADJ]]="3" #FLAG 3 assigned in STEP 4a
        FLAG_412_QC[no.cloudy_412[pts.flag4_412_ADJ]]="3" #FLAG 4 assigned in STEP 4a
        TAB_complete$FLAG_412_QC=FLAG_412_QC
        TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
	type412="2" #catherine  
        TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
        #str(newdata_412)
        # newdata_412$FLAG_412_QC is the right string of characters to be included in NETCDF file
      }
      
      else{ #begin of STEP 4b : TYPE1_FLAG ASSIGNEMENT
        pol_412_3=lm(log(TAB_complete$IRR_412[no.cloudy_412]) ~ TAB_complete$PRES[no.cloudy_412] + I(TAB_complete$PRES[no.cloudy_412]^2) + I(TAB_complete$PRES[no.cloudy_412]^3) + I(TAB_complete$PRES[no.cloudy_412]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_412_3)
        mean_412_3=mean(pol_412_3$residuals, na.rm=T)
        sd_412_3=sd(pol_412_3$residuals, na.rm=T)
        lim_sd1_412_3=1*sd_412_3 # 1 standard deviation #to detect FLAG 2
        lim_sd2_412_3=2*sd_412_3 # 2 standard deviation #to detect FLAG 3
        lim_sd3_412_3=3*sd_412_3 # 3 standard deviation #to detect FLAG 4
        flag2_412_ADJ=pol_412_3$residuals<(mean_412_3-lim_sd1_412_3) | pol_412_3$residuals>(mean_412_3+lim_sd1_412_3) #FLAG 2 points in STEP 4b
        pts.flag2_412_ADJ=which(flag2_412_ADJ==T) #FLAG 2 points in STEP 4b
        flag3_412_ADJ=pol_412_3$residuals<(mean_412_3-lim_sd2_412_3) | pol_412_3$residuals>(mean_412_3+lim_sd2_412_3) #FLAG 3 points in STEP 4b
        pts.flag3_412_ADJ=which(flag3_412_ADJ==T) #FLAG 3 points in STEP 4b
        flag4_412_ADJ=pol_412_3$residuals<(mean_412_3-lim_sd3_412_3) | pol_412_3$residuals>(mean_412_3+lim_sd3_412_3) #FLAG 4 points in STEP 4b
        pts.flag4_412_ADJ=which(flag4_412_ADJ==T) #FLAG 4 points in STEP 4b
        rsquared_412_3=summary(pol_412_3)$r.squared
        #intercept_412_3=pol_412_3$coefficients[1]
        #x1_412_3=pol_412_3$coefficients[2]
        #x2_412_3=pol_412_3$coefficients[3]
        #x3_412_3=pol_412_3$coefficients[4]
        #x4_412_3=pol_412_3$coefficients[5]
        FLAG_412_QC=rep("1", length(TAB_complete$IRR_412)) #ASSIGNEMENT CORRESPONDING to STEP 3
        FLAG_412_QC[(lim412_bis+1):length(FLAG_412_QC)]="3" #dark flags assigned in STEP 1
        FLAG_412_QC[pts.flag3_412]="3" #FLAG 3 assigned in STEP 2a
        FLAG_412_QC[pts.flag4_412]="3" #FLAG 4 assigned in STEP 2a
        FLAG_412_QC[no.cloudy_412[pts.flag2_412_ADJ]]="2" #FLAG 2 assigned in STEP 4b
        FLAG_412_QC[no.cloudy_412[pts.flag3_412_ADJ]]="3" #FLAG 3 assigned in STEP 4b
        FLAG_412_QC[no.cloudy_412[pts.flag4_412_ADJ]]="3" #FLAG 4 assigned in STEP 4b
        TAB_complete$FLAG_412_QC=FLAG_412_QC   
        TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
	type412="1" #catherine  
        TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
        #str(newdata_412)
        # newdata_412$FLAG_412_QC is the right string of characters to be included in NETCDF file
      }}}}}
} 
} else {
	FLAG_412_QC=rep("3", length(TAB_complete$IRR_412)) #ASSIGNEMENT CORRESPONDING to STEP 3
	TAB_complete$FLAG_412_QC=FLAG_412_QC   
	TAB_NA$FLAG_412_QC=rep("NA", length(TAB_NA$PRES))
	type412="3" #catherine 
        TAB_412=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_412 <- TAB_412[order(as.numeric(row.names(TAB_412))),]  #to sort for the initial row order, NA included
        #str(newdata_412)
}
}
# END of RT-QC for the CHANNEL 412 nm  
  
  
#### 5.C _ CHANNEL: 490 nm ####
# y is Irradiance at 490 nm (in log)
# x is PRES (linear)
# STEP 2a: MAJOR CLOUDS IDENTIFICATION
if(OK490){
if (sum(!is.na(TAB_complete$IRR_490[1:lim490_bis]))>5){
pol_490_1=lm(log(TAB_complete$IRR_490[1:lim490_bis]) ~ TAB_complete$PRES[1:lim490_bis] + I(TAB_complete$PRES[1:lim490_bis]^2) + I(TAB_complete$PRES[1:lim490_bis]^3) + I(TAB_complete$PRES[1:lim490_bis]^4)) #not orthogonal fit
summary(pol_490_1)
mean_490=mean(pol_490_1$residuals, na.rm=T)
sd_490=sd(pol_490_1$residuals, na.rm=T)
lim_sd2_490=2*sd_490 # 2 standard deviation #to detect FLAG 3 in STEP 2
lim_sd3_490=3*sd_490 # 3 standard deviation #to detect FLAG 4 in STEP 2
flag3_490=pol_490_1$residuals<(mean_490-lim_sd2_490) | pol_490_1$residuals>(mean_490+lim_sd2_490) # FLAG 3 points in STEP 2
pts.flag3_490=which(flag3_490==T) # FLAG 3 points in STEP 2
flag4_490=pol_490_1$residuals<(mean_490-lim_sd3_490) | pol_490_1$residuals>(mean_490+lim_sd3_490) # FLAG 4 points in STEP 2
pts.flag4_490=which(flag4_490==T) # FLAG 4 points in STEP 2
rsquared_490=summary(pol_490_1)$r.squared
# intercept_490=pol_490_1$coefficients[1]
# x1_490=pol_490_1$coefficients[2]
# x2_490=pol_490_1$coefficients[3]
# x3_490=pol_490_1$coefficients[4]
# x4_490=pol_490_1$coefficients[5]
#END of STEP 2a #REMEMBER! IN THE CODE THE CORRESPONDING STEP 2a FLAGS ARE ASSIGNED TOGETHER THE FLAGS ASSIGNED DURING FOLLOWING STEPS  

no.cloudy_490=which(flag3_490==F) #good records that will be used in STEP 3 
i_depth=which(TAB_complete$PRES>15) #for identifying swinging profiles in STEP 2
flag3plus_490= pol_490_1$residuals[i_depth]>(mean_490+lim_sd2_490)  #for checking sunny points during ovecast profiles below the surface 
pts.flag3plus_490=which(flag3plus_490==T) #for checking sunny points during ovecast profiles below the surface 
flag3minus_490=pol_490_1$residuals[i_depth]<(mean_490-lim_sd2_490) #for checking sunny points during ovecast profiles  below the surface
pts.flag3minus_490=which(flag3minus_490==T) #for checking sunny points during ovecast profiles below the surface

{
  if(length(pts.flag3plus_490)>length(pts.flag3minus_490) & rsquared_490<0.995) { #THIS is STEP 2b:SWINGING PROFILE IDENTIFICATION
    FLAG_490_QC=rep("3", length(TAB_complete$IRR_490)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_490_QC=FLAG_490_QC   
    TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
    type490="3" #catherine  
    TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
    #str(newdata_490) 
    # newdata_490$FLAG_490_QC is the right string of characters to be included in NETCDF file
  } 
  #END of STEP 2b
  
  else { #BEGIN OF STEP 3: TYPE PROFILE IDENTIFICATION
    pol_490_2=lm(log(TAB_complete$IRR_490[no.cloudy_490]) ~ TAB_complete$PRES[no.cloudy_490] + I(TAB_complete$PRES[no.cloudy_490]^2) + I(TAB_complete$PRES[no.cloudy_490]^3) + I(TAB_complete$PRES[no.cloudy_490]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
    summary(pol_490_2)
    rsquared_490_2=summary(pol_490_2)$r.squared
    #     intercept_490_2=pol_490_2$coefficients[1] #intercept
    #     x1_490_2=pol_490_2$coefficients[2] #coefficient
    #     x2_490_2=pol_490_2$coefficients[3] #coefficient
    #     x3_490_2=pol_490_2$coefficients[4] #coefficient
    #     x4_490_2=pol_490_2$coefficients[5] #coefficient
    
{
  if(rsquared_490_2<=0.990) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
    FLAG_490_QC=rep("3", length(TAB_complete$IRR_490)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_490_QC=FLAG_490_QC   
    TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
    type490="3" #catherine  
    TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
    #str(newdata_490) 
    # newdata_490$FLAG_490_QC is the right string of characters to be included in NETCDF file
  }
  
  else { 
    if(rsquared_490_2>0.990 & rsquared_490_2<=0.996) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
      FLAG_490_QC=rep("3", length(TAB_complete$IRR_490)) #stop analysis #dark values are still assigned to FLAG 3 from STEP 1
      FLAG_490_QC[pts.flag3_490]="3" #FLAG 3 pts assigned in STEP 2a #no-useful command but it is kept just to remind the STEP of assignement
      FLAG_490_QC[pts.flag4_490]="3" #FLAG 4 pts assigned in STEP 2a
      TAB_complete$FLAG_490_QC=FLAG_490_QC   
      TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
      type490="3" #catherine  
      TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
      newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
      #str(newdata_490)  
      # newdata_490$FLAG_490_QC is the right string of characters to be included in NETCDF file
    }
    
    else{
      if(rsquared_490_2>0.996 & rsquared_490_2<=0.998) { #begin of STEP 4a: TYPE2_FLAG ASSIGNEMENT
        pol_490_3=lm(log(TAB_complete$IRR_490[no.cloudy_490]) ~ TAB_complete$PRES[no.cloudy_490] + I(TAB_complete$PRES[no.cloudy_490]^2) + I(TAB_complete$PRES[no.cloudy_490]^3) + I(TAB_complete$PRES[no.cloudy_490]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_490_3)
        mean_490_3=mean(pol_490_3$residuals, na.rm=T)
        sd_490_3=sd(pol_490_3$residuals, na.rm=T)
        lim_sd2_490_3=2*sd_490_3 # 2 standard deviation #to detect FLAG 3 in STEP 4a
        lim_sd3_490_3=3*sd_490_3 # 3 standard deviation #to detect FLAG 4 in STEP 4a
        flag3_490_ADJ=pol_490_3$residuals<(mean_490_3-lim_sd2_490_3) | pol_490_3$residuals>(mean_490_3+lim_sd2_490_3) #FLAG 3 points in STEP 4a
        pts.flag3_490_ADJ=which(flag3_490_ADJ==T) #FLAG 3 points in STEP 4a
        flag4_490_ADJ=pol_490_3$residuals<(mean_490_3-lim_sd3_490_3) | pol_490_3$residuals>(mean_490_3+lim_sd3_490_3) #FLAG 4 points in STEP 4a
        pts.flag4_490_ADJ=which(flag4_490_ADJ==T) #FLAG 4 points in STEP 4a
        rsquared_490_3=summary(pol_490_3)$r.squared
        #intercept_490_3=pol_490_3$coefficients[1] 
        #x1_490_3=pol_490_3$coefficients[2]
        #x2_490_3=pol_490_3$coefficients[3]
        #x3_490_3=pol_490_3$coefficients[4]
        #x4_490_3=pol_490_3$coefficients[5]
        FLAG_490_QC=rep("2", length(TAB_complete$IRR_490)) # It is a profile without FLAGS 1 # assignment corresponding to STEP 3
        FLAG_490_QC[(lim490_bis+1):length(FLAG_490_QC)]="3" #dark values assigned in STEP 1
        FLAG_490_QC[pts.flag3_490]="3" #FLAG 3 assigned in STEP 2a
        FLAG_490_QC[pts.flag4_490]="3" #FLAG 4 assigned in STEP 2a
        FLAG_490_QC[no.cloudy_490[pts.flag3_490_ADJ]]="3" #FLAG 3 assigned in STEP 4a
        FLAG_490_QC[no.cloudy_490[pts.flag4_490_ADJ]]="3" #FLAG 4 assigned in STEP 4a
        TAB_complete$FLAG_490_QC=FLAG_490_QC
        TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
        type490="2" #catherine  
        TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
        #str(newdata_490)
        # newdata_490$FLAG_490_QC is the right string of characters to be included in NETCDF file
      }
      
      else{ #begin of STEP 4b : TYPE1_FLAG ASSIGNEMENT
        pol_490_3=lm(log(TAB_complete$IRR_490[no.cloudy_490]) ~ TAB_complete$PRES[no.cloudy_490] + I(TAB_complete$PRES[no.cloudy_490]^2) + I(TAB_complete$PRES[no.cloudy_490]^3) + I(TAB_complete$PRES[no.cloudy_490]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_490_3)
        mean_490_3=mean(pol_490_3$residuals, na.rm=T)
        sd_490_3=sd(pol_490_3$residuals, na.rm=T)
        lim_sd1_490_3=1*sd_490_3 # 1 standard deviation #to detect FLAG 2
        lim_sd2_490_3=2*sd_490_3 # 2 standard deviation #to detect FLAG 3
        lim_sd3_490_3=3*sd_490_3 # 3 standard deviation #to detect FLAG 4
        flag2_490_ADJ=pol_490_3$residuals<(mean_490_3-lim_sd1_490_3) | pol_490_3$residuals>(mean_490_3+lim_sd1_490_3) #FLAG 2 points in STEP 4b
        pts.flag2_490_ADJ=which(flag2_490_ADJ==T) #FLAG 2 points in STEP 4b
        flag3_490_ADJ=pol_490_3$residuals<(mean_490_3-lim_sd2_490_3) | pol_490_3$residuals>(mean_490_3+lim_sd2_490_3) #FLAG 3 points in STEP 4b
        pts.flag3_490_ADJ=which(flag3_490_ADJ==T) #FLAG 3 points in STEP 4b
        flag4_490_ADJ=pol_490_3$residuals<(mean_490_3-lim_sd3_490_3) | pol_490_3$residuals>(mean_490_3+lim_sd3_490_3) #FLAG 4 points in STEP 4b
        pts.flag4_490_ADJ=which(flag4_490_ADJ==T) #FLAG 4 points in STEP 4b
        rsquared_490_3=summary(pol_490_3)$r.squared
        #intercept_490_3=pol_490_3$coefficients[1]
        #x1_490_3=pol_490_3$coefficients[2]
        #x2_490_3=pol_490_3$coefficients[3]
        #x3_490_3=pol_490_3$coefficients[4]
        #x4_490_3=pol_490_3$coefficients[5]
        FLAG_490_QC=rep("1", length(TAB_complete$IRR_490)) #ASSIGNEMENT CORRESPONDING to STEP 3
        FLAG_490_QC[(lim490_bis+1):length(FLAG_490_QC)]="3" #dark flags assigned in STEP 1
        FLAG_490_QC[pts.flag3_490]="3" #FLAG 3 assigned in STEP 2a
        FLAG_490_QC[pts.flag4_490]="3" #FLAG 4 assigned in STEP 2a
        FLAG_490_QC[no.cloudy_490[pts.flag2_490_ADJ]]="2" #FLAG 2 assigned in STEP 4b
        FLAG_490_QC[no.cloudy_490[pts.flag3_490_ADJ]]="3" #FLAG 3 assigned in STEP 4b
        FLAG_490_QC[no.cloudy_490[pts.flag4_490_ADJ]]="3" #FLAG 4 assigned in STEP 4b
        TAB_complete$FLAG_490_QC=FLAG_490_QC   
        TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
	type490="1" #catherine  
        TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
        #str(newdata_490)
        # newdata_490$FLAG_490_QC is the right string of characters to be included in NETCDF file
      }}}}}
} 
} else {
	FLAG_490_QC=rep("3", length(TAB_complete$IRR_490)) #ASSIGNEMENT CORRESPONDING to STEP 3
	TAB_complete$FLAG_490_QC=FLAG_490_QC   
	TAB_NA$FLAG_490_QC=rep("NA", length(TAB_NA$PRES))
	type490="3" #catherine 
        TAB_490=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_490 <- TAB_490[order(as.numeric(row.names(TAB_490))),]  #to sort for the initial row order, NA included
        #str(newdata_490)
}
}
# END of RT-QC for the CHANNEL 490 nm    
  

#### 5.D _ CHANNEL: PAR ####
# y is PAR (in log)
# x is PRES (linear)
# STEP 2a: MAJOR CLOUDS IDENTIFICATION
if(OKPAR){
if (sum(!is.na(TAB_complete$PAR[1:limPAR_bis]))>5 ){
pol_PAR_1=lm(log(TAB_complete$PAR[1:limPAR_bis]) ~ TAB_complete$PRES[1:limPAR_bis] + I(TAB_complete$PRES[1:limPAR_bis]^2) + I(TAB_complete$PRES[1:limPAR_bis]^3) + I(TAB_complete$PRES[1:limPAR_bis]^4)) #not orthogonal fit
summary(pol_PAR_1)
mean_PAR=mean(pol_PAR_1$residuals, na.rm=T)
sd_PAR=sd(pol_PAR_1$residuals, na.rm=T)
lim_sd2_PAR=2*sd_PAR # 2 standard deviation #to detect FLAG 3 in STEP 2
lim_sd3_PAR=3*sd_PAR # 3 standard deviation #to detect FLAG 4 in STEP 2
flag3_PAR=pol_PAR_1$residuals<(mean_PAR-lim_sd2_PAR) | pol_PAR_1$residuals>(mean_PAR+lim_sd2_PAR) # FLAG 3 points in STEP 2
pts.flag3_PAR=which(flag3_PAR==T) # FLAG 3 points in STEP 2
flag4_PAR=pol_PAR_1$residuals<(mean_PAR-lim_sd3_PAR) | pol_PAR_1$residuals>(mean_PAR+lim_sd3_PAR) # FLAG 4 points in STEP 2
pts.flag4_PAR=which(flag4_PAR==T) # FLAG 4 points in STEP 2
rsquared_PAR=summary(pol_PAR_1)$r.squared
# intercept_PAR=pol_PAR_1$coefficients[1]
# x1_PAR=pol_PAR_1$coefficients[2]
# x2_PAR=pol_PAR_1$coefficients[3]
# x3_PAR=pol_PAR_1$coefficients[4]
# x4_PAR=pol_PAR_1$coefficients[5]
#END of STEP 2a #REMEMBER! IN THE CODE THE CORRESPONDING STEP 2a FLAGS ARE ASSIGNED TOGETHER THE FLAGS ASSIGNED DURING FOLLOWING STEPS  

no.cloudy_PAR=which(flag3_PAR==F) #good records that will be used in STEP 3 
i_depth=which(TAB_complete$PRES>15) #for identifying swinging profiles in STEP 2
flag3plus_PAR= pol_PAR_1$residuals[i_depth]>(mean_PAR+lim_sd2_PAR)  #for checking sunny points during ovecast profiles below the surface 
pts.flag3plus_PAR=which(flag3plus_PAR==T) #for checking sunny points during ovecast profiles below the surface 
flag3minus_PAR=pol_PAR_1$residuals[i_depth]<(mean_PAR-lim_sd2_PAR) #for checking sunny points during ovecast profiles  below the surface
pts.flag3minus_PAR=which(flag3minus_PAR==T) #for checking sunny points during ovecast profiles below the surface

{
  if(length(pts.flag3plus_PAR)>length(pts.flag3minus_PAR) & rsquared_PAR<0.995) { #THIS is STEP 2b:SWINGING PROFILE IDENTIFICATION
    FLAG_PAR_QC=rep("3", length(TAB_complete$PAR)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC   
    TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
    typePAR="3" #catherine  
    TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
    #str(newdata_PAR) 
    # newdata_PAR$FLAG_PAR_QC is the right string of characters to be included in NETCDF file
  } 
  #END of STEP 2b
  
  else { #BEGIN OF STEP 3: TYPE PROFILE IDENTIFICATION
    pol_PAR_2=lm(log(TAB_complete$PAR[no.cloudy_PAR]) ~ TAB_complete$PRES[no.cloudy_PAR] + I(TAB_complete$PRES[no.cloudy_PAR]^2) + I(TAB_complete$PRES[no.cloudy_PAR]^3) + I(TAB_complete$PRES[no.cloudy_PAR]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
    summary(pol_PAR_2)
    rsquared_PAR_2=summary(pol_PAR_2)$r.squared
    #     intercept_PAR_2=pol_PAR_2$coefficients[1] #intercept
    #     x1_PAR_2=pol_PAR_2$coefficients[2] #coefficient
    #     x2_PAR_2=pol_PAR_2$coefficients[3] #coefficient
    #     x3_PAR_2=pol_PAR_2$coefficients[4] #coefficient
    #     x4_PAR_2=pol_PAR_2$coefficients[5] #coefficient
    
{
  if(rsquared_PAR_2<=0.990) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
    FLAG_PAR_QC=rep("3", length(TAB_complete$PAR)) #stop analysis #also dark values are degraded to FLAG 4 #ALL POINTS ARE ASSIGNED TO FLAG 4
    TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC   
    TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
    typePAR="3" #catherine  
    TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
    newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
    #str(newdata_PAR) 
    # newdata_PAR$FLAG_PAR_QC is the right string of characters to be included in NETCDF file
  }
  
  else { 
    if(rsquared_PAR_2>0.990 & rsquared_PAR_2<=0.996) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
      FLAG_PAR_QC=rep("3", length(TAB_complete$PAR)) #stop analysis #dark values are still assigned to FLAG 3 from STEP 1
      FLAG_PAR_QC[pts.flag3_PAR]="3" #FLAG 3 pts assigned in STEP 2a #no-useful command but it is kept just to remind the STEP of assignement
      FLAG_PAR_QC[pts.flag4_PAR]="3" #FLAG 4 pts assigned in STEP 2a
      TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC   
      TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
      typePAR="3" #catherine  
      TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
      newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
      #str(newdata_PAR)  
      # newdata_PAR$FLAG_PAR_QC is the right string of characters to be included in NETCDF file
    }
    
    else{
      if(rsquared_PAR_2>0.996 & rsquared_PAR_2<=0.998) { #begin of STEP 4a: TYPE2_FLAG ASSIGNEMENT
        pol_PAR_3=lm(log(TAB_complete$PAR[no.cloudy_PAR]) ~ TAB_complete$PRES[no.cloudy_PAR] + I(TAB_complete$PRES[no.cloudy_PAR]^2) + I(TAB_complete$PRES[no.cloudy_PAR]^3) + I(TAB_complete$PRES[no.cloudy_PAR]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_PAR_3)
        mean_PAR_3=mean(pol_PAR_3$residuals, na.rm=T)
        sd_PAR_3=sd(pol_PAR_3$residuals, na.rm=T)
        lim_sd2_PAR_3=2*sd_PAR_3 # 2 standard deviation #to detect FLAG 3 in STEP 4a
        lim_sd3_PAR_3=3*sd_PAR_3 # 3 standard deviation #to detect FLAG 4 in STEP 4a
        flag3_PAR_ADJ=pol_PAR_3$residuals<(mean_PAR_3-lim_sd2_PAR_3) | pol_PAR_3$residuals>(mean_PAR_3+lim_sd2_PAR_3) #FLAG 3 points in STEP 4a
        pts.flag3_PAR_ADJ=which(flag3_PAR_ADJ==T) #FLAG 3 points in STEP 4a
        flag4_PAR_ADJ=pol_PAR_3$residuals<(mean_PAR_3-lim_sd3_PAR_3) | pol_PAR_3$residuals>(mean_PAR_3+lim_sd3_PAR_3) #FLAG 4 points in STEP 4a
        pts.flag4_PAR_ADJ=which(flag4_PAR_ADJ==T) #FLAG 4 points in STEP 4a
        rsquared_PAR_3=summary(pol_PAR_3)$r.squared
        #intercept_PAR_3=pol_PAR_3$coefficients[1] 
        #x1_PAR_3=pol_PAR_3$coefficients[2]
        #x2_PAR_3=pol_PAR_3$coefficients[3]
        #x3_PAR_3=pol_PAR_3$coefficients[4]
        #x4_PAR_3=pol_PAR_3$coefficients[5]
        FLAG_PAR_QC=rep("2", length(TAB_complete$PAR)) # It is a profile without FLAGS 1 # assignment corresponding to STEP 3
        FLAG_PAR_QC[(limPAR_bis+1):length(FLAG_PAR_QC)]="3" #dark values assigned in STEP 1
	FLAG_PAR_QC[pts.flag3_PAR]="3" #FLAG 3 assigned in STEP 2a
        FLAG_PAR_QC[pts.flag4_PAR]="3" #FLAG 4 assigned in STEP 2a
        FLAG_PAR_QC[no.cloudy_PAR[pts.flag3_PAR_ADJ]]="3" #FLAG 3 assigned in STEP 4a
        FLAG_PAR_QC[no.cloudy_PAR[pts.flag4_PAR_ADJ]]="3" #FLAG 4 assigned in STEP 4a
        TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC
        TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
	typePAR="2" #catherine  
        TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
        #str(newdata_PAR)
        # newdata_PAR$FLAG_PAR_QC is the right string of characters to be included in NETCDF file
      }
      
      else{ #begin of STEP 4b : TYPE1_FLAG ASSIGNEMENT
        pol_PAR_3=lm(log(TAB_complete$PAR[no.cloudy_PAR]) ~ TAB_complete$PRES[no.cloudy_PAR] + I(TAB_complete$PRES[no.cloudy_PAR]^2) + I(TAB_complete$PRES[no.cloudy_PAR]^3) + I(TAB_complete$PRES[no.cloudy_PAR]^4)) #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        summary(pol_PAR_3)
        mean_PAR_3=mean(pol_PAR_3$residuals, na.rm=T)
        sd_PAR_3=sd(pol_PAR_3$residuals, na.rm=T)
        lim_sd1_PAR_3=1*sd_PAR_3 # 1 standard deviation #to detect FLAG 2
        lim_sd2_PAR_3=2*sd_PAR_3 # 2 standard deviation #to detect FLAG 3
        lim_sd3_PAR_3=3*sd_PAR_3 # 3 standard deviation #to detect FLAG 4
        flag2_PAR_ADJ=pol_PAR_3$residuals<(mean_PAR_3-lim_sd1_PAR_3) | pol_PAR_3$residuals>(mean_PAR_3+lim_sd1_PAR_3) #FLAG 2 points in STEP 4b
        pts.flag2_PAR_ADJ=which(flag2_PAR_ADJ==T) #FLAG 2 points in STEP 4b
        flag3_PAR_ADJ=pol_PAR_3$residuals<(mean_PAR_3-lim_sd2_PAR_3) | pol_PAR_3$residuals>(mean_PAR_3+lim_sd2_PAR_3) #FLAG 3 points in STEP 4b
        pts.flag3_PAR_ADJ=which(flag3_PAR_ADJ==T) #FLAG 3 points in STEP 4b
        flag4_PAR_ADJ=pol_PAR_3$residuals<(mean_PAR_3-lim_sd3_PAR_3) | pol_PAR_3$residuals>(mean_PAR_3+lim_sd3_PAR_3) #FLAG 4 points in STEP 4b
        pts.flag4_PAR_ADJ=which(flag4_PAR_ADJ==T) #FLAG 4 points in STEP 4b
        rsquared_PAR_3=summary(pol_PAR_3)$r.squared
        #intercept_PAR_3=pol_PAR_3$coefficients[1]
        #x1_PAR_3=pol_PAR_3$coefficients[2]
        #x2_PAR_3=pol_PAR_3$coefficients[3]
        #x3_PAR_3=pol_PAR_3$coefficients[4]
        #x4_PAR_3=pol_PAR_3$coefficients[5]
        FLAG_PAR_QC=rep("1", length(TAB_complete$PAR)) #ASSIGNEMENT CORRESPONDING to STEP 3
        FLAG_PAR_QC[(limPAR_bis+1):length(FLAG_PAR_QC)]="3" #dark flags assigned in STEP 1
        FLAG_PAR_QC[pts.flag3_PAR]="3" #FLAG 3 assigned in STEP 2a
        FLAG_PAR_QC[pts.flag4_PAR]="3" #FLAG 4 assigned in STEP 2a
        FLAG_PAR_QC[no.cloudy_PAR[pts.flag2_PAR_ADJ]]="2" #FLAG 2 assigned in STEP 4b
        FLAG_PAR_QC[no.cloudy_PAR[pts.flag3_PAR_ADJ]]="3" #FLAG 3 assigned in STEP 4b
        FLAG_PAR_QC[no.cloudy_PAR[pts.flag4_PAR_ADJ]]="3" #FLAG 4 assigned in STEP 4b
        TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC   
        TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
	typePAR="1" #catherine  
        TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
        #str(newdata_PAR)
        # newdata_PAR$FLAG_PAR_QC is the right string of characters to be included in NETCDF file
      }}}}}
} 
} else {
	FLAG_PAR_QC=rep("3", length(TAB_complete$PAR)) #ASSIGNEMENT CORRESPONDING to STEP 3
	TAB_complete$FLAG_PAR_QC=FLAG_PAR_QC   
	TAB_NA$FLAG_PAR_QC=rep("NA", length(TAB_NA$PRES))
	typePAR="3" #catherine 
        TAB_PAR=rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_PAR <- TAB_PAR[order(as.numeric(row.names(TAB_PAR))),]  #to sort for the initial row order, NA included
        #str(newdata_PAR)
}
}
# END of RT-QC for the CHANNEL PAR

return(list("FLAG_380"=newdata_380$FLAG_380_QC, "FLAG_412"=newdata_412$FLAG_412_QC, "FLAG_490"=newdata_490$FLAG_490_QC, 
            "FLAG_PAR"=newdata_PAR$FLAG_PAR_QC, "type380"=type380, "type412"=type412, "type490"=type490, "typePAR"=typePAR))

#print(IDnc)
}
  
#### THE END ####  
