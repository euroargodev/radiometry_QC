require(nortest)# Lilliefors test of normality

RT_QC_radiometry <- function(PRES, IRR_380, IRR_412, IRR_490, PAR) {

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

    TAB = cbind(PRES, IRR_380,IRR_412, IRR_490, PAR) #create a matrix using all the parameters 
    TAB = as.data.frame(TAB) #convert TAB in dataframe
  
    ##################################################
    #### 3. Create a data.frame without NA 
    ##################################################
    TAB_complete = TAB[complete.cases(TAB$IRR_380),]#to remove rows with NA. Same rows for any irradiance and PAR measurements.
    TAB_NA = TAB[!complete.cases(TAB$IRR_380),] #dataframe with only NA values. Useful when assignin flags.
  
    ##################################################
    #### 4. STEP 1: DARK IDENTIFICATION 
    ##################################################
  
    PARAM_NAMES = c("IRR_380", "IRR_412", "IRR_490", "PAR")
    FLAG_NAMES = list("IRR_380"="FLAG_380_QC", "IRR_412"="FLAG_412_QC", "IRR_490"="FLAG_490_QC", "PAR"="FLAG_PAR_QC")
    
    swinging_thresh = 0.995
    thresholds = list("IRR_380" = c(0.990, 0.997, 0.999), 
                      "IRR_412" = c(0.990, 0.997, 0.998), 
                      "IRR_490" = c(0.990, 0.996, 0.998), 
                      "PAR" = c(0.990, 0.996, 0.998))
    
    typeAll = NULL
    newdata_All = NULL
    
    fill_flag_3 <- function(param_name) {
        ### method to fill a parameter flag axis with "3" and flag the profile as "3" as well
        
        FLAG_param_QC = rep("3", length(TAB_complete[[param_name]])) #ASSIGNEMENT CORRESPONDING to STEP 3
        TAB_complete[[ FLAG_NAMES[[param_name]] ]] <<- FLAG_param_QC
        TAB_NA[[ FLAG_NAMES[[param_name]] ]] <<- rep("NA", length(TAB_NA$PRES))
        
        TAB_param = rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_All[[ FLAG_NAMES[[param_name]] ]] <<- TAB_param[order(as.numeric(row.names(TAB_380))),][[ FLAG_NAMES[[param_name]] ]]  #to sort for the initial row order, NA included
    }
  
    for (param_name in PARAM_NAMES) {
        
        #### 4.A Lilliefors Test for detecting lowest good irradiance measurement
        # alfa selected at 0.01
        Lilliefors_param = rep(NA, length(TAB_complete$IRR_380)-4)
        for (l in 1 :(length(TAB_complete$IRR_380)-4)) {
            Lilliefors_param[l] = lillie.test(TAB_complete[[param_name]][ l:length(TAB_complete[[param_name]]) ])[2] #select p-value
        }
        i_param = which(abs(unlist(Lilliefors_param)) > 0.01)
        lim = i_param[1] - 1
        
        
        if (is.na(lim)) { ### There should be no way lim is NA ??
            fill_flag_3(param_name)
            typeAll[[param_name]] = "3" #catherine
            next
        }    
        
        #### 4.B Check for negative values
        neg_param = which(TAB_complete[[param_name]][ 1:lim ] <= 0)
        if (length(neg_param) != 0) {
            lim = length(TAB_complete[[param_name]][ 1:neg_param[1] ]) - 1
        }
            
        #####################################################################################
        #### 5. STEP 2, STEP 3, STEP 4 and FLAG assignment for all radiometric channel
        #####################################################################################
            
        if (sum(!is.na(TAB_complete[[param_name]][ 1:lim ])) <= 5 ) {
                fill_flag_3(param_name)
                typeAll[[param_name]] = "3" #catherine
                next
        }
                
        # STEP 2a: MAJOR CLOUDS IDENTIFICATION
        pol_param_1 = lm(log(TAB_complete[[param_name]][ 1:lim ]) ~ 
                             TAB_complete$PRES[ 1:lim ] + I(TAB_complete$PRES[ 1:lim ]^2) + 
                             I(TAB_complete$PRES[ 1:lim ]^3) + I(TAB_complete$PRES[ 1:lim ]^4)) 
                             # not orthogonal fit
 
        mean_param = mean(pol_param_1$residuals, na.rm=T)
        sd_param = sd(pol_param_1$residuals, na.rm=T)
        lim_sd2_param = 2*sd_param # 2 standard deviation #to detect FLAG 3 in STEP 2
        #lim_sd3_param = 3*sd_param # 3 standard deviation #to detect FLAG 4 in STEP 2
                
        flag3_param = pol_param_1$residuals < (mean_param - lim_sd2_param) | 
                      pol_param_1$residuals > (mean_param + lim_sd2_param) # FLAG 3 points in STEP 2
        pts.flag3_param = which(flag3_param) # FLAG 3 points in STEP 2
        #flag4_param = pol_param_1$residuals < (mean_param-lim_sd3_param) | pol_param_1$residuals > (mean_param+lim_sd3_param) # FLAG 4 points in STEP 2
        #pts.flag4_param = which(flag4_param) # FLAG 4 points in STEP 2
            
        rsquared_param = summary(pol_param_1)$r.squared
        #END of STEP 2a #REMEMBER! IN THE CODE THE CORRESPONDING STEP 2a FLAGS ARE ASSIGNED TOGETHER THE FLAGS ASSIGNED DURING FOLLOWING STEPS  
            
        no.cloudy_param = which(!flag3_param) #good records that will be used in STEP 3 
        i_depth = which(TAB_complete$PRES > 15) #for identifying swinging profiles in STEP 2
        flag3plus_param = pol_param_1$residuals[i_depth] > (mean_param + lim_sd2_param)  #for checking sunny points during ovecast profiles below the surface 
        pts.flag3plus_param = which(flag3plus_param) #for checking sunny points during ovecast profiles below the surface 
        flag3minus_param = pol_param_1$residuals[i_depth] < (mean_param - lim_sd2_param) #for checking sunny points during ovecast profiles  below the surface
        pts.flag3minus_param = which(flag3minus_param) #for checking sunny points during ovecast profiles below the surface
            
        # STEP 2b:SWINGING PROFILE IDENTIFICATION
        if ( (length(pts.flag3plus_param) > length(pts.flag3minus_param)) & (rsquared_param < swinging_thresh) ) { 
            fill_flag_3(param_name)
            typeAll[[param_name]] = "3" #catherine
            next
        }
                
        #BEGIN OF STEP 3: TYPE PROFILE IDENTIFICATION
        pol_param_2 = lm(log(TAB_complete[[param_name]][no.cloudy_param]) ~ 
                             TAB_complete$PRES[no.cloudy_param] + I(TAB_complete$PRES[no.cloudy_param]^2) + 
                             I(TAB_complete$PRES[no.cloudy_param]^3) + I(TAB_complete$PRES[no.cloudy_param]^4)) 
                             #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        rsquared_param_2 = summary(pol_param_2)$r.squared
                
        pol_param_3 = lm(log(TAB_complete[[param_name]][no.cloudy_param]) ~ 
                             TAB_complete$PRES[no.cloudy_param] + I(TAB_complete$PRES[no.cloudy_param]^2) + 
                             I(TAB_complete$PRES[no.cloudy_param]^3) + I(TAB_complete$PRES[no.cloudy_param]^4)) 
                             #not orthogonal fit made without using FLAG 3 and 4 pts form STEP 2a
        mean_param_3 = mean(pol_param_3$residuals, na.rm=T)
        sd_param_3 = sd(pol_param_3$residuals, na.rm=T)
        lim_sd1_param_3 = 1*sd_param_3 # 1 standard deviation #to detect FLAG 2
        lim_sd2_param_3 = 2*sd_param_3 # 2 standard deviation #to detect FLAG 3
        #lim_sd3_param_3 = 3*sd_param_3 # 3 standard deviation #to detect FLAG 4
                
        flag2_param_ADJ = pol_param_3$residuals < (mean_param_3 - lim_sd1_param_3) | 
                          pol_param_3$residuals > (mean_param_3 + lim_sd1_param_3) #FLAG 2 points in STEP 4b
        pts.flag2_param_ADJ = which(flag2_param_ADJ) #FLAG 2 points in STEP 4b
        flag3_param_ADJ = pol_param_3$residuals < (mean_param_3 - lim_sd2_param_3) | 
                          pol_param_3$residuals > (mean_param_3 + lim_sd2_param_3) #FLAG 3 points in STEP 4b
        pts.flag3_param_ADJ = which(flag3_param_ADJ) #FLAG 3 points in STEP 4b
        #flag4_param_ADJ = pol_param_3$residuals < (mean_param_3 - lim_sd3_param_3) | 
        #                  pol_param_3$residuals > (mean_param_3 + lim_sd3_param_3) #FLAG 4 points in STEP 4b
        #pts.flag4_param_ADJ=which(flag4_param_ADJ==T) #FLAG 4 points in STEP 4b
        
        rsquared_param_3 = summary(pol_param_3)$r.squared
                    
        if (rsquared_param_2 <= thresholds[[param_name]][1]) { #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
            fill_flag_3(param_name)
            typeAll[[param_name]] = "3" #catherine
            next
        }
        if (rsquared_param_2 > thresholds[[param_name]][1] & rsquared_param_2 <= thresholds[[param_name]][2]) { 
            #FLAG ASSIGNMENT CORRESPONDING TO STEP 3
            FLAG_param_QC = rep("3", length(TAB_complete[[param_name]])) #stop analysis #dark values are still assigned to FLAG 3 from STEP 1
            FLAG_param_QC[pts.flag3_param] = "3" #FLAG 3 pts assigned in STEP 2a #no-useful command but it is kept just to remind the STEP of assignement
            #FLAG_param_QC[pts.flag4_param]="3" #FLAG 4 pts assigned in STEP 2a
                        
            typeAll[[param_name]] = "3" #catherine
        } else if (rsquared_param_2 > thresholds[[param_name]][2] & rsquared_param_2 <= thresholds[[param_name]][3]) { #begin of STEP 4a: TYPE2_FLAG ASSIGNEMENT
            FLAG_param_QC = rep("2", length(TAB_complete[[param_name]])) # It is a profile without FLAGS 1 # assignment corresponding to STEP 3
            FLAG_param_QC[ (lim+1):length(FLAG_param_QC) ] = "3" #dark values assigned in STEP 1
            FLAG_param_QC[pts.flag3_param] = "3" #FLAG 3 assigned in STEP 2a
            #FLAG_param_QC[pts.flag4_param] = "3" #FLAG 4 assigned in STEP 2a
            FLAG_param_QC[no.cloudy_param[pts.flag3_param_ADJ]] = "3" #FLAG 3 assigned in STEP 4a
            #FLAG_param_QC[no.cloudy_param[pts.flag4_param_ADJ]] = "3" #FLAG 4 assigned in STEP 4a
                        
            typeAll[[param_name]] = "2" #catherine
        } else { #begin of STEP 4b : TYPE1_FLAG ASSIGNEMENT
            FLAG_param_QC = rep("1", length(TAB_complete[[param_name]])) #ASSIGNEMENT CORRESPONDING to STEP 3
            FLAG_param_QC[(lim+1):length(FLAG_param_QC)] = "3" #dark flags assigned in STEP 1
            FLAG_param_QC[pts.flag3_param] = "3" #FLAG 3 assigned in STEP 2a
            #FLAG_param_QC[pts.flag4_param] = "3" #FLAG 4 assigned in STEP 2a
            FLAG_param_QC[no.cloudy_param[pts.flag2_param_ADJ]] = "2" #FLAG 2 assigned in STEP 4b
            FLAG_param_QC[no.cloudy_param[pts.flag3_param_ADJ]] = "3" #FLAG 3 assigned in STEP 4b
            #FLAG_param_QC[no.cloudy_param[pts.flag4_param_ADJ]] = "3" #FLAG 4 assigned in STEP 4b
          
            typeAll[[param_name]] = "1" #catherine
        }
        
        TAB_complete[[ FLAG_NAMES[[param_name]] ]] = FLAG_param_QC   
        TAB_NA[[ FLAG_NAMES[[param_name]] ]] = rep("NA", length(TAB_NA$PRES))
        
        TAB_param = rbind(TAB_complete, TAB_NA) #to create a string of FLAG with same length that NETCDF file
        newdata_All[[ FLAG_NAMES[[param_name]] ]] = TAB_param[order(as.numeric(row.names(TAB_param))),][[ FLAG_NAMES[[param_name]] ]]  #to sort for the initial row order, NA included
        # newdata_All[[ FLAG_NAMES[[param_name]] ]] is the right string of characters to be included in NETCDF file
    }

    return(list("FLAG_380"=newdata_All$FLAG_380_QC, "FLAG_412"=newdata_All$FLAG_412_QC, "FLAG_490"=newdata_All$FLAG_490_QC, 
                "FLAG_PAR"=newdata_All$FLAG_PAR_QC, "type380"=typeAll[["IRR_380"]], "type412"=typeAll[["IRR_412"]], 
                "type490"=typeAll[["IRR_490"]], "typePAR"=typeAll[["PAR"]]))
}
  
#### THE END ####  
