.libPaths(c("C:/Users/rdeo0001/DOcuments/Packages", .libPaths()))
.libPaths()

library(NicheMapR)
library(dplyr)
library(plyr)

####

##adonis
#Ww_g <- 1.806      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 23.72    # LTset (10%)
#T_F_max <- 29.30    # UTset (90%)
#T_B_min <- 23.72    # LTset (10%)
#T_RB_min <- 23.72   # LTset (10%)
#T_pref <- 26.69     # mean selected body temperature in a thermal gradient
#CT_max <- 38.35
#CT_min <- 10.74
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 10


##amicula
#Ww_g <- 1.309      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 25.0777 # LTset (25%)
#T_F_max <- 29.0425    # UTset (90%)
#T_B_min <- 22.5997 # LTset (10%)
#T_RB_min <- 22.5997 # LTset (10%)
#T_pref <- 26.23351     # mean selected body temperature in a thermal gradient
#CT_max <- 38.48947368
#CT_min <- 12.61052632
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2

##caligula
#Ww_g <- 3.011      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 26.4391 # LTset (25%)
#T_F_max <- 29.4097    # UTset (90%)
#T_B_min <- 25.69645 # LTset (10%)
#T_RB_min <- 25.69645 # LTset (10%)
#T_pref <- 27.39243     # mean selected body temperature in a thermal gradient
#CT_max <- 38.8
#CT_min <- 7.87
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2

##couperi
#Ww_g <- 1.3585      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 25.8716 # LTset (25%)
#T_F_max <- 29.19765    # UTset (90%)
#T_B_min <- 25.10405 # LTset (10%)
#T_RB_min <- 25.10405 # LTset (10%)
#T_pref <- 27.09375     # mean selected body temperature in a thermal gradient
#CT_max <- 38.4041
#CT_min <- 11.2171
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2


##delicata Brisbane
Ww_g <- 1.558      # mean weight
shape <- 3          # cylinder body shape
T_F_min <- 24.3596 # LTset (25%)
T_F_max <- 28.9685    # UTset (90%)
T_B_min <- 20.2628 # LTset (10%)
T_RB_min <- 20.2628 # LTset (10%)
T_pref <- 24.6763     # mean selected body temperature in a thermal gradient
CT_max <- 38.89130435
CT_min <- 8.369565217
diurn <- 1
shade_seek <- 1
burrow <- 1
maxdepth <- 2


##delicata nsw
#Ww_g <- 1.15      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 24.59166 # LTset (25%)
#T_F_max <- 30.4441    # UTset (90%)
#T_B_min <- 22.2121 # LTset (10%)
#T_RB_min <- 22.2121 # LTset (10%)
#T_pref <- 26.68834     # mean selected body temperature in a thermal gradient
#CT_max <- 39.1825
#CT_min <- 8.7955
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2

##delicata Townsville
#Ww_g <- 1.229       # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 25.63985 # LTset (25%)
#T_F_max <- 29.22455    # UTset (90%)
#T_B_min <- 24.8717 # LTset (10%)
#T_RB_min <- 24.8717 # LTset (10%)
#T_pref <- 26.61931     # mean selected body temperature in a thermal gradient
#CT_max <- 39.08421053
#CT_min <- 8.115789474
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2


##guichenoti
#Ww_g <- 1.6542      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 24.513 # LTset (25%)
#T_F_max <- 29.9646    # UTset (90%)
#T_B_min <- 21.7872 # LTset (10%)
#T_RB_min <- 21.7872 # LTset (10%)
#T_pref <- 26.24989     # mean selected body temperature in a thermal gradient
#CT_max <- 39.44
#CT_min <- 7.9426
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2

##similis
#Ww_g <- 1.192      # mean weight
#shape <- 3          # cylinder body shape
#T_F_min <- 26.1497 # LTset (25%)
#T_F_max <- 29.22715    # UTset (90%)
#T_B_min <- 25.19115 # LTset (10%)
#T_RB_min <- 25.19115 # LTset (10%)
#T_pref <- 26.72222     # mean selected body temperature in a thermal gradient
#CT_max <- 38.5
#CT_min <- 11.095
#diurn <- 1
#shade_seek <- 1
#burrow <- 1
#maxdepth <- 2


# DEM
library(geodata)
elevation <- geodata::worldclim_global(var = 'elev', res = 2.5, path=tempdir())


##load localities
data <- read.csv("aus_coord_0.25_b.csv", header = TRUE)


#loop through localities

for(i in 1:nrow(data)){
  
  micro1 <- micro_terra(loc = c(data$Long[i], data$Lat[i]),
                        ystart = 2015, yfinish = 2015,
                        scenario = 0, 
                        timeinterval = 365,
                        DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200),
                        Usrhyt = 0.01,
                        elevation = elevation,
                        #terrain = 1, 
                        #elevatr = 1,
                        #runmoist = 1,
                        #microclima.LAI = 0,
                        #soilgrids = 1,
                        minshade = 0, maxshade = 90,runshade = 1)
  
  micro2 <- micro_terra(loc = c(data$Long[i], data$Lat[i]),
                        ystart = 2015, yfinish = 2015,
                        scenario = 2,
                        timeinterval = 365,
                        DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200),
                        Usrhyt = 0.01,
                        elevation = elevation,
                        #terrain = 1,
                        #elevatr = 1,
                        #runmoist = 1,
                        #soilgrids = 1,
                        #microclima.LAI = 0,
                        minshade = 0, maxshade = 90,runshade = 1)
  
  micro3 <- micro_terra(loc = c(data$Long[i], data$Lat[i]),
                        ystart = 2015, yfinish = 2015,
                        scenario = 4,
                        timeinterval = 365,
                        DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200),
                        Usrhyt = 0.01,
                        elevation = elevation,
                        #terrain = 1,
                        #elevatr = 1,
                        #runmoist = 1,
                        #soilgrids = 1,
                        #microclima.LAI = 0,
                        minshade = 0, maxshade = 90,runshade = 1)
  
  
  for(k in 1:3){
    if(k == 1){
      micro <- micro1
    }
    if(k == 2){
      micro <- micro2
    }
    if(k == 3){
      micro <- micro3
    }
    
    # run the ectotherm model
    ecto <- ectotherm(Ww_g = Ww_g, 
                      shape = shape, 
                      CT_max = CT_max, 
                      CT_min = CT_min, 
                      T_F_min = T_F_min, 
                      T_F_max = T_F_max, 
                      T_pref = T_pref, 
                      T_B_min = T_B_min, 
                      T_RB_min = T_RB_min, 
                      diurn = diurn, 
                      shade_seek = shade_seek, 
                      burrow = burrow, 
                      maxdepth = maxdepth)
    
    # output tables
    environ<-as.data.frame(ecto$environ)
    output <- environ
    TC <- environ$TC
    
    # predicted sprint speed
    #output$Speed <- 1.619435 + -0.572094*TC + 0.057061*TC^2 + -0.001105*TC^3 ##caligula
    #output$Speed <- -28.9198779  + 5.2640275*TC + -0.3358946*TC^2 + 0.0098772*TC^3 + -0.0001089*TC^4 ##caligula
    output$Speed <- -8.01963305 + 1.81769600*TC + -0.14154868*TC^2 + 0.00530536*TC^3 ##delicata brisbane
    + -0.00006991*TC^4
    #output$Speed <- -0.9532568  + -0.0238336*TC + 0.0238542*TC^2 + -0.0005552*TC^3 ##delicata nsw
    #output$Speed <- -21.0772442 + 4.7676261*TC + -0.3442419*TC^2 + 0.0111649*TC^3
    #+ -0.0001308*TC^4    ##delicata twv
    #output$Speed <- -2.45864520  + 0.44217708*TC + -0.01942396*TC^2 + 0.00098532*TC^3 
    #+ -0.00001816*TC^4 ##guichenoti
    #output$Speed <- -78.8746526  + 15.2129716*TC + -1.0246042*TC^2 + 0.0299924*TC^3 
    #+ -0.00031846*TC^4 ##similis
    
    
    # predicted evaporative water loss mg H2O/h
    #output$EWL <- -0.90262 + 0.1763*TC   ##caligula
    #output$EWL <- -3.4616 + 0.3061*TC     ##couperi
    output$EWL <- 15.06983 + -1.22115*TC + 0.03309*TC^2  #delicata brisbane
    #output$EWL <- -1.0253 + 0.1734*TC      ## delicata nsw
    #output$EWL <- -3.1709 + 0.3548*TC      ## delicata twv
    #output$EWL <- 1.0615 + 0.1265*TC        ## guichenoti
    #output$EWL <- -3.1405 + 0.3047*TC        ##similis
    
    # predicted dehydration (percentage of water loss in relation to body mass)
    output$dehydration <- (output$EWL*100/Ww_g)/1000
    
    # predicted metabolic rate in vCO2
    #output$MR <- -0.31638 + 0.02738*TC   ##caligula
    #output$MR <- -0.16619 + 0.01692*TC   ##couperi
    output$MR <- -0.20862 + 0.02085*TC
    #output$MR <- -0.1001 + 0.0101*TC
    #output$MR <- -0.16655 + 0.01429*TC    ##delicata twv
    #output$MR <- 0.1821891 + -0.0183657*TC + 0.0007012*TC^2  ##guichenoti
    #output$MR <- -0.17253 + 0.01381*TC    ##similis
    
    # energy consumption
    # based on 1 g of cricket dry mass (67% protein, 10% fat, 18% carbohydrate) - 1 g protein = 4 cal, fat = 9 cal, carb = 4 cal
    # 1 g dry cricket = 4 x 0.67 + 9 x 0.1 + 4 x 0.18 = 4.3 cal
    # 1 cal = 4.18 kJ
    # FMR, in units of ml CO2 day–1, was converted to energy equivalents by multiplication by 25·8 J ml–1 CO2 (Peterson et al. 2002. Func Ecology. Metabolic costs of growth in free-living Garter Snakes and the energy budgets of ectotherms)
    
    
    # number of crickets need per hour
    output$energ_comsump <- (output$MR * 25.8)*2.5/ # energy consumption in kJ/h - Aerobic scope ~4x resting metabolic rate
      17.974               # 1 g of dry cricket in kJ
    
    
    ## save coordinates ##
    output$id <- data$id[i]
    output$Lat <- data$Lat[i]
    output$Long <- data$Long[i]
    
    
    if(i == 1){
      if(k == 1){
        output.frame1 <- output
        
      }
      if(k == 2){
        output.frame2 <- output
        
      }
      if(k == 3){
        output.frame3 <- output
        
      }    
    }else{
      if(k == 1){
        output.frame1 <- rbind(output.frame1, output)
        
      }
      if(k == 2){
        output.frame2 <- rbind(output.frame2, output)
        
      }
      if(k == 3){
        output.frame3 <- rbind(output.frame3, output)
        
      }     
    }
  }
}


## Final dataset current climate - output.frame1 ##

dt.current.total <- ddply(output.frame1, "id", summarise,
                    TC = mean(TC),
                    maxTC = max(TC),
                    speed = mean(Speed),
                    met = mean(MR),
                    ewl = mean(EWL),
                    dehy = mean(dehydration),
                    dehy = max(dehydration),
                    energ_comsump = mean(energ_comsump),
                    shade = mean(SHADE),
                    depth = mean(DEP)
                    )

dt.current.act <- output.frame1 %>% filter(ACT == 2) %>%
                    ddply("id", summarise,
                    TC.act = mean(TC),
                    maxTC.act = max(TC),
                    speed.act = mean(Speed),
                    met.act = mean(MR),
                    ewl.act = mean(EWL),
                    dehy.act = mean(dehydration),
                    energ_comsump.act = mean(energ_comsump),
                    shade.act = mean(SHADE),
                    depth.act = mean(DEP))

dt.h.act.current <- output.frame1 %>%
                      #group_by(id) %>%
                      ddply("id", summarise,
                                total_activity_hours = sum(ACT == 2),
                                dehydrated.10 = sum(dehydration >= 0.20),
                                ctmax = sum(TC >= CT_max))

final.current <- dt.current.total %>%
                  left_join(dt.current.act, by = "id") %>%
                  left_join(dt.h.act.current, by = "id")


## Final dataset +2 C climate - outoutput.frame2 ##

dt.plus2.total <- ddply(output.frame2, "id", summarise,
                        TC = mean(TC),
                        maxTC = max(TC),
                        speed = mean(Speed),
                        met = mean(MR),
                        ewl = mean(EWL),
                        dehy = mean(dehydration),
                        dehy = max(dehydration),
                        energ_comsump = mean(energ_comsump),
                        shade = mean(SHADE),
                        depth = mean(DEP)
)

dt.plus2.act <- output.frame2 %>% filter(ACT == 2) %>%
                ddply("id", summarise,
                      TC.act = mean(TC),
                      maxTC.act = max(TC),
                      speed.act = mean(Speed),
                      met.act = mean(MR),
                      ewl.act = mean(EWL),
                      dehy.act = mean(dehydration),
                      energ_comsump.act = mean(energ_comsump),
                      shade.act = mean(SHADE),
                      depth.act = mean(DEP))

dt.h.act.plus2 <- output.frame2 %>%
                    #group_by(id) %>%
                    ddply("id", summarise,
                              total_activity_hours = sum(ACT == 2),
                              dehydrated.10 = sum(dehydration >= 0.20),
                              ctmax = sum(TC >= CT_max))

final.plus2 <- dt.plus2.total %>%
  left_join(dt.plus2.act, by = "id") %>%
  left_join(dt.h.act.plus2, by = "id")


## Final dataset + 4 C climate - outoutput.frame3 ##

dt.plus4.total <- ddply(output.frame3, "id", summarise,
                              TC.act = mean(TC),
                              maxTC.act = max(TC),
                              speed.act = mean(Speed),
                              met.act = mean(MR),
                              ewl.act = mean(EWL),
                              dehy.act = mean(dehydration),
                              energ_comsump.act = mean(energ_comsump),
                              shade.act = mean(SHADE),
                              depth.act = mean(DEP))
                        
dt.plus4.act <- output.frame3 %>% filter(ACT == 2) %>%
                  ddply("id", summarise,
                        TC.act = mean(TC),
                        maxTC.act = max(TC),
                        speed.act = mean(Speed),
                        met.act = mean(MR),
                        ewl.act = mean(EWL),
                        dehy.act = mean(dehydration),
                        energ_comsump.act = mean(energ_comsump),
                        shade.act = mean(SHADE),
                        depth.act = mean(DEP))

dt.h.act.plus4 <- output.frame3 %>%
                      #group_by(id) %>%
                      ddply("id", summarise,
                                total_activity_hours = sum(ACT == 2),
                                dehydrated.10 = sum(dehydration >= 0.20),
                                ctmax = sum(TC >= CT_max))

final.plus4 <- dt.plus4.total %>%
                  left_join(dt.plus4.act, by = "id") %>%
                  left_join(dt.h.act.plus4, by = "id")

write.csv(final.current, "final_current.csv")
write.csv(final.plus2, "final_plus2.csv")
write.csv(final.plus4, "final_plus4.csv")

