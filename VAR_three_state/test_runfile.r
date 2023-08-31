library(RcppArmadillo)
library(RcppDist)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# load('test_post.rda')
# load('test_list.rda')
Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)

simulation = T
data_format = NULL
if(simulation) {
    load('Data/use_data1_1.rda')
    data_format = use_data
    trialNum = 1
} else {
    load('Data/data_format_new.rda')
    pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
                417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
                732525, 758025, 763050, 843000)
    data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
    trialNum = 7 # CHANGE THIS EVERY TIME **********************
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

otype = !is.na(Y[, c('hemo','hr','map','lactate')])
colnames(otype) = c('hemo','hr','map','lactate')

B = list()
for(i in EIDs){
    
    if(simulation) {
        B[[i]] = data_format[data_format[,'EID']==as.numeric(i), "b_true", drop=F]
    } else {
        b_temp = matrix( 1, sum(Y[,'EID']==as.numeric(i)), 1)
        # b_length = nrow(b_temp)
        # b_temp[(b_length-5):b_length, ] = 1
        B[[i]] = b_temp
    }
}

# nu_R = prop_R_test$nu_R
# psi_R = prop_R_test$psi_R
# Y = prop_R_test$Y
# Dn = prop_R_test$Dn
# Xn = prop_R_test$Xn
# par = prop_R_test$par
# par_index = prop_R_test$par_index
# EIDs = prop_R_test$EIDs


test_fnc()
