load("Data/long_data_agg.rda") # long_data_agg
load('Data/all_keys.rda')      # all_keys
load('Data/cov_info.rda')      # cov_info
load('Data/timing_issues.rda') # timing_issues

# ------------------------------------------------------------------------------
# (1) Filtering people based on undesirable characteristics --------------------
# ------------------------------------------------------------------------------

# Removing time issues
all_keys_temp = all_keys[!(all_keys %in% timing_issues)]

# Removing pacing patients
# pace_info = read.csv('Data/_raw_data_new/jw_pacemaker.csv')
# no_pace_patient = pace_info$key[is.na(pace_info$pacemaker_attention) 
#                                 & is.na(pace_info$pacemaker_present)]
# all_keys_temp = all_keys_temp[all_keys_temp %in% no_pace_patient]

# Removing patients who died in the ICU
deaths = cov_info[cov_info[,"icu_death"] == 1, "key"]
all_keys_temp = all_keys_temp[!(all_keys_temp %in% deaths)]

# Isolating focus to patients who have had a long enough stay
max_times = matrix(nrow = length(all_keys_temp), ncol = 2)
max_times[,1] = all_keys_temp

for(i in 1:length(long_data_agg)) {
    if(long_data_agg[[i]]$key %in% all_keys_temp) {
        if(nrow(long_data_agg[[i]]$covariates) > 1) {
            max_times[max_times[,1] == long_data_agg[[i]]$key, 2] = 
                max(long_data_agg[[i]]$covariates[, "time"], na.rm = T)
        } else {
            max_times[max_times[,1] == long_data_agg[[i]]$key, 2] = 0
        }
    }
}

# Convert to hours and only take patients with >= 10 and <= 48 length stay
max_times[,2] = max_times[,2] / 60
max_times = max_times[max_times[,2] >= 10 & max_times[,2] <= 48, ]
all_keys_temp = c(max_times[,1])

# Temporarily subsetting long_data_agg -----
long_data_agg_sub = vector(mode = 'list', length = length(all_keys_temp))
ldas = 1
for(i in 1:length(long_data_agg)) {
    if(long_data_agg[[i]]$key %in% all_keys_temp) {
        long_data_agg_sub[[ldas]] = long_data_agg[[i]]
        ldas = ldas + 1
    }
}
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key != all_keys_temp[i]) print("wrong order")
}
# ------------------------------------------

# Choose patients with enough observations of hr and map 
enough_dat = rep(0, length(all_keys_temp))
for(i in 1:length(enough_dat)) {
    temp = long_data_agg_sub[[i]]$covariates
    hr_ratio = sum(!is.na(temp[,"hr"])) / nrow(temp)
    map_ratio = sum(!is.na(temp[,"map"])) / nrow(temp)
    
    if(hr_ratio > 0.6 & map_ratio > 0.6) {
        enough_dat[i] = 1
    }
}

all_keys_temp = all_keys_temp[enough_dat == 1]

# ------------------------------------------------------------------------------
# (2) Adding new variables such as RBC rule ------------------------------------
# ------------------------------------------------------------------------------

# Get patients that had at least 3 RBC at some point in their stay
min_three_RBC = NULL
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key %in% all_keys_temp) {
        if(tail(long_data_agg_sub[[i]]$covariates[,"n_RBC"], 1) >= 3) {
            min_three_RBC = c(min_three_RBC, 1)
        } else {
            min_three_RBC = c(min_three_RBC, 0)
        }
    }
}

bleed_pat = all_keys_temp[min_three_RBC == 1]
stable_pat = all_keys_temp[min_three_RBC == 0]

# Formatting the existing data into one data set
data_format = NULL
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key %in% all_keys_temp) {
        print(long_data_agg_sub[[i]]$key)
        data_format = rbind(data_format, long_data_agg_sub[[i]]$covariates)
    }
}

# Add a column to the df
data_format = cbind(data_format, rep(0, nrow(data_format)), rep(0, nrow(data_format)))
colnames(data_format) = c('EID', 'time', 'temperature', 'hemo', 'map', 'hr', 'lactate', 'RBC',
                          'n_labs', 'n_RBC', 'RBC_rule', 'clinic_rule')

# Adding the rule that a patient is bleeding if >= 3 in 12hrs or >= 6 in 24hrs
for(i in 1:length(bleed_pat)) {
    sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
    
    # Check in any 12 hour period
    max_time = tail(sub_dat[,"time"], 1)
    when_rbc = c(1, which(diff(sub_dat[,"n_RBC"]) != 0))
    
    for(j in 1:length(when_rbc)) {
        s_time = sub_dat[when_rbc[j], "time"]
        e_time_12 = s_time + 720
        e_time_24 = s_time + 1440
        RBC_diff_12 = RBC_diff_24 = 0
        
        if (e_time_12 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
            RBC_diff_12 = sub_dat[ind_12, "n_RBC"] - sub_dat[s_ind, "n_RBC"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_12 = sub_dat[e_ind, "n_RBC"] - sub_dat[s_ind, "n_RBC"]
        }
        if (e_time_24 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_24 = order(abs(sub_dat[,"time"] - e_time_24))[1]
            RBC_diff_24 = sub_dat[ind_24, "n_RBC"] - sub_dat[s_ind, "n_RBC"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_24 = sub_dat[e_ind, "n_RBC"] - sub_dat[s_ind, "n_RBC"]
        }
        
        if(RBC_diff_12 >=3 | RBC_diff_24 >= 6) {
            data_format[data_format[,"EID"] == bleed_pat[i], "RBC_rule"] = 1
            print(paste0("Bleeding: ", i))
            break
        }
        
    }
    
}

# Adding the clinical rule
clinical_lab = c(109125, 111375, 133750, 156325, 165725, 195475, 198975, 208100, 327375, 360525,
                 431400, 467200, 494300, 531650, 533825, 543625, 588100, 622450, 633100, 697600,
                 727675, 750900, 758025, 781875, 785950, 801300, 820775, 827350, 841925, 843000)

data_format[data_format[,"EID"] %in% c(111375, 133750, 327375, 431400, 
                                       531650, 633100, 697600, 727675,
                                       758025, 820775, 827350, 841925, 843000), "clinic_rule"] = 1 # bleed
data_format[data_format[,"EID"] %in% c(109125, 156325, 165725, 195475, 198975,
                                       208100, 360525, 467200, 494300, 533825,
                                       543625, 588100, 622450, 750900, 781875,
                                       785950, 801300), "clinic_rule"] = -1 # no bleed

# ------------------------------------------------------------------------------
# (3) Add the new RBC transfusions ---------------------------------------------
# ------------------------------------------------------------------------------
transfusions = read.csv('Data/_raw_data_new/jw_transfusions.csv')
transfusions = transfusions[transfusions$OrderedProduct == "RBC", ]

times = transfusions[,c("key",
                        "transfus_order_datetime", 
                        "transfus_admin_start_datetime", 
                        "OrderedProduct", 
                        "Order_Description")]

data_keys = unique(data_format[,"EID"])
times = times[times$key %in% data_keys, ]
times = times[times$transfus_order_datetime > 0, ]
times$transfus_order_datetime = times$transfus_order_datetime / 60
times$transfus_admin_start_datetime = times$transfus_admin_start_datetime / 60
times = times[!is.na(times$key), ]

RBC_ordered = 0
RBC_admin = 0
data_format = cbind(cbind(data_format, RBC_ordered), RBC_admin)

for(i in 1:length(data_keys)) {
    print(i)
    if(data_keys[i] %in% times$key) {
        sub_dat = data_format[data_format[,"EID"] == data_keys[i], ]
        sub_rbc = times[times$key == data_keys[i], ]
        for(j in 1:(nrow(sub_dat) - 1)) {
            timed_rbc_order = which(sub_rbc$transfus_order_datetime > sub_dat[j,"time"] & sub_rbc$transfus_order_datetime <= sub_dat[j+1,"time"])
            if(length(timed_rbc_order) > 0) {
                sub_dat[j+1, "RBC_ordered"] = length(timed_rbc_order)
            }
            
            timed_rbc_admin = which(sub_rbc$transfus_admin_start_datetime > sub_dat[j,"time"] & sub_rbc$transfus_admin_start_datetime <= sub_dat[j+1,"time"])
            if(length(timed_rbc_admin) > 0) {
                sub_dat[j+1, "RBC_admin"] = length(timed_rbc_admin)
            }
        }
        
        data_format[data_format[,"EID"] == data_keys[i], ] = sub_dat
    }
}


# Adding n_RBC_ordered and n_RBC_admin
n_RBC_ordered = 0
n_RBC_admin = 0
data_format = cbind(cbind(data_format, n_RBC_ordered), n_RBC_admin)

for(i in unique(data_format[,'EID'])){
    data_format[data_format[,'EID'] == i, 'n_RBC_ordered'] = cumsum(data_format[data_format[,'EID'] == i, 'RBC_ordered'])
    data_format[data_format[,'EID'] == i, 'n_RBC_admin'] = cumsum(data_format[data_format[,'EID'] == i, 'RBC_admin'])
}

# ------------------------------------------------------------------------------
# Removing pacing patients based on own heuristic ------------------------------
# ------------------------------------------------------------------------------
pace_ids = NULL
for(i in unique(data_format[,"EID"])) {
    sub_id_hr = data_format[data_format[,"EID"] == i, "hr"]
    diff_hr = diff(sub_id_hr)
    conseq_same = which(diff_hr == 0)
    if(sum(diff(conseq_same) == 1) > 6) {
        pace_ids = c(pace_ids, i)
    }
}

data_format = data_format[!(data_format[,"EID"] %in% pace_ids), ]

# ------------------------------------------------------------------------------
# (4) Saving the data_format with subset of IDs --------------------------------
# ------------------------------------------------------------------------------
set.seed(2023)
load('Data/curr_id.rda')
rbc_rule = unique(data_format[data_format[,"RBC_rule"] == 1, "EID"])
no_rbc_rule = unique(data_format[data_format[,"RBC_rule"] == 0, "EID"])
clinic_rule = unique(data_format[data_format[,"clinic_rule"] != 0, "EID"])

curr_id = curr_id[curr_id %in% data_format[,"EID"]]

curr_id = unique(c(curr_id, rbc_rule, clinic_rule))

add_id  = no_rbc_rule[!(no_rbc_rule %in% curr_id)]

curr_id = c(curr_id, sample(add_id, 400 - length(curr_id)))
curr_id = sort(unique(curr_id))

data_format = data_format[data_format[,"EID"] %in% curr_id, ]

# ------------------------------------------------------------------------------
# (5) Remove any eroneous data points ------------------------------------------
# ------------------------------------------------------------------------------
for(i in 1:nrow(data_format)) {
    # hemo
    if(!is.na(data_format[i,'hemo'])) {
        if(data_format[i,'hemo'] < 0 | data_format[i,'hemo'] > 20) {
            data_format[i, 'hemo'] = NA
        } 
    }
    # hr
    if(!is.na(data_format[i,'hr'])) {
        if(data_format[i,'hr'] < 20 | data_format[i,'hr'] > 200) {
            data_format[i, 'hr'] = NA
        } 
    }
    # map
    if(!is.na(data_format[i,'map'])) {
        if(data_format[i,'map'] < 25 | data_format[i,'map'] > 150) {
            data_format[i, 'map'] = NA
        } 
    }

}

save(data_format, file = 'Data/data_format_new2.rda')
# Quick check for no pacing
pdf("Data/quick_check.pdf")
par(mfrow=c(4,1))
for(i in unique(data_format[,"EID"])) {
    plot(data_format[data_format[,"EID"] == i, "time"], data_format[data_format[,"EID"] == i, "hr"],
         main = i)
}
dev.off()

# ------------------------------------------------------------------------------
# (6) Filter based on level of care --------------------------------------------
# ------------------------------------------------------------------------------
# load("Data/data_format_FULL_48hr_update_RBC_sub.rda")
# # Removing pacing patients
# pace_id = c(18075, 108825, 110750, 125025, 173750, 260100, 304700, 307225, 310100,
#             382450, 429375, 516150, 533075, 666750, 677225, 732525, 763050, 767500,
#             769025, 777175, 794900, 799125, 819225)
# data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
# curr_id = unique(data_format[,'EID'])
# save(curr_id, file = 'Data/curr_id.rda')

level_of_care = read.csv('Data/_raw_data_new/jw_patient_level_of_care.csv')
load('Data/data_format_new2.rda')
care_time_df = matrix(nrow = nrow(data_format), ncol = 3)
care_time_df[,1] = data_format[,'EID']
care_time_df[,2] = data_format[,'time']

level_of_care_patient = level_of_care[level_of_care$key %in% data_format[,"EID"],]
level_of_care_patient$level_of_care_datetime = level_of_care_patient$level_of_care_datetime / 60

for(i in unique(care_time_df[,1])) {
    sub_level = level_of_care_patient[level_of_care_patient$key == i, ]
    sub_df    = data_format[data_format[,"EID"] == i, ]
    sub_care  = care_time_df[care_time_df[,1] == i, ]
    
    level_times = sub_level$level_of_care_datetime
    
    for(j in 1:nrow(sub_df)) {
        time_j = sub_df[j,"time"]
        
        if(time_j < level_times[1]) {
            ind = which(level_times == level_times[1])
            sub_care[j,3] = sub_level$patient_level_of_care[ind]
        } else if(time_j > tail(level_times,1)) {
            ind = which(level_times == tail(level_times,1))
            sub_care[j,3] = sub_level$patient_level_of_care[ind]
        } else {
            end_pts = c(max(level_times[level_times < time_j]), 
                        min(level_times[level_times > time_j]))
            ind = max(which(level_times == end_pts[1]))
            sub_care[j,3] = sub_level$patient_level_of_care[ind]
        }
    }
    care_time_df[care_time_df[,1] == i, 3] = sub_care[,3]
}
save(care_time_df, file = 'Data/care_time_df.rda')

