jw15 = read.csv("Data/_raw_data_new/jw15b.csv")
jw16 = read.csv("Data/_raw_data_new/jw16b.csv")
jw17 = read.csv("Data/_raw_data_new/jw17b.csv")
jw18 = read.csv("Data/_raw_data_new/jw18b.csv")
jw19 = read.csv("Data/_raw_data_new/jw19b.csv")

load('Data/data_format_new.rda')
pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
            417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
            732525, 758025, 763050, 843000)
data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
select_id = unique(data_format[,"EID"])
med_select_id = jw15[jw15$key %in% select_id, ]
med_select_id = rbind(med_select_id, jw16[jw16$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw17[jw17$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw18[jw18$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw19[jw19$key %in% select_id, ])
med_select_id$administered_dtm = med_select_id$administered_dtm / 60
med_select_id = med_select_id[!is.na(med_select_id[,1]),]
rownames(med_select_id) = NULL

# Summary Status
all_status = unique(c(jw15$Status, jw16$Status, jw17$Status, jw18$Status, jw19$Status))
status_update = vector(mode = 'list', length = 4)
names(status_update) = c("Start", "Stop", "Continue", "Changed")

status_update[["Start"]]    = c("Given", "New Bag", "Restarted", "MAR Unhold", 
                                "Unheld by provider", "Started During Downtime", 
                                "Override Pull", "Override pull for Anesthesia",
                                "Bolus from Bag", "Medication Applied", "Given by Other",
                                "Given During Downtime", "Bolus")
status_update[["Stop"]]     = c("Stopped", "Held by provider", "Automatically Held", 
                                "MAR Hold", "Anesthesia Discontinued", "Medication Removed")
status_update[["Continue"]] = c("Rate Verify", "Continued from OR", "Handoff",
                                "Continue to Inpatient Floor", "Continue to External Healthcare Facility",
                                "Subsequent Bag")
status_update[["Changed"]]  = c("Rate Change", "Anesthesia Volume Adjustment")


# Simplifying the medication names **** DOUBLE CHECK to make sure we are labelled the correct one
med_key = read.csv('Med_chart.csv', na.strings = "")
colnames(med_key) = c('id', 'hr', 'map', 'onset', 'offset', 'time_check', 'X1')
med_key$id = paste0(" ", med_key$id)
library(stringr)

med_name_simple = paste0(" ", med_select_id$Med_Name_Desc)
hr_map = matrix(0, ncol = 2, nrow = length(med_name_simple))
colnames(hr_map) = c("hr", "map")
for(i in 1:nrow(med_key)) {
    ind = grep(med_key$id[i], med_name_simple)
    if(length(ind) > 0) {
        med_name_simple[ind] = med_key$id[i]
        hr_map[ind, "hr"] = med_key[i, "hr"]
        hr_map[ind, "map"] = med_key[i, "map"]
    }
}

hr_map[hr_map[,1] == "Up", 1] = 1
hr_map[hr_map[,1] == "Down", 1] = -1
hr_map[is.na(hr_map[,1]), 1] = 0
hr_map[hr_map[,2] == "Up", 2] = 1
hr_map[hr_map[,2] == "Down", 2] = -1
hr_map[is.na(hr_map[,2]), 2] = 0

med_select_id = cbind(med_name_simple, hr_map, med_select_id)

med_select_id_sub = med_select_id[,c("key", "administered_dtm", "med_name_simple", "hr", "map",
                                     "Status", "Frequency", "Dose", "Dose_Units",
                                     "Strength", "Med_Name_Desc")]
# Putting in chronological order
for(i in unique(med_select_id_sub$key)) {
    sub_ind = med_select_id_sub[med_select_id_sub$key == i, , drop =F]
    re_ordered_sub = sub_ind[order(sub_ind$administered_dtm), ]
    med_select_id_sub[med_select_id_sub$key == i, ] = re_ordered_sub
}

# making the strength numeric
Strength_num = as.numeric(gsub("\\D", "", med_select_id_sub$Strength))
med_select_id_sub = cbind(med_select_id_sub, Strength_num)
med_select_id_sub$Dose[is.na(med_select_id_sub$Dose)] = 0
med_select_id_sub$Strength_num[is.na(med_select_id_sub$Strength_num)] = 0

# Verifying that the dose * strength = 0 for "Stopped"
print(unique(med_select_id_sub$Dose[med_select_id_sub$Status == "Stopped"]))

# First, make a binary matrix of when a hr or map medication is active
hr_med  = med_select_id_sub[med_select_id_sub$hr  != 0, ]
map_med = med_select_id_sub[med_select_id_sub$map != 0, ] 

hr_continuous   = hr_med[hr_med$Frequency   == "Continuous", ]
hr_single_dose  = hr_med[hr_med$Frequency   != "Continuous", ]
map_continuous  = map_med[map_med$Frequency == "Continuous", ]
map_single_dose = map_med[map_med$Frequency != "Continuous", ]

hr_binary_c = matrix(0, nrow = nrow(data_format), 
                   ncol = length(unique(hr_continuous$med_name_simple)))
colnames(hr_binary_c) = unique(hr_continuous$med_name_simple)
hr_binary_s = matrix(0, nrow = nrow(data_format), 
                     ncol = length(unique(hr_single_dose$med_name_simple)))
colnames(hr_binary_s) = unique(hr_single_dose$med_name_simple)

map_binary_c = matrix(0, nrow = nrow(data_format), 
                     ncol = length(unique(map_continuous$med_name_simple)))
colnames(map_binary_c) = unique(map_continuous$med_name_simple)
map_binary_s = matrix(0, nrow = nrow(data_format), 
                     ncol = length(unique(map_single_dose$med_name_simple)))
colnames(map_binary_s) = unique(map_single_dose$med_name_simple)

# Continuous: Heart Rate Medication --------------------------------------------
for(i in 1:ncol(hr_binary_c)) {
    for(j in unique(hr_continuous$key)) {
        sub_dat = hr_continuous[hr_continuous$key == j, , drop = F]
        med_specific = sub_dat[sub_dat$med_name_simple == colnames(hr_binary_c)[i], , drop = F]
        doses = as.matrix(table(med_specific$Dose_Units))
        dose_name = rownames(doses)[which.max(doses)]
        
        if(length(unique(med_specific$Dose_Units)) > 1) {
            print(paste0(i, ": ", colnames(hr_binary_c)[i], " Patient: ", j))
            print(dose_name)
            print(unique(med_specific$Dose_Units))
            
            # INITIAL FIX! Come back later! ***********************************
            # mean_dose          = mean(med_specific$Dose[med_specific$Dose_Units == dose_name])
            # mean_dose_strength = mean(med_specific$Strength_num[med_specific$Dose_Units == dose_name])
            # med_specific$Dose[med_specific$Dose_Units != dose_name] = mean_dose
            # med_specific$Strength_num[med_specific$Dose_Units != dose_name] = mean_dose_strength
            med_specific$Dose[med_specific$Dose_Units != dose_name] = 0
        }
        
        slot = matrix(0, ncol = 2, nrow = sum(data_format[,"EID"] == j))
        slot[,1] = data_format[data_format[,"EID"] == j, "time"]
        
        if(nrow(med_specific > 0)) {
            for (k in 1:nrow(slot)) {
                main_effect = 0
                if(k == 1) {
                    # First entry
                    sub_med = med_specific[med_specific$administered_dtm >= 0 & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                } else {
                    # Last entry
                    sub_med = med_specific[med_specific$administered_dtm >= slot[k-1] & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                }
                slot[k, 2] = main_effect
            }
        }
        
        hr_binary_c[data_format[,"EID"] == j, i] = c(slot[,2])
        
    }
}

# Continuous: MAP Medication --------------------------------------------------
for(i in 1:ncol(map_binary_c)) {
    for(j in unique(map_continuous$key)) {
        sub_dat = map_continuous[map_continuous$key == j, , drop = F]
        med_specific = sub_dat[sub_dat$med_name_simple == colnames(map_binary_c)[i], , drop = F]
        doses = as.matrix(table(med_specific$Dose_Units))
        dose_name = rownames(doses)[which.max(doses)]
        
        if(length(unique(med_specific$Dose_Units)) > 1) {
            print(paste0(i, ": ", colnames(map_binary_c)[i], " Patient: ", j))
            print(dose_name)
            print(unique(med_specific$Dose_Units))
            
            # INITIAL FIX! Come back later! ***********************************
            # mean_dose          = mean(med_specific$Dose[med_specific$Dose_Units == dose_name])
            # mean_dose_strength = mean(med_specific$Strength_num[med_specific$Dose_Units == dose_name])
            # med_specific$Dose[med_specific$Dose_Units != dose_name] = mean_dose
            # med_specific$Strength_num[med_specific$Dose_Units != dose_name] = mean_dose_strength
            med_specific$Dose[med_specific$Dose_Units != dose_name] = 0
        }
        
        slot = matrix(0, ncol = 2, nrow = sum(data_format[,"EID"] == j))
        slot[,1] = data_format[data_format[,"EID"] == j, "time"]
        
        if(nrow(med_specific > 0)) {
            for (k in 1:nrow(slot)) {
                main_effect = 0
                if(k == 1) {
                    # First entry
                    sub_med = med_specific[med_specific$administered_dtm >= 0 & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                } else {
                    # Last entry
                    sub_med = med_specific[med_specific$administered_dtm >= slot[k-1] & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                }
                slot[k, 2] = main_effect
            }
        }
        
        map_binary_c[data_format[,"EID"] == j, i] = c(slot[,2])
        
    }
}

# Single Dose: Heart Rate Medication -------------------------------------------
for(i in 1:ncol(hr_binary_s)) {
    for(j in unique(hr_single_dose$key)) {
        sub_dat = hr_single_dose[hr_single_dose$key == j, , drop = F]
        med_specific = sub_dat[sub_dat$med_name_simple == colnames(hr_binary_s)[i], , drop = F]
        doses = as.matrix(table(med_specific$Dose_Units))
        dose_name = rownames(doses)[which.max(doses)]
        
        if(length(unique(med_specific$Dose_Units)) > 1) {
            print(paste0(i, ": ", colnames(hr_binary_s)[i], " Patient: ", j))
            print(dose_name)
            print(unique(med_specific$Dose_Units))
            
            # INITIAL FIX! Come back later! ***********************************
            # mean_dose          = mean(med_specific$Dose[med_specific$Dose_Units == dose_name])
            # mean_dose_strength = mean(med_specific$Strength_num[med_specific$Dose_Units == dose_name])
            # med_specific$Dose[med_specific$Dose_Units != dose_name] = mean_dose
            # med_specific$Strength_num[med_specific$Dose_Units != dose_name] = mean_dose_strength
            med_specific$Dose[med_specific$Dose_Units != dose_name] = 0
        }
        
        slot = matrix(0, ncol = 2, nrow = sum(data_format[,"EID"] == j))
        slot[,1] = data_format[data_format[,"EID"] == j, "time"]
        
        if(nrow(med_specific > 0)) {
            for (k in 1:nrow(slot)) {
                main_effect = 0
                if(k == 1) {
                    # First entry
                    sub_med = med_specific[med_specific$administered_dtm >= 0 & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                } else {
                    # Last entry
                    sub_med = med_specific[med_specific$administered_dtm >= slot[k-1] & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                }
                slot[k, 2] = main_effect
            }
        }
        
        hr_binary_s[data_format[,"EID"] == j, i] = c(slot[,2])
        
    }
}

# Single Dose: MAP Medication --------------------------------------------------
for(i in 1:ncol(map_binary_s)) {
    for(j in unique(map_single_dose$key)) {
        sub_dat = map_single_dose[map_single_dose$key == j, , drop = F]
        med_specific = sub_dat[sub_dat$med_name_simple == colnames(map_binary_s)[i], , drop = F]
        doses = as.matrix(table(med_specific$Dose_Units))
        dose_name = rownames(doses)[which.max(doses)]
        
        if(length(unique(med_specific$Dose_Units)) > 1) {
            print(paste0(i, ": ", colnames(map_binary_s)[i], " Patient: ", j))
            print(dose_name)
            print(unique(med_specific$Dose_Units))
            
            # INITIAL FIX! Come back later! ***********************************
            # mean_dose          = mean(med_specific$Dose[med_specific$Dose_Units == dose_name])
            # mean_dose_strength = mean(med_specific$Strength_num[med_specific$Dose_Units == dose_name])
            # med_specific$Dose[med_specific$Dose_Units != dose_name] = mean_dose
            # med_specific$Strength_num[med_specific$Dose_Units != dose_name] = mean_dose_strength
            med_specific$Dose[med_specific$Dose_Units != dose_name] = 0
        }
        
        slot = matrix(0, ncol = 2, nrow = sum(data_format[,"EID"] == j))
        slot[,1] = data_format[data_format[,"EID"] == j, "time"]
        
        if(nrow(med_specific > 0)) {
            for (k in 1:nrow(slot)) {
                main_effect = 0
                if(k == 1) {
                    # First entry
                    sub_med = med_specific[med_specific$administered_dtm >= 0 & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                } else {
                    # Last entry
                    sub_med = med_specific[med_specific$administered_dtm >= slot[k-1] & med_specific$administered_dtm < slot[k], , drop = F]
                    if(nrow(sub_med) > 0) {
                        time_scale = c(sub_med$administered_dtm, slot[k, 1])
                        time_scale = diff(time_scale)
                        main_effect = time_scale * (sub_med$Dose * sub_med$Strength_num)
                        main_effect = mean(main_effect)
                    }
                }
                slot[k, 2] = main_effect
            }
        }
        
        map_binary_s[data_format[,"EID"] == j, i] = c(slot[,2])
        
    }
}








