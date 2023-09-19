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

med_select_id_sub = med_select_id_sub[!is.na(med_select_id_sub$administered_dtm), ]

# making the strength numeric
Strength_num = as.numeric(gsub("\\D", "", med_select_id_sub$Strength))
med_select_id_sub = cbind(med_select_id_sub, Strength_num)
med_select_id_sub$Dose[is.na(med_select_id_sub$Dose)] = 0
med_select_id_sub$Strength_num[is.na(med_select_id_sub$Strength_num)] = 0

# Combine Med Name & med administration to fully characterize
continuous_app = c("Code/trauma/sedation medication", 
                   "Continuous Infusion: Per Instructions PRN",
                   "Continuous")
continuous_med = rep(0, nrow(med_select_id_sub))
continuous_med[med_select_id_sub$Frequency %in% continuous_app] = 1

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

status_med = rep(0, nrow(med_select_id_sub))
status_med[med_select_id_sub$Status %in% status_update[["Start"]]] = "Start"
status_med[med_select_id_sub$Status %in% status_update[["Stop"]]] = "Stop"
status_med[med_select_id_sub$Status %in% status_update[["Continue"]]] = "Continue"
status_med[med_select_id_sub$Status %in% status_update[["Changed"]]] = "Changed"

med_name_admin = paste0(med_select_id_sub$med_name_simple, continuous_med)

instance_num = rep(0, nrow(med_select_id_sub))
med_select_FINAL = cbind(med_select_id_sub, med_name_admin, status_med, continuous_med, instance_num)

# Verifying that the dose * strength = 0 for "Stopped"
print(unique(med_select_FINAL$Dose[med_select_FINAL$Status == "Stopped"]))
print(table(med_select_FINAL$status_med[med_select_FINAL$Dose == 0]))
dose_0_instance = med_select_FINAL[med_select_FINAL$Dose == 0, ]
ex_patient = med_select_FINAL[med_select_FINAL$key == 114500, ]

# Create a column with the "instance" number of that med to help determine if
#       the patient is just starting the medicine or if they've been using it
sub_med_list = NULL
first_id = med_select_FINAL$key[1]
for(i in 1:nrow(med_select_FINAL)) {
    if(med_select_FINAL$key[i] != first_id) {
        sub_med_list = NULL
        first_id = med_select_FINAL$key[i]
    }
    
    curr_med = med_select_FINAL$med_name_admin[i]
    sub_med_list = c(sub_med_list, curr_med)
    med_select_FINAL$instance_num[i] = sum(sub_med_list == curr_med)
    
}

hr_medications = med_select_FINAL[med_select_FINAL$hr != 0, ]
map_medications = med_select_FINAL[med_select_FINAL$map != 0, ]

# HR continuous
hr_cont = hr_medications[hr_medications$continuous_med == 1, ,drop=F]
hr_cont_names = unique(hr_cont$med_name_admin)
# HR discrete
hr_disc = hr_medications[hr_medications$continuous_med == 0, ,drop=F]
hr_disc_names = unique(hr_disc$med_name_admin)
# MAP continuous
map_cont = map_medications[map_medications$continuous_med == 1, ,drop=F]
map_cont_names = unique(map_cont$med_name_admin)
# MAP discrete
map_disc = map_medications[map_medications$continuous_med == 0, ,drop=F]
map_disc_names = unique(map_disc$med_name_admin)

EIDs = unique(data_format[,"EID"])
hr_cont_design = hr_disc_design = map_cont_design = map_disc_design = vector(mode = 'list', length = length(EIDs))

for(j in 1:length(EIDs)) { 
    if(EIDs[j] %in% med_select_FINAL$key) {
        print(paste0(j, ", ", EIDs[j]))
        
        sub_data = data_format[data_format[,'EID'] == EIDs[j], ]
        hr_cont_design[[j]] = matrix(0, ncol = length(unique(hr_cont$med_name_admin)),
                                     nrow = nrow(sub_data))
        colnames(hr_cont_design[[j]]) = hr_cont_names
        hr_disc_design[[j]] = matrix(0, ncol = length(unique(hr_disc$med_name_admin)),
                                     nrow = nrow(sub_data))
        colnames(hr_disc_design[[j]]) = hr_disc_names
        map_cont_design[[j]] = matrix(0, ncol = length(unique(map_cont$med_name_admin)),
                                      nrow = nrow(sub_data))
        colnames(map_cont_design[[j]]) = map_cont_names
        map_disc_design[[j]] = matrix(0, ncol = length(unique(map_disc$med_name_admin)),
                                      nrow = nrow(sub_data))
        colnames(map_disc_design[[j]]) = map_disc_names
        
        hr_cont_j = hr_cont[hr_cont$key == EIDs[j], ,drop=F]
        hr_disc_j = hr_disc[hr_disc$key == EIDs[j], ,drop=F]
        map_cont_j = map_cont[map_cont$key == EIDs[j], ,drop=F]
        map_disc_j = map_disc[map_disc$key == EIDs[j], ,drop=F]
        
        time_1 = time_2 = 0
        for(k in 1:nrow(sub_data)) {
            if(k == 1) {
                time_1 = 0; time_2 = sub_data[k,"time"]
            } else {
                time_1 = sub_data[k-1,"time"]; time_2 = sub_data[k,"time"]
            }
            
            # HR continuous --------------------------------------------------------
            hr_cont_meds = hr_cont_j[hr_cont_j$administered_dtm <= time_2 & 
                                         hr_cont_j$administered_dtm > time_1, , drop = F]
            if(nrow(hr_cont_meds) > 0) {
                for(jj in 1:nrow(hr_cont_meds)) {
                    # Need to look at instance_num and status_med
                    if(hr_cont_meds$status_med[jj] == "Stop") {
                        hr_cont_design[[j]][k:nrow(sub_data),hr_cont_meds$med_name_admin[jj]] = 0
                    } else if(hr_cont_meds$status_med[jj] == "Start") {
                        dosage = hr_cont_meds$Dose[jj] * hr_cont_meds$Strength_num[jj]
                        hr_cont_design[[j]][k:nrow(sub_data),hr_cont_meds$med_name_admin[jj]] = dosage
                    } else if(hr_cont_meds$status_med[jj] == "Continue") {
                        dosage = hr_cont_meds$Dose[jj] * hr_cont_meds$Strength_num[jj]
                        hr_cont_design[[j]][k:nrow(sub_data),hr_cont_meds$med_name_admin[jj]] = dosage
                        # This means we do not see the start of this medication
                        if(hr_cont_meds$instance_num[jj] == 1) {
                            hr_cont_design[[j]][1:k,hr_cont_meds$med_name_admin[jj]] = dosage
                        }
                    } else {
                        dosage = hr_cont_meds$Dose[jj] * hr_cont_meds$Strength_num[jj]
                        hr_cont_design[[j]][k:nrow(sub_data),hr_cont_meds$med_name_admin[jj]] = dosage
                        # This means we do not see the start of this medication
                        # (ASSUME the dose is the same beforehand)
                        if(hr_cont_meds$instance_num[jj] == 1) {
                            hr_cont_design[[j]][1:k,hr_cont_meds$med_name_admin[jj]] = dosage
                        }
                    }
                }
            }
            
            # MAP continuous -------------------------------------------------------
            map_cont_meds = map_cont_j[map_cont_j$administered_dtm <= time_2 & 
                                           map_cont_j$administered_dtm > time_1, , drop = F]
            if(nrow(map_cont_meds) > 0) {
                for(jj in 1:nrow(map_cont_meds)) {
                    # Need to look at instance_num and status_med
                    if(map_cont_meds$status_med[jj] == "Stop") {
                        map_cont_design[[j]][k:nrow(sub_data),map_cont_meds$med_name_admin[jj]] = 0
                    } else if(map_cont_meds$status_med[jj] == "Start") {
                        dosage = map_cont_meds$Dose[jj] * map_cont_meds$Strength_num[jj]
                        map_cont_design[[j]][k:nrow(sub_data),map_cont_meds$med_name_admin[jj]] = dosage
                    } else if(map_cont_meds$status_med[jj] == "Continue") {
                        dosage = map_cont_meds$Dose[jj] * map_cont_meds$Strength_num[jj]
                        map_cont_design[[j]][k:nrow(sub_data),map_cont_meds$med_name_admin[jj]] = dosage
                        # This means we do not see the start of this medication
                        if(map_cont_meds$instance_num[jj] == 1) {
                            map_cont_design[[j]][1:k,map_cont_meds$med_name_admin[jj]] = dosage
                        }
                    } else {
                        dosage = map_cont_meds$Dose[jj] * map_cont_meds$Strength_num[jj]
                        map_cont_design[[j]][k:nrow(sub_data),map_cont_meds$med_name_admin[jj]] = dosage
                        # This means we do not see the start of this medication
                        # (ASSUME the dose is the same beforehand)
                        if(map_cont_meds$instance_num[jj] == 1) {
                            map_cont_design[[j]][1:k,map_cont_meds$med_name_admin[jj]] = dosage
                        }
                    }
                }
            }
            
            # HR discrete ----------------------------------------------------------
            hr_disc_meds = hr_disc_j[hr_disc_j$administered_dtm <= time_2 & 
                                         hr_disc_j$administered_dtm > time_1, , drop = F]
            if(nrow(hr_disc_meds) > 0) {
                for(jj in 1:nrow(hr_disc_meds)) {
                    # Need to look at instance_num and status_med
                    if(hr_disc_meds$status_med[jj] == "Stop") {
                        hr_disc_design[[j]][k:nrow(sub_data),hr_disc_meds$med_name_admin[jj]] = 0
                    } else if(hr_disc_meds$status_med[jj] == "Start") {
                        dosage = hr_disc_meds$Dose[jj] * hr_disc_meds$Strength_num[jj]
                        total_time = med_key$onset[med_key$id == hr_disc_meds$med_name_simple[jj]] + 
                            med_key$offset[med_key$id == hr_disc_meds$med_name_simple[jj]]
                        
                        time_ind = which(sub_data[,"time"] <= hr_disc_meds$administered_dtm[jj] + total_time)
                        if(length(time_ind) > 0) {
                            max_ind = max(time_ind)
                        } else {
                            max_ind = k
                        }
                        # Medications are additive
                        hr_disc_design[[j]][k:max_ind,hr_disc_meds$med_name_admin[jj]] = 
                            hr_disc_design[[j]][k:max_ind,hr_disc_meds$med_name_admin[jj]] + dosage
                    } else if(hr_disc_meds$status_med[jj] == "Continue") {
                        dosage = hr_disc_meds$Dose[jj] * hr_disc_meds$Strength_num[jj]
                        total_time = med_key$onset[med_key$id == hr_disc_meds$med_name_simple[jj]] + 
                            med_key$offset[med_key$id == hr_disc_meds$med_name_simple[jj]]
                        time_ind = which(sub_data[,"time"] <= hr_disc_meds$administered_dtm[jj] + total_time)
                        if(length(time_ind) > 0) {
                            max_ind = max(time_ind)
                        } else {
                            max_ind = k
                        }
                        # Medications are additive
                        hr_disc_design[[j]][k:max_ind,hr_disc_meds$med_name_admin[jj]] = 
                            hr_disc_design[[j]][k:max_ind,hr_disc_meds$med_name_admin[jj]] + dosage
                    } else {
                        dosage = hr_disc_meds$Dose[jj] * hr_disc_meds$Strength_num[jj]
                        total_time = med_key$onset[med_key$id == hr_disc_meds$med_name_simple[jj]] + 
                            med_key$offset[med_key$id == hr_disc_meds$med_name_simple[jj]]
                        time_ind = which(sub_data[,"time"] <= hr_disc_meds$administered_dtm[jj] + total_time)
                        if(length(time_ind) > 0) {
                            max_ind = max(time_ind)
                        } else {
                            max_ind = k
                        }
                        # Medications are additive
                        hr_disc_design[[j]][k:max_ind,hr_disc_meds$med_name_admin[jj]] = 
                            hr_disc_design[[j]][k:max_ind,hr_disc_meds$med_name_admin[jj]] + dosage
                    }
                }
            }
            
            # MAP discrete ---------------------------------------------------------
            map_disc_meds = map_disc_j[map_disc_j$administered_dtm <= time_2 & 
                                           map_disc_j$administered_dtm > time_1, , drop = F]
            if(nrow(map_disc_meds) > 0) {
                for(jj in 1:nrow(map_disc_meds)) {
                    # Need to look at instance_num and status_med
                    if(map_disc_meds$status_med[jj] == "Stop") {
                        map_disc_design[[j]][k:nrow(sub_data),map_disc_meds$med_name_admin[jj]] = 0
                    } else if(map_disc_meds$status_med[jj] == "Start") {
                        dosage = map_disc_meds$Dose[jj] * map_disc_meds$Strength_num[jj]
                        total_time = med_key$onset[med_key$id == map_disc_meds$med_name_simple[jj]] + 
                            med_key$offset[med_key$id == map_disc_meds$med_name_simple[jj]]
                        time_ind = which(sub_data[,"time"] <= map_disc_meds$administered_dtm[jj] + total_time)
                        if(length(time_ind) > 0) {
                            max_ind = max(time_ind)
                        } else {
                            max_ind = k
                        }
                        # Medications are additive
                        map_disc_design[[j]][k:max_ind,map_disc_meds$med_name_admin[jj]] = 
                            map_disc_design[[j]][k:max_ind,map_disc_meds$med_name_admin[jj]] + dosage
                    } else if(map_disc_meds$status_med[jj] == "Continue") {
                        dosage = map_disc_meds$Dose[jj] * map_disc_meds$Strength_num[jj]
                        total_time = med_key$onset[med_key$id == map_disc_meds$med_name_simple[jj]] + 
                            med_key$offset[med_key$id == map_disc_meds$med_name_simple[jj]]
                        time_ind = which(sub_data[,"time"] <= map_disc_meds$administered_dtm[jj] + total_time)
                        if(length(time_ind) > 0) {
                            max_ind = max(time_ind)
                        } else {
                            max_ind = k
                        }
                        # Medications are additive
                        map_disc_design[[j]][k:max_ind,map_disc_meds$med_name_admin[jj]] = 
                            map_disc_design[[j]][k:max_ind,map_disc_meds$med_name_admin[jj]] + dosage
                    } else {
                        dosage = map_disc_meds$Dose[jj] * map_disc_meds$Strength_num[jj]
                        total_time = med_key$onset[med_key$id == map_disc_meds$med_name_simple[jj]] + 
                            med_key$offset[med_key$id == map_disc_meds$med_name_simple[jj]]
                        time_ind = which(sub_data[,"time"] <= map_disc_meds$administered_dtm[jj] + total_time)
                        if(length(time_ind) > 0) {
                            max_ind = max(time_ind)
                        } else {
                            max_ind = k
                        }
                        # Medications are additive
                        map_disc_design[[j]][k:max_ind,map_disc_meds$med_name_admin[jj]] = 
                            map_disc_design[[j]][k:max_ind,map_disc_meds$med_name_admin[jj]] + dosage
                    }
                }
            }
            
        }  
    } else {
        print(paste0(j, ", ", EIDs[j], " no meds"))
    }
    
}
