# Load the medication data
jw15 = read.csv("Data/_raw_data_new/jw15b.csv")
jw16 = read.csv("Data/_raw_data_new/jw16b.csv")
jw17 = read.csv("Data/_raw_data_new/jw17b.csv")
jw18 = read.csv("Data/_raw_data_new/jw18b.csv")
jw19 = read.csv("Data/_raw_data_new/jw19b.csv")

data_num = 1

# Load the existing patient information
load('Data_updates/data_format.rda')
select_id = unique(data_format[,"EID"])
med_select_id = jw15[jw15$key %in% select_id, ]
med_select_id = rbind(med_select_id, jw16[jw16$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw17[jw17$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw18[jw18$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw19[jw19$key %in% select_id, ])
med_select_id$administered_dtm = med_select_id$administered_dtm / 60
med_select_id = med_select_id[!is.na(med_select_id[,1]),]
rownames(med_select_id) = NULL
med_select_id = med_select_id[med_select_id$Med_Name_Desc != "", ]

med_key = read.csv('Med_chart.csv', na.strings = "")
colnames(med_key) = c('id', 'hr', 'map', 'onset', 'offset', 'time_check', 'X1')
med_key$id = paste0(" ", med_key$id)
med_key$id[med_key$id == " Metoprolol IV"] = " METOPROLOL"
library(stringr)
med_name_simple = paste0(" ", med_select_id$Med_Name_Desc)

unique_meds = unique(med_name_simple)
med_list = matrix(0, nrow = length(unique_meds), ncol = nrow(med_key))
for(i in 1:nrow(med_key)) {
    ind = grep(med_key$id[i], unique_meds)
    if(length(ind) > 0) {
        med_list[ind, i] = med_list[ind, i] +1
    }
}
rownames(med_list) = unique_meds
colnames(med_list) = med_key$id
check_vec = apply(med_list, 1, sum)
watch_out = which(check_vec != 1)
names(watch_out) = NULL
print(length(watch_out))

med_name_simple_mat = matrix(nrow = length(unique_meds), ncol = 2)
med_name_simple_mat[,1] = unique_meds
for(i in 1:nrow(med_name_simple_mat)) {
    if(!(i %in% watch_out)) {
        i_name = colnames(med_list)[which(med_list[i, ] == 1)]
        med_name_simple_mat[i,2] = i_name
    }
}

print(paste0("Number of mismatches before: ", sum(is.na(med_name_simple_mat))))
watch_out_meds = data.frame("ind" = watch_out,
                            "med" = rep(NA, length(watch_out)))
for(i in 1:length(watch_out)) {
    full_term = unique_meds[watch_out[i]]
    test_phrase = med_key$id[med_list[watch_out[i],] != 0]
    test_phrase_pos = NULL
    for(j in test_phrase) {
        test_phrase_pos = c(test_phrase_pos, stringi::stri_locate_first_fixed(full_term, j)[, 1])
    }
    
    if(length(which(test_phrase_pos == min(test_phrase_pos))) > 1) {
        n_char_phrase = NULL
        for(j in test_phrase) {
            n_char_phrase = c(n_char_phrase, nchar(j))
        }
        watch_out_meds[i, "med"] = test_phrase[which.max(n_char_phrase)]   
        med_name_simple_mat[watch_out[i], 2] = test_phrase[which.max(n_char_phrase)]   
    } else {
        watch_out_meds[i, "med"] = test_phrase[which.min(test_phrase_pos)]   
        med_name_simple_mat[watch_out[i], 2] = test_phrase[which.min(test_phrase_pos)]   
    }
}
print("Original")
print(names(check_vec[watch_out]))
print("Corrected")
print(watch_out_meds)
print(paste0("Number of mismatches after: ", sum(is.na(med_name_simple_mat))))

hr_map = matrix(0, ncol = 2, nrow = length(med_name_simple))
colnames(hr_map) = c("hr", "map")
for(i in 1:nrow(med_key)) {
    simple_name_ind = med_name_simple_mat[med_name_simple_mat[,2] == med_key$id[i], 1]
    ind = which(med_name_simple %in% simple_name_ind)
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

# Combine Med Name & med administration to fully characterize
# Do we count administration at a frequency of "every 15 min" or faster as continuous?
# temp = med_select_id_sub[med_select_id_sub$Frequency %in% c("Every 5 min PRN", "Every 10 min PRN",
#                                                             "Every 2 min PRN", "Every 15 min"), ]

continuous_app = c("Continuous Infusion: Per Instructions PRN",
                   "Continuous")
names_of_meds_cont = unique(med_select_id_sub$Med_Name_Desc[med_select_id_sub$Frequency %in% continuous_app])
continuous_med = rep(0, nrow(med_select_id_sub))
continuous_med[med_select_id_sub$Frequency %in% continuous_app] = 1
continuous_med[med_select_id_sub$Med_Name_Desc %in% names_of_meds_cont] = 1

med_select_id_sub = cbind(med_select_id_sub, continuous_med)

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
print("Current Med Status Values:")
print(unique(status_med))

# Looking through if things are too simplified ---------------------------------
meds_to_check = NULL
total_simple_meds = sort(unique(med_select_id_sub$med_name_simple))
for(i in 1:length(total_simple_meds)) {
    med_d = med_select_id_sub[med_select_id_sub$med_name_simple == total_simple_meds[i], ,drop = F] 
    if(sum(med_d$continuous_med == 1) != nrow(med_d) & sum(med_d$continuous_med == 1) != 0) {
        meds_to_check = c(meds_to_check, total_simple_meds[i])   
    }
}

# List: meds_to_check
print("medications that have both continuous and discrete forms of admin")
print(meds_to_check)

prev_meds_to_check = c(" ALBUMIN", " AMIODARONE ", " CALCIUM CHLORIDE", 
                       " CLEVIDIPINE ", " DEXMEDETOMIDINE ", " DILTIAZEM ", 
                       " EPINEPHRINE ", " ESMOLOL ", " KETAMINE ", 
                       " LABETALOL ", " MILRINONE ", " NITROGLYCERIN ", 
                       " PHENYLEPHRINE", " PROPOFOL ", " SODIUM BICARBONATE ", 
                       " VASOPRESSIN ")

print(paste0("The previous med_check list contains ", sum(meds_to_check %in% prev_meds_to_check), 
      " of the total ", length(meds_to_check), " current meds to check"))

med_d = med_select_id_sub[med_select_id_sub$med_name_simple == prev_meds_to_check[11], ,drop = F]
unique(med_d$Med_Name_Desc[med_d$continuous_med == 1])
unique(med_d$Med_Name_Desc[med_d$continuous_med == 0])

med_name_simple_new = paste0(med_select_id_sub$med_name_simple, continuous_med)

# med_name_admin = paste0(med_select_id_sub$Med_Name_Desc, continuous_med)
med_name_admin = med_name_simple_new
unique_med_names = sort(unique(med_name_admin)); 
print(paste0("Num. meds broken up based on continuous vs. discrete admin: ", 
             length(unique_med_names)))
# ------------------------------------------------------------------------------

instance_num = rep(0, nrow(med_select_id_sub))

# FINAL MEDICATION DF *********************************************************
med_select_update = cbind(med_select_id_sub, med_name_admin, status_med, instance_num)
med_select_FINAL = matrix(nrow = 1, ncol = ncol(med_select_update))
colnames(med_select_FINAL) = colnames(med_select_update)
for(i in unique(med_select_update$key)){
    max_time = max(data_format[data_format[,"EID"] == i, "time"])
    sub_group = med_select_update[med_select_update$key == i, ,drop = F]
    ind_grab = which(sub_group$administered_dtm <= max_time)
    if(length(ind_grab) > 0) {
        med_select_FINAL = rbind(med_select_FINAL, sub_group[ind_grab, ,drop = F])   
    }
}
med_select_FINAL = med_select_FINAL[-1, ,drop=F]; rownames(med_select_FINAL) = NULL
# med_select_FINAL = med_select_FINAL[med_select_FINAL$Dose != 0, ]
# med_select_FINAL = med_select_FINAL[!is.na(med_select_FINAL$Dose), ]
# med_select_FINAL$continuous_med[med_select_FINAL$med_name_admin == "ALBUMIN_0"] = 0
# *****************************************************************************

# First, check all time instances when the medication is stopped and set
# the dose and strength to 0.
med_select_FINAL[med_select_FINAL$status_med == "Stop", "Dose"] = 0
med_select_FINAL[med_select_FINAL$status_med == "Stop", "Strength"] = ""

# Second, when is Dose NA
na_dose = med_select_FINAL[is.na(med_select_FINAL$Dose), ]

# making the strength numeric
all_med_types = unique(med_select_FINAL$med_name_admin)
strength_units = dose_units = vector(mode = 'list', length = length(all_med_types))
for(i in 1:length(strength_units)) {
    strength_units[[i]] = unique(med_select_FINAL$Strength[med_select_FINAL$med_name_admin == all_med_types[i]])
    dose_units[[i]] = unique(med_select_FINAL$Dose_Units[med_select_FINAL$med_name_admin == all_med_types[i]])
}
names(strength_units) = names(dose_units) = all_med_types

med_names_alphabet = sort(all_med_types)

# ******************************************************************************
# Manually check the strength units and dose units to see which has multiple ***
# ******************************************************************************
# 1. " ADENOSINE 0" ------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ADENOSINE 0"] = 1
med_select_FINAL$Dose[med_select_FINAL$med_name_admin == " ADENOSINE 0" & 
                          is.na(med_select_FINAL$Dose)] = 0
# 2. " ALBUMIN1" ---------------------------------------------------------------
# mL -> g (divide by 10)
alb_ind = which(med_select_FINAL$med_name_admin == " ALBUMIN1" & 
                    med_select_FINAL$Dose_Units == "mL")
med_select_FINAL[alb_ind, "Dose"] = med_select_FINAL[alb_ind, "Dose"] / 10
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ALBUMIN1" & 
                              med_select_FINAL$Strength == "25 %"] = 0.25
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ALBUMIN1" & 
                              med_select_FINAL$Strength == "5 %"] = 0.05
alb_ind_miss = which(med_select_FINAL$med_name_admin == " ALBUMIN1" & 
                         med_select_FINAL$Strength == "")
med_select_FINAL[alb_ind_miss, "Dose"] = med_select_FINAL[alb_ind_miss, "Strength"] = 0

# 3. " AMIODARONE 0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " AMIODARONE 0")
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], " ")[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " AMIODARONE 0"] = 1

# 4. " AMIODARONE 1" -----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " AMIODARONE 1"] = 1

# 5. " AMLODIPINE 0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " AMLODIPINE 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], " ")[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " AMLODIPINE 0"] = 1

# 6. " ANGIOTENSIN 1" ----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ANGIOTENSIN 1"] = 1

# 7. " ATENOLOL 0" -------------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " ATENOLOL 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[- ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ATENOLOL 0"] = 1

# 8. " ATROPINE 0" -------------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " ATROPINE 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ATROPINE 0"] = 1

# 9. " BISOPROLOL0" ------------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " BISOPROLOL0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " BISOPROLOL0"] = 1

# 10. " CALCIUM CHLORIDE0" -----------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " CALCIUM CHLORIDE0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

cc_ind = which(med_select_FINAL$med_name_admin == " CALCIUM CHLORIDE0" &  
                   (med_select_FINAL$Dose_Units %in% c("mg", "mg/hr")))
med_select_FINAL$Dose[cc_ind] = med_select_FINAL$Dose[cc_ind] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CALCIUM CHLORIDE0"] = 1

# 11. " CALCIUM CHLORIDE1" -----------------------------------------------------
cc_ind = which(med_select_FINAL$med_name_admin == " CALCIUM CHLORIDE1" &  
                   (med_select_FINAL$Dose_Units %in% c("g/day")))
med_select_FINAL$Dose[cc_ind] = med_select_FINAL$Dose[cc_ind] * (1000 / 24)
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CALCIUM CHLORIDE1"] = 1

# 12. " CALCIUM GLUCONATE0" ----------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

cal_g_ind = which(med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0" & 
                      med_select_FINAL$Dose_Units == "mg")
med_select_FINAL[cal_g_ind, "Dose"] = med_select_FINAL[cal_g_ind, "Dose"] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0" & 
                              med_select_FINAL$Strength == "10%"] = 5
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0" &
                              med_select_FINAL$Strength != 5] = 1

# 13. " CANDESARTAN 0" ---------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " CANDESARTAN 0", ]

# 14. " CAPTOPRIL 0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " CAPTOPRIL 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CAPTOPRIL 0"] = 1

# 15. " CARVEDILOL0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " CARVEDILOL0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Med_Name_Desc[na_dose[d]], '[ ]+')[[1]][2])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CARVEDILOL0"] = 1

# 16. " CLEVIDIPINE 0" ---------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CLEVIDIPINE 0"] = 1

# 17. " CLEVIDIPINE 1" -------------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " CLEVIDIPINE 1")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CLEVIDIPINE 1"] = 1

# 18. " CLONIDINE 0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " CLONIDINE 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " CLONIDINE 0", ]

# 19. " DEXMEDETOMIDINE 0" -----------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " DEXMEDETOMIDINE 0", ]

# 20. " DEXMEDETOMIDINE 1" -----------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " DEXMEDETOMIDINE 1" ] = 1

# 21. " DIGOXIN 0" -------------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " DIGOXIN 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]
d_0 = which(med_select_FINAL$med_name_admin == " DIGOXIN 0" & med_select_FINAL$Dose_Units == "mg")
med_select_FINAL = med_select_FINAL[-d_0, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " DIGOXIN 0" ] = 1

# 22. " DILTIAZEM 0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " DILTIAZEM 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " DILTIAZEM 0"] = 1

# 23. " DILTIAZEM 1" -----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " DILTIAZEM 1"] = 1

# 24. " DOBUTAMINE 1" ----------------------------------------------------------
dob_1 = which(med_select_FINAL$med_name_admin == " DOBUTAMINE 1" & med_select_FINAL$Dose_Units == "mg")
med_select_FINAL = med_select_FINAL[-dob_1, ]
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " DOBUTAMINE 1"] = 1

# 25. " DOPAMINE 1" ------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " DOPAMINE 1"] = 1

# 26. " ENALAPRIL0" ------------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " ENALAPRIL0", ]

# 27. " ENALAPRILAT 0" ---------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " ENALAPRILAT 0", ]

# 28. " EPHEDRINE0" ------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " EPHEDRINE0"] = 1

# 29. " EPINEPHRINE 0" ---------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " EPINEPHRINE 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

ep_ind = which(med_select_FINAL$med_name_admin == " EPINEPHRINE 0" &  med_select_FINAL$Dose_Units %in% c("mcg","mcg/kg/min"))
med_select_FINAL$Dose[ep_ind] = med_select_FINAL$Dose[ep_ind] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " EPINEPHRINE 0"] = 1

# 30. " EPINEPHRINE 1" ---------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " EPINEPHRINE 1")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " EPINEPHRINE 1"] = 1

# 31. " ESMOLOL 0" -------------------------------------------------------------  
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " ESMOLOL 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ESMOLOL 0"] = 1

# 32. " ESMOLOL 1" ------------------------------------------------------------- 
t_strength = which((med_select_FINAL$med_name_admin == " ESMOLOL 1") & (med_select_FINAL$Strength == "20 mg/mL"))
med_select_FINAL$Strength[t_strength] = 2
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ESMOLOL 1" &
                              med_select_FINAL$Strength != "20 mg/mL"] = 1

# 33. " IRBESARTAN 0" ----------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " IRBESARTAN 0", ]

# 34. " ISOPROTERENOL 1" -------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ISOPROTERENOL 1"] = 1

# 35. " ISOSORBIDE DINITRATE0" -------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " ISOSORBIDE DINITRATE0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ISOSORBIDE DINITRATE0"] = 1

# 36. " ISOSORBIDE MONONITRATE0" -----------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " ISOSORBIDE MONONITRATE0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " ISOSORBIDE MONONITRATE0"] = 1

# 37. " KETAMINE 0" ------------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " KETAMINE 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " KETAMINE 0"] = 1

# 38. " KETAMINE 1" ------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " KETAMINE 1"] = 1

# 39. " LABETALOL 0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " LABETALOL 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " LABETALOL 0"] = 1

# 40. " LABETALOL 1" -----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " LABETALOL 1"] = 1

# 41. " LISINOPRIL 0" ----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " LISINOPRIL 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Med_Name_Desc[na_dose[d]], '[ ]+')[[1]][2])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " LISINOPRIL 0"] = 1

# 42. " LOSARTAN 0"  -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " LOSARTAN 0" )
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Med_Name_Desc[na_dose[d]], '[ ]+')[[1]][2])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " LOSARTAN 0" ] = 1

# 43. " LOSARTAN-HYDROCHLOROTHIAZIDE 0"-----------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " LOSARTAN-HYDROCHLOROTHIAZIDE 0", ]

# 44. " METOPROLOL SUCCINATE0"--------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " METOPROLOL SUCCINATE0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " METOPROLOL SUCCINATE0"] = 1

# 45. " METOPROLOL TARTRATE0" --------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " METOPROLOL TARTRATE0" )
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Med_Name_Desc[na_dose[d]], '[ ]+')[[1]][3])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " METOPROLOL TARTRATE0" ] = 1

# 46. " METOPROLOL0" -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " METOPROLOL0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " METOPROLOL0"] = 1

# 47. " MILRINONE 0" -----------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " MILRINONE 0", ]

# 48. " MILRINONE 1" -----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " MILRINONE 1"] = 1

# 49. " NADOLOL 0" -------------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " NADOLOL 0", ]

# 50. " NICARDIPINE 1"----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NICARDIPINE 1"] = 1

# 51. " NIFEDIPINE0" -----------------------------------------------------------
app_nif = which(med_select_FINAL$med_name_admin == " NIFEDIPINE0" & 
                    med_select_FINAL$Dose_Units == "application")
med_select_FINAL = med_select_FINAL[-app_nif, ]
app_nif = which(med_select_FINAL$med_name_admin == " NIFEDIPINE0" & 
                    med_select_FINAL$Strength == "0.1 %")
med_select_FINAL = med_select_FINAL[-app_nif, ]

na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " NIFEDIPINE0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NIFEDIPINE0"] = 1

# 52. " NIMODIPINE 0"-----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NIMODIPINE 0"] = 1

# 53. " NITROGLYCERIN 0"--------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " NITROGLYCERIN 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

nitro_inch = which(med_select_FINAL$med_name_admin == " NITROGLYCERIN 0" & 
                       med_select_FINAL$Dose_Units == "inch")
med_select_FINAL = med_select_FINAL[-nitro_inch, ]

med_select_FINAL$Dose[med_select_FINAL$med_name_admin == " NITROGLYCERIN 0" & 
                          med_select_FINAL$Dose_Units == "mcg"] = med_select_FINAL$Dose[med_select_FINAL$med_name_admin == " NITROGLYCERIN 0" & 
                                                                                            med_select_FINAL$Dose_Units == "mcg"] / 1000

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NITROGLYCERIN 0"] = 1

# 54. " NITROGLYCERIN 1" -------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " NITROGLYCERIN 1")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

mcg_ind = which(med_select_FINAL$med_name_admin == " NITROGLYCERIN 1" & 
                    med_select_FINAL$Dose_Units %in% c("mcg", "mcg/kg/min"))

med_select_FINAL$Dose[mcg_ind] = med_select_FINAL$Dose[mcg_ind] / 1000

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NITROGLYCERIN 1"] = 1

# 55. " NITROPRUSSIDE 1" -------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NITROPRUSSIDE 1"] = 1

# 56. " NOREPINEPHRINE1" -------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " NOREPINEPHRINE1")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

new_bag_nor = which(med_select_FINAL$med_name_admin == " NOREPINEPHRINE1" & med_select_FINAL$Dose > 16)
med_select_FINAL = med_select_FINAL[-new_bag_nor, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NOREPINEPHRINE1"] = 1

# 57. " OLMESARTAN 0"  ---------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " OLMESARTAN 0", ]

# 58. " PHENYLEPHRINE0" --------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " PHENYLEPHRINE0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

p_ind = which(med_select_FINAL$med_name_admin == " PHENYLEPHRINE0" &  med_select_FINAL$Dose_Units == "mcg")
med_select_FINAL$Dose[p_ind] = med_select_FINAL$Dose[p_ind] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " PHENYLEPHRINE0"] = 1

# 59. " PHENYLEPHRINE1"---------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " PHENYLEPHRINE1")
med_select_FINAL = med_select_FINAL[-na_dose, ]
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " PHENYLEPHRINE1"] = 1

# 60. " PROPOFOL 0"    ---------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " PROPOFOL 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

prop_mcg = which(med_select_FINAL$med_name_admin == " PROPOFOL 0" &  med_select_FINAL$Dose_Units == "mcg/kg/min")
med_select_FINAL$Dose[prop_mcg] = med_select_FINAL$Dose[prop_mcg] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " PROPOFOL 0"] = 1

# 61. " PROPOFOL 1"  -----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " PROPOFOL 1")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " PROPOFOL 1"] = 1

# 62. " SILDENAFIL 0"   --------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " SILDENAFIL 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " SILDENAFIL 0"] = 1

# 63. " SODIUM BICARBONATE 0" --------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " SODIUM BICARBONATE 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " SODIUM BICARBONATE 0"] = 1

# 64. " SODIUM BICARBONATE 1" --------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " SODIUM BICARBONATE 1")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " SODIUM BICARBONATE 1"] = 1

# 65. " SOTALOL 0"   -----------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " SOTALOL 0", ]

# 66. " SPIRONOLACTONE 0"   ----------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " SPIRONOLACTONE 0")
na_med = med_select_FINAL[na_dose, ]
dose_nums = rep(0, length(na_dose))
for(d in 1:length(dose_nums)) {
    dose_nums[d] = as.numeric(strsplit(med_select_FINAL$Strength[na_dose[d]], '[ ]+')[[1]][1])
}
med_select_FINAL$Dose[na_dose] = dose_nums

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " SPIRONOLACTONE 0"] = 1

# 67. " VALSARTAN 0"     -------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " VALSARTAN 0", ]

# 68. " VASOPRESSIN 0"   -------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " VASOPRESSIN 0"] = 1

# 69. " VASOPRESSIN 1"   -------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " VASOPRESSIN 1")
na_keys = med_select_FINAL$key[na_dose]
temp = med_select_FINAL[med_select_FINAL$key %in% na_keys & med_select_FINAL$med_name_admin == " VASOPRESSIN 1", ]

med_select_FINAL$Dose[na_dose] = 0.04

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " VASOPRESSIN 1"] = 1

# 70. " VERAPAMIL 0"  ----------------------------------------------------------
na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == " VERAPAMIL 0")
med_select_FINAL = med_select_FINAL[-na_dose, ]

med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " VERAPAMIL 0"] = 1

# # ==============================================================================
# m_id = 70
# print(med_names_alphabet[m_id])
# specific_med = med_select_FINAL[med_select_FINAL$med_name_admin == med_names_alphabet[m_id], ]
# apply(specific_med[,c("Dose", "Dose_Units", "Strength")], 2, table)
# na_dose = which(is.na(med_select_FINAL$Dose) & med_select_FINAL$med_name_admin == med_names_alphabet[m_id])
# na_med = med_select_FINAL[na_dose, ]
# # ==============================================================================

# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


Strength_num = as.numeric(med_select_FINAL$Strength)
med_select_FINAL = cbind(med_select_FINAL, Strength_num)
med_select_FINAL$Dose[is.na(med_select_FINAL$Dose)] = 0
med_select_FINAL$Strength_num[is.na(med_select_FINAL$Strength_num)] = 0

save(med_select_FINAL, file = paste0("Data_updates/med_select_FINAL", data_num, ".rda"))

# SCALING DOSE! ************************************************************
med_dose_unique = unique(med_select_FINAL$med_name_admin)
med_dose_scale_factor = matrix(nrow=length(med_dose_unique), ncol = 2)
med_dose_scale_factor[,1] = med_dose_unique
colnames(med_dose_scale_factor) = c('name', 'mean')
for(i in 1:length(med_dose_unique)){
    med_i = which(med_select_FINAL$med_name_admin == med_dose_unique[i])
    mean_i = mean(med_select_FINAL$Dose[med_i])
    med_dose_scale_factor[i, 2] = mean_i
    med_select_FINAL$Dose[med_i] = med_select_FINAL$Dose[med_i] / mean_i
}
save(med_dose_scale_factor, file = paste0('Data_updates/med_dose_scale_factor', data_num, '.rda'))
# ******************************************************************************

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

Dn_omega_names = c(hr_cont_names, hr_disc_names, map_cont_names, map_disc_names)
save(Dn_omega_names, file = paste0("Data_updates/Dn_omega_names", data_num, ".rda"))
hr_map_names = c(rep("hr_cont", length(hr_cont_names)),
                 rep("hr_disc", length(hr_disc_names)),
                 rep("map_cont", length(map_cont_names)),
                 rep("map_disc", length(map_disc_names)))
save(hr_map_names, file = paste0('Data_updates/hr_map_names', data_num, '.rda'))

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

Dn_omega = vector(mode = 'list', length = length(EIDs))
hr_cont_num = length(hr_cont_names)
hr_disc_num = length(hr_disc_names)
map_cont_num = length(map_cont_names)
map_disc_num = length(map_disc_names)
total_cols_hr = hr_cont_num + hr_disc_num
total_cols_map = map_cont_num + map_disc_num
for(i in 1:length(Dn_omega)) {
    Dn_omega[[i]] = vector(mode = 'list', length = sum(data_format[,"EID"] == EIDs[i]))
    for(j in 1:length(Dn_omega[[i]])) {
        med_mat = matrix(0, nrow = 4, ncol = total_cols_hr + total_cols_map)
        hr_info = map_info = NULL
        
        if(is.null(hr_cont_design[[i]])) {
            hr_info = rep(0, hr_cont_num)    
        } else {
            hr_info = hr_cont_design[[i]][j,]
        }
        
        if(is.null(hr_disc_design[[i]])) {
            hr_info = c(hr_info, rep(0, hr_disc_num))
        } else {
            hr_info = c(hr_info, hr_disc_design[[i]][j,])
        }
        names(hr_info) = NULL
        
        if(is.null(map_cont_design[[i]])) {
            map_info = rep(0, map_cont_num)    
        } else {
            map_info = map_cont_design[[i]][j,]
        }
        
        if(is.null(map_disc_design[[i]])) {
            map_info = c(map_info, rep(0, map_disc_num))
        } else {
            map_info = c(map_info, map_disc_design[[i]][j,])
        }
        names(map_info) = NULL
        
        med_mat[2, 1:total_cols_hr] = hr_info
        med_mat[3, (total_cols_hr+1):(total_cols_hr+total_cols_map)] = map_info
        
        Dn_omega[[i]][[j]] = med_mat
    }
}

save(Dn_omega, file = paste0('Data_updates/Dn_omega', data_num, '.rda'))

# Understanding what the mean of Dn_omega should be
upp_down_omega = matrix(nrow = length(Dn_omega_names), ncol = 2)
upp_down_omega[,1] = Dn_omega_names
ind = 1
for(i in 1:length(hr_cont_names)) {
    if(upp_down_omega[ind,1] != hr_cont_names[i]) {
        print(hr_cont_names[i])
    } else{
        ef = unique(med_select_FINAL$hr[med_select_FINAL$med_name_admin == hr_cont_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
for(i in 1:length(hr_disc_names)) {
    if(upp_down_omega[ind,1] != hr_disc_names[i]) {
        print(hr_disc_names[i])
    } else{
        ef = unique(med_select_FINAL$hr[med_select_FINAL$med_name_admin == hr_disc_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
for(i in 1:length(map_cont_names)) {
    if(upp_down_omega[ind,1] != map_cont_names[i]) {
        print(map_cont_names[i])
    } else{
        ef = unique(med_select_FINAL$map[med_select_FINAL$med_name_admin == map_cont_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
for(i in 1:length(map_disc_names)) {
    if(upp_down_omega[ind,1] != map_disc_names[i]) {
        print(map_disc_names[i])
    } else{
        ef = unique(med_select_FINAL$map[med_select_FINAL$med_name_admin == map_disc_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}

mean_dn_omega = as.numeric(c(upp_down_omega[,2]))
# mean_dn_omega = 4 * mean_dn_omega
print(c(mean_dn_omega))
