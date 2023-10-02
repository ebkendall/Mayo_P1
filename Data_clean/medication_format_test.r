jw15 = read.csv("Data/_raw_data_new/jw15b.csv")
jw16 = read.csv("Data/_raw_data_new/jw16b.csv")
jw17 = read.csv("Data/_raw_data_new/jw17b.csv")
jw18 = read.csv("Data/_raw_data_new/jw18b.csv")
jw19 = read.csv("Data/_raw_data_new/jw19b.csv")

# load('Data/data_format_new.rda')
# pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
#             417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
#             732525, 758025, 763050, 843000)
# data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
load('Data/data_format_new2.rda')
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
print(watch_out)

med_name_simple_mat = matrix(nrow = length(unique_meds), ncol = 2)
med_name_simple_mat[,1] = unique_meds
for(i in 1:nrow(med_name_simple_mat)) {
    if(!(i %in% watch_out)) {
        i_name = colnames(med_list)[which(med_list[i, ] == 1)]
        med_name_simple_mat[i,2] = i_name
    }
}
med_name_simple_mat[1,2] = " METOPROLOL TARTRATE" 
med_name_simple_mat[5,2] = " METOPROLOL SUCCINATE"
med_name_simple_mat[8,2] = " EPINEPHRINE "
med_name_simple_mat[20,2] = " PHENYLEPHRINE"
med_name_simple_mat[22,2] = " METOPROLOL TARTRATE"
med_name_simple_mat[25,2] = " CALCIUM GLUCONATE"
med_name_simple_mat[26,2] = " METOPROLOL TARTRATE"
med_name_simple_mat[30,2] = " METOPROLOL TARTRATE"
med_name_simple_mat[32,2] = " CALCIUM GLUCONATE"
med_name_simple_mat[36,2] = " DEXMEDETOMIDINE "
med_name_simple_mat[47,2] = " METOPROLOL TARTRATE"
med_name_simple_mat[76,2] = " EPHEDRINE"
med_name_simple_mat[83,2] = " ESMOLOL "
med_name_simple_mat[92,2] = " METOPROLOL SUCCINATE"
med_name_simple_mat[117,2] = " METOPROLOL SUCCINATE"
print(paste0("Number of mismatches: ", sum(is.na(med_name_simple_mat))))

hr_map = matrix(0, ncol = 2, nrow = length(med_name_simple))
colnames(hr_map) = c("hr", "map")
for(i in 1:nrow(med_key)) {
    simple_name_ind = med_name_simple_mat[med_name_simple_mat[,2] == med_key$id[i], 1]
    # print(med_key$id[i])
    # print(simple_name_ind)
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
print(meds_to_check); print(13)
# " ALBUMIN", " AMIODARONE ", " CALCIUM CHLORIDE", " CLEVIDIPINE ", " DILTIAZEM "
# " EPINEPHRINE ", " ESMOLOL " " KETAMINE ", " NITROGLYCERIN ", " PHENYLEPHRINE", " PROPOFOL "
# " SODIUM BICARBONATE ", " VASOPRESSIN "
# med_d = med_select_id_sub[med_select_id_sub$med_name_simple == " ESMOLOL ", ,drop = F]
# unique(med_d$Med_Name_Desc[med_d$continuous_med == 1])
# unique(med_d$Med_Name_Desc[med_d$continuous_med == 0])

med_name_simple_new = paste0(med_select_id_sub$med_name_simple, continuous_med)

med_name_simple_new[med_select_id_sub$med_name_simple == " ALBUMIN"] = "ALBUMIN_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " AMIODARONE " &
                    med_select_id_sub$continuous_med == 1] = "AMIODARONE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " AMIODARONE " &
                    med_select_id_sub$continuous_med == 0] = "AMIODARONE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " CALCIUM CHLORIDE" &
                        med_select_id_sub$continuous_med == 1] = "CALCIUM CHLORIDE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " CALCIUM CHLORIDE" &
                        med_select_id_sub$continuous_med == 0] = "CALCIUM CHLORIDE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " CLEVIDIPINE " &
                        med_select_id_sub$continuous_med == 1] = "CLEVIDIPINE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " CLEVIDIPINE " &
                        med_select_id_sub$continuous_med == 0] = "CLEVIDIPINE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " DILTIAZEM " &
                        med_select_id_sub$continuous_med == 1] = "DILTIAZEM_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " DILTIAZEM " &
                        med_select_id_sub$continuous_med == 0] = "DILTIAZEM_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " EPINEPHRINE " &
                        med_select_id_sub$continuous_med == 1] = "EPINEPHRINE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " EPINEPHRINE " &
                        med_select_id_sub$continuous_med == 0] = "EPINEPHRINE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " ESMOLOL " &
                        med_select_id_sub$continuous_med == 1] = "ESMOLOL_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " ESMOLOL " &
                        med_select_id_sub$continuous_med == 0] = "ESMOLOL_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " KETAMINE " &
                        med_select_id_sub$continuous_med == 1] = "KETAMINE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " KETAMINE " &
                        med_select_id_sub$continuous_med == 0] = "KETAMINE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " NITROGLYCERIN " &
                        med_select_id_sub$continuous_med == 1] = "NITROGLYCERIN_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " NITROGLYCERIN " &
                        med_select_id_sub$continuous_med == 0] = "NITROGLYCERIN_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " PHENYLEPHRINE" &
                        med_select_id_sub$continuous_med == 1] = "PHENYLEPHRINE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " PHENYLEPHRINE" &
                        med_select_id_sub$continuous_med == 0] = "PHENYLEPHRINE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " PROPOFOL " &
                        med_select_id_sub$continuous_med == 1] = "PROPOFOL_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " PROPOFOL " &
                        med_select_id_sub$continuous_med == 0] = "PROPOFOL_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " SODIUM BICARBONATE " &
                        med_select_id_sub$continuous_med == 1] = "SODIUM BICARBONATE_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " SODIUM BICARBONATE " &
                        med_select_id_sub$continuous_med == 0] = "SODIUM BICARBONATE_0"

med_name_simple_new[med_select_id_sub$med_name_simple == " VASOPRESSIN " &
                        med_select_id_sub$continuous_med == 1] = "VASOPRESSIN_1"
med_name_simple_new[med_select_id_sub$med_name_simple == " VASOPRESSIN " &
                        med_select_id_sub$continuous_med == 0] = "VASOPRESSIN_0"

# med_name_admin = paste0(med_select_id_sub$Med_Name_Desc, continuous_med)
med_name_admin = med_name_simple_new
unique_med_names = sort(unique(med_name_admin)); print(length(unique_med_names))
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
med_select_FINAL = med_select_FINAL[med_select_FINAL$Dose != 0, ]
med_select_FINAL = med_select_FINAL[!is.na(med_select_FINAL$Dose), ]
# *****************************************************************************

# making the strength numeric
all_med_types = unique(med_select_FINAL$med_name_admin)
strength_units = dose_units = vector(mode = 'list', length = length(all_med_types))
for(i in 1:length(strength_units)) {
    strength_units[[i]] = unique(med_select_FINAL$Strength[med_select_FINAL$med_name_admin == all_med_types[i]])
    dose_units[[i]] = unique(med_select_FINAL$Dose_Units[med_select_FINAL$med_name_admin == all_med_types[i]])
}
names(strength_units) = names(dose_units) = all_med_types

# looking at which strength_units and dose_units have more than 1 per drug
strength_check = c("ALBUMIN_0", " CALCIUM GLUCONATE0"," LABETALOL 0",
                   "AMIODARONE_0","DILTIAZEM_0"," NIFEDIPINE0",
                   " NIMODIPINE 0")
dose_check = c("NITROGLYCERIN_1", " NOREPINEPHRINE1","PROPOFOL_1", "PROPOFOL_0",
               "CALCIUM CHLORIDE_0", " CALCIUM GLUCONATE0",
               "SODIUM BICARBONATE_0", "CLEVIDIPINE_1", "PHENYLEPHRINE_0"," NIFEDIPINE0")
dose_strength_combo = unique(c(strength_check, dose_check))
print(dose_strength_combo)

med_select_FINAL$Strength[!(med_select_FINAL$med_name_admin %in% dose_strength_combo)] = 1
# **************************************************************** 
# FINISH RIGHT HERE!!!!!!!!!
# **************************************************************** 
# specific_med = med_select_FINAL[med_select_FINAL$med_name_admin =="PHENYLEPHRINE_0" , ]
# apply(specific_med[,c("Dose", "Dose_Units", "Strength")], 2, table)
# ALBUMIN ---------------------------------------------------------------------
# mL -> g (divide by 10)
alb_ind = which(med_select_FINAL$med_name_admin == "ALBUMIN_0" & 
                    med_select_FINAL$Dose_Units == "mL")
med_select_FINAL[alb_ind, "Dose"] = med_select_FINAL[alb_ind, "Dose"] / 10
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "ALBUMIN_0" & 
                    med_select_FINAL$Strength == "25 %"] = 0.25
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "ALBUMIN_0" & 
                              med_select_FINAL$Strength == "5 %"] = 0.05
# CALCIUM GLUCONATE0 ----------------------------------------------------------
cal_g_ind = which(med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0" & 
                      med_select_FINAL$Dose_Units == "mg")
med_select_FINAL[cal_g_ind, "Dose"] = med_select_FINAL[cal_g_ind, "Dose"] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0" & 
                              med_select_FINAL$Strength == "10%"] = 5
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " CALCIUM GLUCONATE0" &
                              med_select_FINAL$Strength != 5] = 1
# LABETALOL -------------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " LABETALOL 0"] = 1
# AMIODARONE -------------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "AMIODARONE_0"] = 1
# DILTIAZEM -------------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "DILTIAZEM_0"] = 1
# NIFEDIPINE ------------------------------------------------------------------
med_select_FINAL = med_select_FINAL[med_select_FINAL$med_name_admin != " NIFEDIPINE0", ]
# NIMODIPINE ------------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NIMODIPINE 0"] = 1
# NITROGLYCERIN ---------------------------------------------------------------
med_select_FINAL = med_select_FINAL[!(med_select_FINAL$med_name_admin == "NITROGLYCERIN_1" &
                                        med_select_FINAL$Dose_Units == "mcg"), ]
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "NITROGLYCERIN_1"] = 1
# NOREPINEPHRINE --------------------------------------------------------------
med_select_FINAL = med_select_FINAL[!(med_select_FINAL$med_name_admin == " NOREPINEPHRINE1" &
                                          med_select_FINAL$Dose == 4000), ]
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == " NOREPINEPHRINE1"] = 1
# PROPOFOL 1 ------------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "PROPOFOL_1"] = 1
# PROPOFOL 0 ------------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "PROPOFOL_0"] = 1
# CALCIUM CHLORIDE ------------------------------------------------------------
cc_ind = which(med_select_FINAL$med_name_admin == "CALCIUM CHLORIDE_0" &  med_select_FINAL$Dose_Units == "mg")
med_select_FINAL$Dose[cc_ind] = med_select_FINAL$Dose[cc_ind] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "CALCIUM CHLORIDE_0"] = 1
# SODIUM BICARBONATE ----------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "SODIUM BICARBONATE_0"] = 1
# CLEVIDIPINE -----------------------------------------------------------------
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "CLEVIDIPINE_1"] = 1
# PHENYLEPHRINE ---------------------------------------------------------------
p_ind = which(med_select_FINAL$med_name_admin == "PHENYLEPHRINE_0" &  med_select_FINAL$Dose_Units == "mcg")
med_select_FINAL$Dose[p_ind] = med_select_FINAL$Dose[p_ind] / 1000
med_select_FINAL$Strength[med_select_FINAL$med_name_admin == "PHENYLEPHRINE_0"] = 1

# making sure the units are the same within one medication
Strength_num = rep(0, nrow(med_select_FINAL))
for(i in names(strength_units)) {
    print(i)
    if(length(strength_units[[i]]) > 1) {
        # Currently only 1 is and that is "CARVEDILOL 3.125 MG TABLET0"
        print(i)
        Strength_num[med_select_FINAL$med_name_admin == i] = 
            as.numeric(strsplit(strength_units[[i]][1], " ")[[1]][1])
    } else {
        # First check for %
        if(length(strsplit(strength_units[[i]][1], '%')[[1]]) == 1 & str_detect(strength_units[[i]][1], '%')) {
            # This means the only unit is %
            perc = as.numeric(strsplit(strength_units[[i]][1], '%')[[1]])[1]
            perc = perc / 100
            Strength_num[med_select_FINAL$med_name_admin == i] = perc
        } else if(length(strsplit(strength_units[[i]][1], " ")[[1]]) == 0) {
            # This means no units are provided
            Strength_num[med_select_FINAL$med_name_admin == i] = 1
            
        } else if(strsplit(strength_units[[i]][1], " ")[[1]][1] == "IP") {
            # This means the unit is IP ONLY
            Strength_num[med_select_FINAL$med_name_admin == i] = 1
            
        } else {
            # Standard unit
            Strength_num[med_select_FINAL$med_name_admin == i] = 
                as.numeric(strsplit(strength_units[[i]][1], " ")[[1]][1])
        }
    }
}

med_select_FINAL = cbind(med_select_FINAL, Strength_num)
med_select_FINAL$Dose[is.na(med_select_FINAL$Dose)] = 0
med_select_FINAL$Strength_num[is.na(med_select_FINAL$Strength_num)] = 0

strength_units_update = vector(mode = 'list', length = length(all_med_types))
for(i in 1:length(strength_units_update)) {
    strength_units_update[[i]] = unique(med_select_FINAL$Strength_num[med_select_FINAL$med_name_admin == all_med_types[i]])
}
names(strength_units_update) = all_med_types


# Verifying that the dose * strength = 0 for "Stopped"
print(unique(med_select_FINAL$Dose[med_select_FINAL$Status == "Stopped"]))
print(table(med_select_FINAL$status_med[med_select_FINAL$Dose == 0]))
dose_0_instance = med_select_FINAL[med_select_FINAL$Dose == 0, ]

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

save(Dn_omega, file = 'Data/Dn_omega.rda')