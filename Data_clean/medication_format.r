jw15 = read.csv("Data/_raw_data_new/jw15b.csv")
jw16 = read.csv("Data/_raw_data_new/jw16b.csv")
jw17 = read.csv("Data/_raw_data_new/jw17b.csv")
jw18 = read.csv("Data/_raw_data_new/jw18b.csv")
jw19 = read.csv("Data/_raw_data_new/jw19b.csv")

med_key = read.csv('Med_chart.csv', na.strings = "")
colnames(med_key) = c('id', 'hr', 'map', 'onset', 'offset', 'time_check', 'X1')

load('Data/data_format_FULL_48hr_update_RBC_sub.rda')
pace_id = c(18075, 108825, 110750, 125025, 173750, 260100, 304700, 307225, 310100,
            382450, 429375, 516150, 533075, 666750, 677225, 732525, 763050, 767500, 
            769025, 777175, 794900, 799125, 819225)
data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
select_id = unique(data_format[,"EID"])
med_select_id = jw15[jw15$key %in% select_id, ]
med_select_id = rbind(med_select_id, jw16[jw16$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw17[jw17$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw18[jw18$key %in% select_id, ])
med_select_id = rbind(med_select_id, jw19[jw19$key %in% select_id, ])
med_select_id$administered_dtm = med_select_id$administered_dtm / 60
med_select_id = med_select_id[!is.na(med_select_id[,1]),]


# cleaning this information
library(stringr)
clean_meds = data.frame( "id" = med_select_id$key,
                         "med" = med_select_id$Med_Name_Desc,
                         "time" = med_select_id$administered_dtm,
                         "hr" = rep(0, length(med_select_id$Med_Name_Desc)),
                         "map" = rep(0, length(med_select_id$Med_Name_Desc)),
                         "onset" = rep(0, length(med_select_id$Med_Name_Desc)),
                         "offset" = rep(0, length(med_select_id$Med_Name_Desc)),
                         "row_med_key" = rep(0, length(med_select_id$Med_Name_Desc)),
                         "T_1_max" = rep(0, length(med_select_id$Med_Name_Desc)),
                         "T_2_max" = rep(0, length(med_select_id$Med_Name_Desc)))
for(i in 1:nrow(med_key)) {
    ind = grep(med_key$id[i], clean_meds$med)
    if(length(ind) > 0) {
        clean_meds[ind, c('hr', 'map')] = med_key[i, c('hr', 'map')]
        clean_meds[ind, c('onset', 'offset')] = med_key[i, c("onset", "offset")]
        clean_meds[ind, 'row_med_key'] = med_key[i, 'id']
    }
}

# Double checking no missingness
print(sum(is.na(clean_meds$onset)))

clean_meds$T_1_max = clean_meds$time + clean_meds$onset
clean_meds$T_2_max = clean_meds$time + clean_meds$onset + clean_meds$offset

# Reordering to be in numerical order
for(i in unique(clean_meds$id)) {
    sub_ind = clean_meds[clean_meds$id == i, , drop =F]
    re_ordered_sub = sub_ind[order(sub_ind$time), ]
    clean_meds[clean_meds$id == i, ] = re_ordered_sub
}

# Removing the words for HR and MAP
clean_meds$hr[is.na(clean_meds$hr)] = 0
clean_meds$hr[clean_meds$hr == "Up"] = 1
clean_meds$hr[clean_meds$hr == "Down"] = -1
clean_meds$map[is.na(clean_meds$map)] = 0
clean_meds$map[clean_meds$map == "Up"] = 1
clean_meds$map[clean_meds$map == "Down"] = -1

save(clean_meds, file = 'Data/clean_meds.rda')

# Medication formatting for the final dataset
med_format = vector(mode = "list", length = length(select_id))
for(i in 1:length(med_format)) {
    print(paste0("ID: ", i))
    med_format[[i]] = matrix(0, nrow = sum(data_format[,'EID'] == select_id[i]), ncol = 8)
    sub_data = data_format[data_format[,'EID'] == select_id[i], ]
    
    clean_meds_sub_pre = clean_meds[clean_meds[,'id'] == select_id[i], ]
    if(nrow(clean_meds_sub_pre) >= 1 & clean_meds_sub_pre$med[1] != "") {
        # Remove any medication who's total effect time is done by t_1
        clean_meds_sub_pre = clean_meds_sub_pre[clean_meds_sub_pre[,"T_2_max"] > sub_data[1, "time"], ]
        
        for(j in 1:nrow(sub_data)) {
            
            clean_meds_sub = clean_meds_sub_pre[clean_meds_sub_pre[,'time'] <= sub_data[j, "time"], ]
            
            # subsetting for the particular type of medication
            hr_uppers = clean_meds_sub[clean_meds_sub$hr == 1, , drop = F]
            hr_downer = clean_meds_sub[clean_meds_sub$hr == -1, , drop = F]
            map_uppers = clean_meds_sub[clean_meds_sub$map == 1, , drop = F]
            map_downer = clean_meds_sub[clean_meds_sub$map == -1, , drop = F]
            
            if(nrow(hr_uppers) > 0) {
                diff1 = hr_uppers$T_1_max - hr_uppers$time
                diff2 = sub_data[j,"time"] - hr_uppers$time
                med_format[[i]][j, 1] = sum(apply(cbind(diff1, diff2), 1, min))
                
                check_ind = which(hr_uppers$T_1_max < sub_data[j,"time"])
                if(sum(hr_uppers$T_1_max < sub_data[j,"time"]) > 0) {
                    diff1_off = hr_uppers$T_2_max[check_ind] - hr_uppers$T_1_max[check_ind]
                    diff2_off = sub_data[j,"time"] - hr_uppers$T_1_max[check_ind]
                    med_format[[i]][j, 2] = sum(apply(cbind(diff1_off, diff2_off), 1, min))
                }
            }
            
            if(nrow(hr_downer) > 0) {
                diff1 = hr_downer$T_1_max - hr_downer$time
                diff2 = sub_data[j,"time"] - hr_downer$time
                med_format[[i]][j, 3] = sum(apply(cbind(diff1, diff2), 1, min))
                
                check_ind = which(hr_downer$T_1_max < sub_data[j,"time"])
                if(length(check_ind) > 0) {
                    diff1_off = hr_downer$T_2_max[check_ind] - hr_downer$T_1_max[check_ind]
                    diff2_off = sub_data[j,"time"] - hr_downer$T_1_max[check_ind]
                    med_format[[i]][j, 4] = sum(apply(cbind(diff1_off, diff2_off), 1, min))
                }
            }
            
            if(nrow(map_uppers) > 0) {
                diff1 = map_uppers$T_1_max - map_uppers$time
                diff2 = sub_data[j,"time"] - map_uppers$time
                med_format[[i]][j, 5] = sum(apply(cbind(diff1, diff2), 1, min))
                
                check_ind = which(map_uppers$T_1_max < sub_data[j,"time"])
                if(length(check_ind) > 0) {
                    diff1_off = map_uppers$T_2_max[check_ind] - map_uppers$T_1_max[check_ind]
                    diff2_off = sub_data[j,"time"] - map_uppers$T_1_max[check_ind]
                    med_format[[i]][j, 6] = sum(apply(cbind(diff1_off, diff2_off), 1, min))
                }
            }
            
            if(nrow(map_downer) > 0) {
                diff1 = map_downer$T_1_max - map_downer$time
                diff2 = sub_data[j,"time"] - map_downer$time
                med_format[[i]][j, 7] = sum(apply(cbind(diff1, diff2), 1, min))
                
                check_ind = which(map_downer$T_1_max < sub_data[j,"time"])
                if(length(check_ind) > 0) {
                    diff1_off = map_downer$T_2_max[check_ind] - map_downer$T_1_max[check_ind]
                    diff2_off = sub_data[j,"time"] - map_downer$T_1_max[check_ind]
                    med_format[[i]][j, 8] = sum(apply(cbind(diff1_off, diff2_off), 1, min))
                }
            }
        }
    } else print(paste0("ID: ", i, " no meds"))
}

save(med_format, file = 'Data/med_format.rda')

# Structuring this for the runfile
Dn_omega = vector(mode = 'list', length = length(select_id))
for (i in 1:length(Dn_omega)) {
    print(paste0("ID: ", i))
    n_i = sum(data_format[,'EID'] == select_id[i])
    Dn_omega[[i]] = matrix(0, nrow = 4*n_i, ncol = 8)
    Dn_omega[[i]][(n_i+1):(2*n_i), 1:4] = med_format[[i]][,1:4]
    Dn_omega[[i]][(2*n_i+1):(3*n_i), 1:4] = med_format[[i]][,5:8]
}

save(Dn_omega, file = 'Data/Dn_omega.rda')
