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

# Medication formatting for the final dataset
med_format = vector(mode = "list", length = length(select_id))
for(i in 1:length(med_format)) {
    print(paste0("ID: ", i))
    med_format[[i]] = matrix(0, nrow = sum(data_format[,'EID'] == select_id[i]), ncol = 8)
    sub_data = data_format[data_format[,'EID'] == select_id[i], ]
    
    clean_meds_sub_pre = clean_meds[clean_meds[,'id'] == select_id[i], ]
    # Remove any medication who's total effect time is done by t_1
    clean_meds_sub_pre = clean_meds_sub_pre[clean_meds_sub_pre[,"T_2_max"] > sub_data[1, "time"], ]
    
    for(j in 1:nrow(sub_data)) {
        
        clean_meds_sub = clean_meds_sub_pre[clean_meds_sub_pre[,'time'] <= sub_data[j, "time"], ]
        
        # subsetting for the particular type of medication
        hr_uppers = clean_meds_sub[clean_meds_sub$hr == 1, , drop = F]; print(nrow(hr_uppers))
        hr_downer = clean_meds_sub[clean_meds_sub$hr == -1, , drop = F]; print(nrow(hr_downer))
        map_uppers = clean_meds_sub[clean_meds_sub$map == 1, , drop = F]; print(nrow(map_uppers))
        map_downer = clean_meds_sub[clean_meds_sub$map == -1, , drop = F]; print(nrow(map_downer))
        
        # First component ------------------------------------------------------
        # Onset
        hr_uppers_sub_on = which(hr_uppers$time < sub_data[j,"time"] & 
                                      (hr_uppers$time + hr_uppers$onset) >= sub_data[j,"time"])
        hr_downer_sub_on = which(hr_downer$time < sub_data[j,"time"] & 
                                  (hr_downer$time + hr_downer$onset) >= sub_data[j,"time"])
        map_uppers_sub_on = which(map_uppers$time < sub_data[j,"time"] & 
                                  (map_uppers$time + map_uppers$onset) >= sub_data[j,"time"])
        map_downer_sub_on = which(map_downer$time < sub_data[j,"time"] & 
                                  (map_downer$time + map_downer$onset) >= sub_data[j,"time"])
        
        # Offset
        hr_uppers_sub_off = which((hr_uppers$time + hr_uppers$onset) < sub_data[j,"time"] & 
                                     (hr_uppers$time + hr_uppers$onset + hr_uppers$offset) >= sub_data[j,"time"])
        hr_downer_sub_off = which((hr_downer$time + hr_downer$onset) < sub_data[j,"time"] & 
                                      (hr_downer$time + hr_downer$onset + hr_downer$offset) >= sub_data[j,"time"])
        map_uppers_sub_off = which((map_uppers$time + map_uppers$onset) < sub_data[j,"time"] & 
                                       (map_uppers$time + map_uppers$onset + map_uppers$offset) >= sub_data[j,"time"])
        map_downer_sub_off = which((map_downer$time + map_downer$onset) < sub_data[j,"time"] & 
                                       (map_downer$time + map_downer$onset + map_downer$offset) >= sub_data[j,"time"])
        
        # Onset ---------------------------------------------------------------
        if(length(hr_uppers_sub_on) > 0) {
            if(j == 1) {
                med_format[[i]][j, 1] = sum(sub_data[j,"time"] - hr_uppers$time[hr_uppers_sub_on])
            } else {
                admin_diff = sub_data[j,"time"] - hr_uppers$time[hr_uppers_sub_on]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 1] = med_format[[i]][j-1, 1] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 1] = med_format[[i]][j-1, 1]}
        
        if(length(hr_downer_sub_on) > 0) {
            if(j == 1) {
                med_format[[i]][j, 3] = sum(sub_data[j,"time"] - hr_downer$time[hr_downer_sub_on])
            } else {
                admin_diff = sub_data[j,"time"] - hr_downer$time[hr_downer_sub_on]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 3] = med_format[[i]][j-1, 3] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 3] = med_format[[i]][j-1, 3]}
            
        if(length(map_uppers_sub_on) > 0) {
            if(j == 1) {
                med_format[[i]][j, 5] = sum(sub_data[j,"time"] - map_uppers$time[map_uppers_sub_on])
            } else {
                admin_diff = sub_data[j,"time"] - map_uppers$time[map_uppers_sub_on]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 5] = med_format[[i]][j-1, 5] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 5] = med_format[[i]][j-1, 5]}
            
        if(length(map_downer_sub_on) > 0) {
            if(j == 1) {
                med_format[[i]][j, 7] = sum(sub_data[j,"time"] - map_downer$time[map_downer_sub_on])
            } else {
                admin_diff = sub_data[j,"time"] - map_downer$time[map_downer_sub_on]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 7] = med_format[[i]][j-1, 7] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 7] = med_format[[i]][j-1, 7]}
        
        # Offset ---------------------------------------------------------------
        if(length(hr_uppers_sub_off) > 0) {
            if(j == 1) {
                med_format[[i]][j, 2] = sum(sub_data[j,"time"] - hr_uppers$T_1_max[hr_uppers_sub_off])
            } else {
                admin_diff = sub_data[j,"time"] - hr_uppers$T_1_max[hr_uppers_sub_off]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 2] = med_format[[i]][j-1, 2] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 2] = med_format[[i]][j-1, 2]}
        
        if(length(hr_downer_sub_off) > 0) {
            if(j == 1) {
                med_format[[i]][j, 4] = sum(sub_data[j,"time"] - hr_downer$T_1_max[hr_downer_sub_off])
            } else {
                admin_diff = sub_data[j,"time"] - hr_downer$T_1_max[hr_downer_sub_off]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 4] = med_format[[i]][j-1, 4] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 4] = med_format[[i]][j-1, 4]}
        
        if(length(map_uppers_sub_off) > 0) {
            if(j == 1) {
                med_format[[i]][j, 6] = sum(sub_data[j,"time"] - map_uppers$T_1_max[map_uppers_sub_off])
            } else {
                admin_diff = sub_data[j,"time"] - map_uppers$T_1_max[map_uppers_sub_off]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 6] = med_format[[i]][j-1, 6] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 6] = med_format[[i]][j-1, 6]}
        
        if(length(map_downer_sub_off) > 0) {
            if(j == 1) {
                med_format[[i]][j, 8] = sum(sub_data[j,"time"] - map_downer$T_1_max[map_downer_sub_off])
            } else {
                admin_diff = sub_data[j,"time"] - map_downer$T_1_max[map_downer_sub_off]
                temp_diff  = rep(sub_data[j,"time"] - sub_data[j-1,"time"], length(admin_diff))
                med_format[[i]][j, 8] = med_format[[i]][j-1, 8] + sum(apply(cbind(admin_diff, temp_diff), 1, min))
            }
        } else {if(j > 1) med_format[[i]][j, 8] = med_format[[i]][j-1, 8]}
        
        
        # Second component -----------------------------------------------------
        if(j == 1) {
            hr_up_max_on    = which(hr_uppers$T_1_max < sub_data[j,"time"] & hr_uppers$T_1_max >= 0)
            hr_down_max_on  = which(hr_downer$T_1_max < sub_data[j,"time"] & hr_downer$T_1_max >= 0)
            map_up_max_on   = which(map_uppers$T_1_max < sub_data[j,"time"] & map_uppers$T_1_max >= 0)
            map_down_max_on = which(map_downer$T_1_max < sub_data[j,"time"] & map_downer$T_1_max >= 0)
            # Don't need to worry about the case of offset times finishing before time 1
            hr_up_max_off = hr_down_max_off = map_up_max_off = map_down_max_off = integer(0)
        } else {
            hr_up_max_on    = which(hr_uppers$T_1_max < sub_data[j,"time"] & hr_uppers$T_1_max > sub_data[j-1,"time"])
            hr_down_max_on  = which(hr_downer$T_1_max < sub_data[j,"time"] & hr_downer$T_1_max > sub_data[j-1,"time"])
            map_up_max_on   = which(map_uppers$T_1_max < sub_data[j,"time"] & map_uppers$T_1_max > sub_data[j-1,"time"])
            map_down_max_on = which(map_downer$T_1_max < sub_data[j,"time"] & map_downer$T_1_max > sub_data[j-1,"time"])
            
            hr_up_max_off    = which(hr_uppers$T_2_max < sub_data[j,"time"] & hr_uppers$T_2_max > sub_data[j-1,"time"])
            hr_down_max_off  = which(hr_downer$T_2_max < sub_data[j,"time"] & hr_downer$T_2_max > sub_data[j-1,"time"])
            map_up_max_off   = which(map_uppers$T_2_max < sub_data[j,"time"] & map_uppers$T_2_max > sub_data[j-1,"time"])
            map_down_max_off = which(map_downer$T_2_max < sub_data[j,"time"] & map_downer$T_2_max > sub_data[j-1,"time"])
        }
        
        # Onset ---------------------------------------------------------------
        if(length(hr_up_max_on) > 0) {
            if(j == 1) med_format[[i]][j, 1] = med_format[[i]][j, 1] + sum(hr_uppers$T_1_max[hr_up_max_on])
            else {
                t_1_diff = hr_uppers$T_1_max[hr_up_max_on] - hr_uppers$time[hr_up_max_on]
                temp_diff  = hr_uppers$T_1_max[hr_up_max_on] - sub_data[j-1,"time"]
                med_format[[i]][j, 1] = med_format[[i]][j, 1] + sum(apply(cbind(t_1_diff, temp_diff), 1, min))
            }
        }
        if(length(hr_down_max_on) > 0) {
            if(j == 1) med_format[[i]][j, 3] = med_format[[i]][j, 3] + sum(hr_downer$T_1_max[hr_down_max_on])
            else {
                t_1_diff = hr_downer$T_1_max[hr_down_max_on] - hr_downer$time[hr_down_max_on]
                temp_diff  = hr_downer$T_1_max[hr_down_max_on] - sub_data[j-1,"time"]
                med_format[[i]][j, 3] = med_format[[i]][j, 3] + sum(apply(cbind(t_1_diff, temp_diff), 1, min))
            }
        }
        if(length(map_up_max_on) > 0) {
            if(j == 1) med_format[[i]][j, 5] = med_format[[i]][j, 5] + sum(map_uppers$T_1_max[map_up_max_on])
            else {
                t_1_diff = map_uppers$T_1_max[map_up_max_on] - map_uppers$time[map_up_max_on]
                temp_diff  = map_uppers$T_1_max[map_up_max_on] - sub_data[j-1,"time"]
                med_format[[i]][j, 5] = med_format[[i]][j, 5] + sum(apply(cbind(t_1_diff, temp_diff), 1, min))
            }
        }
        if(length(map_down_max_on) > 0) {
            if(j == 1) med_format[[i]][j, 7] = med_format[[i]][j, 7] + sum(map_downer$T_1_max[map_down_max_on])
            else {
                t_1_diff = map_downer$T_1_max[map_down_max_on] - map_downer$time[map_down_max_on]
                temp_diff  = map_downer$T_1_max[map_down_max_on] - sub_data[j-1,"time"]
                med_format[[i]][j, 7] = med_format[[i]][j, 7] + sum(apply(cbind(t_1_diff, temp_diff), 1, min))
            }
        }
        
        # Offset ---------------------------------------------------------------
        if(length(hr_up_max_off) > 0) {
            if(j == 1) med_format[[i]][j, 2] = med_format[[i]][j, 2] + sum(hr_uppers$T_2_max[hr_up_max_off])
            else {
                t_2_diff = hr_uppers$T_2_max[hr_up_max_off] - hr_uppers$time[hr_up_max_off]
                temp_diff  = hr_uppers$T_2_max[hr_up_max_off] - sub_data[j-1,"time"]
                med_format[[i]][j, 2] = med_format[[i]][j, 2] + sum(apply(cbind(t_2_diff, temp_diff), 1, min))
            }
        }
        if(length(hr_down_max_off) > 0) {
            if(j == 1) med_format[[i]][j, 4] = med_format[[i]][j, 4] + sum(hr_downer$T_2_max[hr_down_max_off])
            else {
                t_2_diff = hr_downer$T_2_max[hr_down_max_off] - hr_downer$time[hr_down_max_off]
                temp_diff  = hr_downer$T_2_max[hr_down_max_off] - sub_data[j-1,"time"]
                med_format[[i]][j, 4] = med_format[[i]][j, 4] + sum(apply(cbind(t_2_diff, temp_diff), 1, min))
            }
        }
        if(length(map_up_max_off) > 0) {
            if(j == 1) med_format[[i]][j, 6] = med_format[[i]][j, 6] + sum(map_uppers$T_2_max[map_up_max_off])
            else {
                t_2_diff = map_uppers$T_2_max[map_up_max_off] - map_uppers$time[map_up_max_off]
                temp_diff  = map_uppers$T_2_max[map_up_max_off] - sub_data[j-1,"time"]
                med_format[[i]][j, 6] = med_format[[i]][j, 6] + sum(apply(cbind(t_2_diff, temp_diff), 1, min))
            }
        }
        if(length(map_down_max_off) > 0) {
            if(j == 1) med_format[[i]][j, 8] = med_format[[i]][j, 8] + sum(map_downer$T_2_max[map_down_max_off])
            else {
                t_2_diff = map_downer$T_2_max[map_down_max_off] - map_downer$time[map_down_max_off]
                temp_diff  = map_downer$T_2_max[map_down_max_off] - sub_data[j-1,"time"]
                med_format[[i]][j, 8] = med_format[[i]][j, 8] + sum(apply(cbind(t_2_diff, temp_diff), 1, min))
            }
        }
    }
}

save(med_format, file = 'Data/med_format.rda')

med_format = vector(mode = "list", length = length(select_id))
for(i in 1:length(med_format)) {
    print(paste0("ID: ", i))
    med_format[[i]] = matrix(0, nrow = sum(data_format[,'EID'] == select_id[i]), ncol = 8)
    sub_data = data_format[data_format[,'EID'] == select_id[i], ]
    
    clean_meds_sub_pre = clean_meds[clean_meds[,'id'] == select_id[i], ]
    # Remove any medication who's total effect time is done by t_1
    clean_meds_sub_pre = clean_meds_sub_pre[clean_meds_sub_pre[,"T_2_max"] > sub_data[1, "time"], ]
    
    for(j in 1:nrow(sub_data)) {
        
        clean_meds_sub = clean_meds_sub_pre[clean_meds_sub_pre[,'time'] <= sub_data[j, "time"], ]
        
        # subsetting for the particular type of medication
        hr_uppers = clean_meds_sub[clean_meds_sub$hr == 1, , drop = F]; print(nrow(hr_uppers))
        hr_downer = clean_meds_sub[clean_meds_sub$hr == -1, , drop = F]; print(nrow(hr_downer))
        map_uppers = clean_meds_sub[clean_meds_sub$map == 1, , drop = F]; print(nrow(map_uppers))
        map_downer = clean_meds_sub[clean_meds_sub$map == -1, , drop = F]; print(nrow(map_downer))
        
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
}
# save(med_key, file = 'Data/med_key.rda')
# if(j==1) {
#     med_before_j = which(clean_meds_sub$time <= sub_data[j,"time"] & clean_meds_sub$time > 0)
# } else {
#     med_before_j = which(clean_meds_sub$time <= sub_data[j,"time"] & clean_meds_sub$time > sub_data[j-1,"time"])
# }
# if(length(med_before_j) > 0) {
#     med_format[[i]][[j]] = matrix(c(clean_meds_sub$row_med_key[med_before_j], clean_meds_sub$time[med_before_j]), ncol = 2)
#     colnames(med_format[[i]][[j]]) = c('id', 't_star')
# } else {
#     med_format[[i]][[j]] = NA
# }
# ITEMS NEEDED for routine: med_format, data_format, and med_key


# difference in order time and admin time
# what is meant by onset and offset time

