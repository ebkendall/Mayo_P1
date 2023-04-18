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
                         "last_time" = rep(0, length(med_select_id$Med_Name_Desc)))
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

clean_meds$last_time = clean_meds$time + clean_meds$onset + clean_meds$offset

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
    clean_meds_sub = clean_meds[clean_meds[,'id'] == select_id[i], ]
    
    # subsetting for the particular type of medication
    hr_uppers = clean_meds_sub[clean_meds_sub$hr == 1, , drop = F]; print(nrow(hr_uppers))
    hr_downer = clean_meds_sub[clean_meds_sub$hr == -1, , drop = F]; print(nrow(hr_downer))
    map_uppers = clean_meds_sub[clean_meds_sub$map == 1, , drop = F]; print(nrow(map_uppers))
    map_downer = clean_meds_sub[clean_meds_sub$map == -1, , drop = F]; print(nrow(map_downer))
    for(j in 1:nrow(sub_data)) {
        
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
                                     (hr_uppers$time + hr_uppers$onset + hr_uppers$offset) > sub_data[j,"time"])
        hr_downer_sub_off = which((hr_downer$time + hr_downer$onset) < sub_data[j,"time"] & 
                                      (hr_downer$time + hr_downer$onset + hr_downer$offset) > sub_data[j,"time"])
        map_uppers_sub_off = which((map_uppers$time + map_uppers$onset) < sub_data[j,"time"] & 
                                       (map_uppers$time + map_uppers$onset + map_uppers$offset) > sub_data[j,"time"])
        map_downer_sub_off = which((map_downer$time + map_downer$onset) < sub_data[j,"time"] & 
                                       (map_downer$time + map_downer$onset + map_downer$offset) > sub_data[j,"time"])
        
        
        if(length(hr_uppers_sub_on) > 0) med_format[[i]][j, 1] = sum(sub_data[j,"time"] - hr_uppers$time[hr_uppers_sub_on])
        if(length(hr_downer_sub_on) > 0) med_format[[i]][j, 3] = sum(sub_data[j,"time"] - hr_downer$time[hr_downer_sub_on])
        if(length(map_uppers_sub_on) > 0) med_format[[i]][j, 5] = sum(sub_data[j,"time"] - map_uppers$time[map_uppers_sub_on])
        if(length(map_downer_sub_on) > 0) med_format[[i]][j, 7] = sum(sub_data[j,"time"] - map_downer$time[map_downer_sub_on])
        
        if(length(hr_uppers_sub_off) > 0) {
            med_format[[i]][j, 1] = med_format[[i]][j, 1] + sum(hr_uppers$onset[hr_uppers_sub_off])
            med_format[[i]][j, 2] = sum(sub_data[j,"time"] - hr_uppers$time[hr_uppers_sub_off] - hr_uppers$onset[hr_uppers_sub_off])
        }
        if(length(hr_downer_sub_off) > 0) {
            med_format[[i]][j, 3] = med_format[[i]][j, 3] + sum(hr_downer$onset[hr_downer_sub_off])
            med_format[[i]][j, 4] = sum(sub_data[j,"time"] - hr_downer$time[hr_downer_sub_off] - hr_downer$onset[hr_downer_sub_off])
        }
        if(length(map_uppers_sub_off) > 0) {
            med_format[[i]][j, 5] = med_format[[i]][j, 5] + sum(map_uppers$onset[map_uppers_sub_off])
            med_format[[i]][j, 6] = sum(sub_data[j,"time"] - map_uppers$time[map_uppers_sub_off] - map_uppers$onset[map_uppers_sub_off])
        }
        if(length(map_downer_sub_off) > 0) {
            med_format[[i]][j, 7] = med_format[[i]][j, 7] + sum(map_downer$onset[map_downer_sub_off])
            med_format[[i]][j, 8] = sum(sub_data[j,"time"] - map_downer$time[map_downer_sub_off] - map_downer$onset[map_downer_sub_off])
        }
    }
}

save(med_format, file = 'Data/med_format.rda')
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

