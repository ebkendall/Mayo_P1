# This is the Rscript that will format the data from all of the csv files into
# two lists. One list contains the covariate information for each patient, and
# the other list contains the longitudinal information about each patient

# (1) Load all of the data -----------------------------------------------------
# This data is all covariates and vitals, NO MEDICATIONS
print('loading data')

jw1 = read.csv("Data/_raw_data_new/jw1b.csv")
jw2 = read.csv("Data/_raw_data_new/jw2b.csv")
jw3 = read.csv("Data/_raw_data_new/jw3b.csv")
jw4 = read.csv("Data/_raw_data_new/jw4b.csv")
jw5 = read.csv("Data/_raw_data_new/jw5b.csv")
jw6 = read.csv("Data/_raw_data_new/jw6b.csv")
jw7 = read.csv("Data/_raw_data_new/jw7b.csv")
jw8 = read.csv("Data/_raw_data_new/jw8b.csv")
jw9 = read.csv("Data/_raw_data_new/jw9b.csv")
jw10 = read.csv("Data/_raw_data_new/jw10b.csv")
jw11 = read.csv("Data/_raw_data_new/jw11b.csv")
jw12 = read.csv("Data/_raw_data_new/jw12b.csv")
jw13 = read.csv("Data/_raw_data_new/jw13b.csv")
jw14 = read.csv("Data/_raw_data_new/jw14b.csv")
jw_t = read.csv("Data/_raw_data_new/jw_transfusions.csv")

# (2) Define the training and testing data -------------------------------------
all_keys = jw1$key
load('Data/test_keys.rda')
all_keys = all_keys[-which(all_keys %in% test_keys)]
save(all_keys, file = 'Data/all_keys.rda')

# (3) Get the baseline info ----------------------------------------------------
print('getting baseline covariates')

rowInd = which(jw1$key %in% all_keys)

cov_info = jw1[rowInd, c('key', 'age', 'cardio_shock')]
cov_info = cbind(cov_info, jw7[rowInd, c('Gender', 'hemorrhage')])
cov_info = cbind(cov_info, jw8[rowInd, c('hypovo_shock', 'icu_death')])
cov_info = cbind(cov_info, jw11$los_icu[rowInd])
cov_info = cbind(cov_info, jw2[rowInd, c('other_shock', 'Race_Name')])
cov_info = cbind(cov_info, jw6$septic_shock[rowInd])
cov_info = cbind(cov_info, jw14[rowInd, c("admit_service", "icu_type")])

colnames(cov_info)[c(8, 11)] = c('los_icu', 'septic_shock') # LOS = "length of stay"

save(cov_info, file = "Data/cov_info.rda")

# (4) Get the longitudinal info ------------------------------------------------
set.seed(2022)

long_data_rd = vector(mode = "list", length = length(all_keys))
list_ind = 1

timing_issues = NULL

for (key_num in all_keys) {
    
    print(paste0(key_num, ", ind ", list_ind))
    
    sub_long_list = vector(mode = 'list', length = 24)
    main_ind = 1
    
    # jw1 & jw13 -----------------------------------------------------------------
    time_names = c('albumin_datetime1', 'albumin_datetime57', 
                   'alp_datetime1', 'alp_datetime56',
                   'alt_datetime1', 'alt_datetime57',
                   'ast_datetime1', 'ast_datetime57',
                   'bili_datetime1', 'bili_datetime54')
    result_names = c('resultn_albumin1', 'resultn_albumin57',
                     'resultn_alp1', 'resultn_alp56',
                     'resultn_alt1', 'resultn_alt57', 
                     'resultn_ast1', 'resultn_ast57',
                     'resultn_bili1', 'resultn_bili54')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw1) == time_names[i])
        e_t_ind = which(colnames(jw1) == time_names[i+1])
        s_ind = which(colnames(jw13) == result_names[i])
        e_ind = which(colnames(jw13) == result_names[i+1])
        
        info_result = cbind(t(jw1[jw1$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw13[jw13$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw7 & jw13 -----------------------------------------------------------------
    time_names = c('creat_datetime1', 'creat_datetime140',
                   'etco_datetime1', 'etco_datetime1716')
    result_names = c('resultn_creat1', 'resultn_creat140',
                     'resultn_etco1', 'resultn_etco1716')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw7) == time_names[i])
        e_t_ind = which(colnames(jw7) == time_names[i+1])
        s_ind = which(colnames(jw13) == result_names[i])
        e_ind = which(colnames(jw13) == result_names[i+1])
        
        info_result = cbind(t(jw7[jw7$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw13[jw13$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw7 & jw3 ------------------------------------------------------------------
    time_names = c('hemo_datetime1', 'hemo_datetime94')
    result_names = c('resultn_hemo1', 'resultn_hemo94')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw7) == time_names[i])
        e_t_ind = which(colnames(jw7) == time_names[i+1])
        s_ind = which(colnames(jw3) == result_names[i])
        e_ind = which(colnames(jw3) == result_names[i+1])
        
        info_result = cbind(t(jw7[jw7$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw3[jw3$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw7 & jw2 ------------------------------------------------------------------
    time_names = c('fio2_datetime1', 'fio2_datetime4278')
    result_names = c('resultn_fio1', 'resultn_fio4278')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw7) == time_names[i])
        e_t_ind = which(colnames(jw7) == time_names[i+1])
        s_ind = which(colnames(jw2) == result_names[i])
        e_ind = which(colnames(jw2) == result_names[i+1])
        
        info_result = cbind(t(jw7[jw7$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw2[jw2$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw8 & jw3 ------------------------------------------------------------------
    time_names = c('inr_datetime1', 'inr_datetime45')
    result_names = c('resultn_inr1', 'resultn_inr45')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw8) == time_names[i])
        e_t_ind = which(colnames(jw8) == time_names[i+1])
        s_ind = which(colnames(jw3) == result_names[i])
        e_ind = which(colnames(jw3) == result_names[i+1])
        
        info_result = cbind(t(jw8[jw8$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw3[jw3$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw8 & jw2 ------------------------------------------------------------------
    time_names = c('hr_datetime1', 'hr_datetime4612')
    result_names = c('resultn_hr1', 'resultn_hr4612')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw8) == time_names[i])
        e_t_ind = which(colnames(jw8) == time_names[i+1])
        s_ind = which(colnames(jw2) == result_names[i])
        e_ind = which(colnames(jw2) == result_names[i+1])
        
        info_result = cbind(t(jw8[jw8$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw2[jw2$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw11 & jw3 -----------------------------------------------------------------
    time_names = c('map_datetime1', 'map_datetime4575')
    result_names = c('resultn_map1', 'resultn_map4575')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw11) == time_names[i])
        e_t_ind = which(colnames(jw11) == time_names[i+1])
        s_ind = which(colnames(jw3) == result_names[i])
        e_ind = which(colnames(jw3) == result_names[i+1])
        
        info_result = cbind(t(jw11[jw11$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw3[jw3$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw2 & jw4 ------------------------------------------------------------------
    time_names = c('oxyamount_datetime1','oxyamount_datetime1816',
                   'pao2_datetime1', 'pao2_datetime141',
                   'plate_datetime1', 'plate_datetime63')
    result_names = c('resultn_oxyamount1', 'resultn_oxyamount1816',
                     'resultn_pao1', 'resultn_pao141',
                     'resultn_plate1', 'resultn_plate63')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw2) == time_names[i])
        e_t_ind = which(colnames(jw2) == time_names[i+1])
        s_ind = which(colnames(jw4) == result_names[i])
        e_ind = which(colnames(jw4) == result_names[i+1])
        
        info_result = cbind(t(jw2[jw2$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw4[jw4$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw9 & jw4 ------------------------------------------------------------------
    time_names = c('resprate_datetime1', 'resprate_datetime4682')
    result_names = c('resultn_resprate1', 'resultn_resprate4682')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw9) == time_names[i])
        e_t_ind = which(colnames(jw9) == time_names[i+1])
        s_ind = which(colnames(jw4) == result_names[i])
        e_ind = which(colnames(jw4) == result_names[i+1])
        
        info_result = cbind(t(jw9[jw9$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw4[jw4$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw6 & jw5 ------------------------------------------------------------------
    time_names = c('spo2_datetime1', 'spo2_datetime4588')
    result_names = c('resultn_spo1', 'resultn_spo4588')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw6) == time_names[i])
        e_t_ind = which(colnames(jw6) == time_names[i+1])
        s_ind = which(colnames(jw5) == result_names[i])
        e_ind = which(colnames(jw5) == result_names[i+1])
        
        info_result = cbind(t(jw6[jw6$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw5[jw5$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw12 & jw5 -----------------------------------------------------------------
    time_names = c('temp_datetime1', 'temp_datetime1181',
                   'tidalex_datetime1', 'tidalex_datetime1425')
    result_names = c('resultn_temp1', 'resultn_temp1181',
                     'resultn_tidalex1', 'resultn_tidalex1425')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw12) == time_names[i])
        e_t_ind = which(colnames(jw12) == time_names[i+1])
        s_ind = which(colnames(jw5) == result_names[i])
        e_ind = which(colnames(jw5) == result_names[i+1])
        
        info_result = cbind(t(jw12[jw12$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw5[jw5$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw10 & jw5 -----------------------------------------------------------------
    time_names = c('tidalin_datetime1', 'tidalin_datetime1425',
                   'troponin_datetime1', 'troponin_datetime11',
                   'vent_datetime1', 'vent_datetime877')
    result_names = c('resultn_tidalin1', 'resultn_tidalin1425',
                     'resultn_troponin1', 'resultn_troponin11',
                     'resultn_vent1', 'resultn_vent877')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw10) == time_names[i])
        e_t_ind = which(colnames(jw10) == time_names[i+1])
        s_ind = which(colnames(jw5) == result_names[i])
        e_ind = which(colnames(jw5) == result_names[i+1])
        
        info_result = cbind(t(jw10[jw10$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw5[jw5$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw10 & jw10 ----------------------------------------------------------------
    time_names = c('tranfus_rbc_datetime1', 'tranfus_rbc_datetime18')
    result_names = c('resultn_tranfus_rbc1', 'resultn_tranfus_rbc18')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw10) == time_names[i])
        e_t_ind = which(colnames(jw10) == time_names[i+1])
        s_ind = which(colnames(jw10) == result_names[i])
        e_ind = which(colnames(jw10) == result_names[i+1])
        
        info_result = cbind(t(jw10[jw10$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw10[jw10$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ] # NOTE: SWAP rows 299 and 298
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    # jw14 & jw14 ----------------------------------------------------------------
    time_names = c('lactate_datetime1', 'lactate_datetime93')
    result_names = c('resultn_lactate1', 'resultn_lactate93')
    
    for(i in seq(1, length(time_names), 2)) {
        s_t_ind = which(colnames(jw14) == time_names[i])
        e_t_ind = which(colnames(jw14) == time_names[i+1])
        s_ind = which(colnames(jw14) == result_names[i])
        e_ind = which(colnames(jw14) == result_names[i+1])
        
        info_result = cbind(t(jw14[jw14$key == key_num, s_t_ind:e_t_ind]), 
                            t(jw14[jw14$key == key_num, s_ind:e_ind]))
        info_result = data.frame(info_result)
        colnames(info_result) = c('time', result_names[i])
        rownames(info_result) = NULL
        info_result = info_result[order(info_result$time), ]
        
        sub_long_list[[main_ind]] = info_result
        main_ind = main_ind + 1
    }
    
    
    TIME = c()
    column_names = c('TIME')
    for(i in 1:length(sub_long_list)) {
        info_ind = which(!is.na(sub_long_list[[i]][,1]))
        TIME = c(TIME, sub_long_list[[i]][info_ind,1])
        
        column_names = c(column_names, colnames(sub_long_list[[i]])[2])
    }
    TIME = sort(unique(TIME))
    
    main_df = matrix(nrow = length(TIME), ncol = length(sub_long_list) + 1)
    colnames(main_df) = column_names
    main_df[,1] = TIME
    for(i in 1:length(sub_long_list)) {
        ind_1 = which(TIME %in% sub_long_list[[i]]$time)
        ind_2 = which(sub_long_list[[i]]$time %in% TIME) # We have repeat measures
        
        if(length(ind_1) != length(ind_2)) {
            print(paste0("Timing Issue: ", key_num))
            timing_issues = c(timing_issues, key_num)
            next
        }
        
        main_df[ind_1, i+1] = sub_long_list[[i]][ind_2, 2]
        
        
    }
    
    long_data_rd[[list_ind]] = main_df
    list_ind = list_ind + 1
}


save(long_data_rd, file = "Data/long_data_rd.rda")
save(timing_issues, file = "Data/timing_issues.rda")


# (5) Cleaning and reorganizing ------------------------------------------------
long_data_clean <- vector(mode = 'list', length = length(long_data_rd))

for(i in 1:length(long_data_rd)) {
    covariates = long_data_rd[[i]][, c("TIME", "resultn_temp1", "resultn_hemo1", 
                                       "resultn_map1", "resultn_hr1", "resultn_lactate1",
                                       "resultn_tranfus_rbc1", "resultn_creat1", 
                                       "resultn_plate1", "resultn_inr1", 
                                       "resultn_hemo1", "resultn_lactate1"), drop = F]
    key = all_keys[i]
    
    long_data_clean[[i]] <- list("key" = key, "covariates" = covariates)
}

for (i in 1:length(long_data_clean)) {
    # Change the time variable into minutes
    if(nrow(long_data_clean[[i]]$covariates) != 0) {
        temp = long_data_clean[[i]]$covariates[, "TIME"]
        long_data_clean[[i]]$covariates[, "TIME"] = temp / 60
    }
    
    # Filter out any extreme values
    # print(paste0(i, ": ", sum(long_data_clean[[i]]$covariates[,"resultn_map1"] > 150, na.rm = T)))
    # long_data_clean[[i]]$covariates[which(long_data_clean[[i]]$covariates[,"resultn_map1"] > 150), 
    #                                 "resultn_map1"] = NA
    # print(paste0(i, ": ", sum(long_data_clean[[i]]$covariates[,"resultn_map1"] < 25, na.rm = T)))
    # long_data_clean[[i]]$covariates[which(long_data_clean[[i]]$covariates[,"resultn_map1"] < 25), 
    #                                 "resultn_map1"] = NA
    # print(paste0(i, ": ", sum(long_data_clean[[i]]$covariates[,"resultn_hr1"] > 200, na.rm = T)))
    # long_data_clean[[i]]$covariates[which(long_data_clean[[i]]$covariates[,"resultn_hr1"] > 200), 
    #                                 "resultn_hr1"] = NA
    # print(paste0(i, ": ", sum(long_data_clean[[i]]$covariates[,"resultn_hr1"] < 20, na.rm = T)))
    # long_data_clean[[i]]$covariates[which(long_data_clean[[i]]$covariates[,"resultn_hr1"] < 20), 
    #                                 "resultn_hr1"] = NA
    # print(paste0(i, ": ", sum(long_data_clean[[i]]$covariates[,"resultn_hemo1"] > 20, na.rm = T)))
    # long_data_clean[[i]]$covariates[which(long_data_clean[[i]]$covariates[,"resultn_hemo1"] > 20), 
    #                                 "resultn_hemo1"] = NA
    # print(paste0(i, ": ", sum(long_data_clean[[i]]$covariates[,"resultn_hemo1"] < 0, na.rm = T)))
    # long_data_clean[[i]]$covariates[which(long_data_clean[[i]]$covariates[,"resultn_hemo1"] < 0), 
    #                                 "resultn_hemo1"] = NA
}

# Add number of n_labs up to that time point
long_data_agg = vector(mode = 'list', length = length(long_data_rd))
for (i in 1:length(long_data_clean)) {
    print(i)
    key = long_data_clean[[i]]$key
    covariates = NULL
    
    if(nrow(long_data_clean[[i]]$covariates) > 1) {
        time_seq = seq(min(long_data_clean[[i]]$covariates[,"TIME"]),
                       max(long_data_clean[[i]]$covariates[,"TIME"]),
                       by = 15)
        
        # Note: Check if we missed the last point
        if(max(time_seq) < max(long_data_clean[[i]]$covariates[,"TIME"])) {
            time_seq = seq(min(long_data_clean[[i]]$covariates[,"TIME"]),
                           max(long_data_clean[[i]]$covariates[,"TIME"]) + 15,
                           by = 15)
        }
        
        covariates = matrix(nrow = length(time_seq), ncol = 9)
        colnames(covariates) = c(colnames(long_data_clean[[i]]$covariates)[1:7], "n_labs", "n_RBC")
        
        n_labs = sum(!is.na(long_data_clean[[i]]$covariates[1, c("resultn_creat1", "resultn_plate1", 
                                                                 "resultn_inr1", "resultn_hemo1", "resultn_lactate1")]))
        
        n_RBC = sum(!is.na(long_data_clean[[i]]$covariates[1, "resultn_tranfus_rbc1"]))
        
        
        covariates[1, ] = c(long_data_clean[[i]]$covariates[1, 1:7], n_labs, n_RBC)
        covariates[ , "TIME"] = time_seq
        
        for(j in 2:length(time_seq)) {
            
            indices = which(long_data_clean[[i]]$covariates[, "TIME"] <= time_seq[j] &
                                long_data_clean[[i]]$covariates[, "TIME"] > time_seq[j-1])
            if(length(indices) != 0) {
                
                # temp -----------------------------------------------------------------
                temp = long_data_clean[[i]]$covariates[indices, "resultn_temp1"]
                if(sum(is.na(temp)) < length(temp)) {
                    temp = mean(temp, na.rm = T)
                } else {temp = NA}
                
                # hemo -----------------------------------------------------------------
                hemo = long_data_clean[[i]]$covariates[indices, "resultn_hemo1"]
                if(sum(is.na(hemo)) < length(hemo)) {
                    hemo = mean(hemo, na.rm = T)
                } else {hemo = NA}
                
                # map  -----------------------------------------------------------------
                map = long_data_clean[[i]]$covariates[indices, "resultn_map1"]
                if(sum(is.na(map)) < length(map)) {
                    map = mean(map, na.rm = T)
                } else {map = NA}
                
                # hr   -----------------------------------------------------------------
                hr = long_data_clean[[i]]$covariates[indices, "resultn_hr1"]
                if(sum(is.na(hr)) < length(hr)) {
                    hr = mean(hr, na.rm = T)
                } else {hr = NA}
                
                # lactate --------------------------------------------------------------
                lact = long_data_clean[[i]]$covariates[indices, "resultn_lactate1"]
                if(sum(is.na(lact)) < length(lact)) {
                    lact = mean(lact, na.rm = T)
                } else {lact = NA}
                
                # rbc  -----------------------------------------------------------------
                rbc = long_data_clean[[i]]$covariates[indices, "resultn_tranfus_rbc1"]
                if(sum(is.na(rbc)) < length(rbc)) {
                    rbc = mean(rbc, na.rm = T)
                } else {rbc = NA}
                
                # n_labs  --------------------------------------------------------------
                n_labs = sum(!is.na(long_data_clean[[i]]$covariates[indices, c("resultn_creat1", 
                                                                               "resultn_plate1", 
                                                                               "resultn_inr1", 
                                                                               "resultn_hemo1",
                                                                               "resultn_lactate1")]))
                
                # n_RBC  ---------------------------------------------------------------
                n_RBC = sum(!is.na(long_data_clean[[i]]$covariates[indices, "resultn_tranfus_rbc1"]))
                
                
                covariates[j, 2:7] = c(temp, hemo, map, hr, lact, rbc)
                covariates[j, "n_labs"] = covariates[j-1, "n_labs"] + n_labs
                covariates[j, "n_RBC"] = covariates[j-1, "n_RBC"] + n_RBC
                
            } else {
                covariates[j, "n_labs"] = covariates[j-1, "n_labs"]
                covariates[j, "n_RBC"] = covariates[j-1, "n_RBC"]
            }
        }
        
    } else {
        if(nrow(long_data_clean[[i]]$covariates) == 0) {
            covariates = matrix(nrow = 1, ncol = 9)
            colnames(covariates) = c(colnames(long_data_clean[[i]]$covariates)[1:7], "n_labs", "n_RBC")
        } else {
            covariates = matrix(c(long_data_clean[[i]]$covariates[,1:7], rep(NA, 2)),nrow = 1, ncol = 9)
            colnames(covariates) = c(colnames(long_data_clean[[i]]$covariates)[1:7], "n_labs", "n_RBC")
        }
    }
    
    covariates = cbind(rep(key, nrow(covariates)), covariates)
    colnames(covariates) = c('EID', 'time', 'temp', 'hemo', 'map', 'hr', 'lactate', 'RBC', 'n_labs', 'n_RBC')
    
    long_data_agg[[i]] <- list("key" = key, "covariates" = covariates)
    
}

# Filter out the patients that had time conflicts
timing_issues = unique(timing_issues)
t_i_ind = which(all_keys %in% timing_issues)
for(i in 1:length(long_data_agg)) {
    if (i %in% t_i_ind) {
        print("Time issue")
        long_data_agg[[i]]$time_flag = 1
    } else {
        long_data_agg[[i]]$time_flag = 0
    }
}

save(long_data_agg, file = "Data/long_data_agg.rda")



