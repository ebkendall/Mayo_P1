set.seed(2024)

# Loading the old data to keep the same IDs -----------------
load(paste0('../Data/data_format_new', 3, '.rda'))
old_EIDs = unique(data_format[,"EID"])
old_time = data_format[,"time"]
old_ids = data_format[,"EID"]

# Data with 7,488 patients to choose from ------------------
load('Data_updates/data_format.rda')
all_EIDs = unique(data_format[,'EID'])
save(all_EIDs, file = 'Data_updates/all_EIDs.rda')


data_format_big = data_format


# Splitting the data into 10 separate datasets of size 748
# (a) Old EIDs 
EID_a = all_EIDs[all_EIDs %in% old_EIDs]
EID_a = EID_a[-c(1:5)]

data_format = data_format_big[data_format_big[,"EID"] %in% EID_a, ]
save(data_format, file = 'Data_updates/data_format_1.rda')

used_EIDs = EID_a
for(i in 2:10) {
    possible_EIDs = all_EIDs[!(all_EIDs %in% used_EIDs)]
    
    EID_i = sample(possible_EIDs, size = 748, replace = F)
    
    data_format = data_format_big[data_format_big[,"EID"] %in% EID_i, ]
    save(data_format, file = paste0('Data_updates/data_format_', i, '.rda'))
    
    used_EIDs = c(used_EIDs, EID_i)
}


# Double checking there is no overlap
all_EIDs = NULL
for(i in 1:10) {
    load(paste0('Data_updates/data_format_', i, '.rda'))
    all_EIDs = c(all_EIDs, unique(data_format[,"EID"]))
}
length(all_EIDs)
length(unique(all_EIDs))
