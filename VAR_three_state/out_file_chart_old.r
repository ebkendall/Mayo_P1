library(matrixStats)
library(plotrix)
args <- commandArgs(TRUE)
set.seed(args[1])

trialNum = 2
itNum = 5
data_num = 3
simulation = F

Dir = 'Model_out/'

if(simulation) {
    load(paste0(Dir,'mcmc_out_interm_',toString(args[1]),'_', trialNum,'it',itNum,'_sim.rda'))
} else {
    load(paste0(Dir,'mcmc_out_interm_',toString(args[1]),'_', trialNum,'it',itNum,'.rda'))
}
mcmc_out_temp$B_chain = mcmc_out_temp$B_chain[300:1000, ]
mcmc_out_temp$hr_chain = mcmc_out_temp$hr_chain[300:1000, ]
mcmc_out_temp$bp_chain = mcmc_out_temp$bp_chain[300:1000, ]
mcmc_out_temp$hc_chain = mcmc_out_temp$hc_chain[300:1000, ]
mcmc_out_temp$la_chain = mcmc_out_temp$la_chain[300:1000, ]

if(simulation) {
    load(paste0('Data/use_data1_', data_num, '.rda'))
} else {
    load('Data/data_format_new.rda')
    pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
                417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
                732525, 758025, 763050, 843000, 117525)
    data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
    use_data = data_format   
}

load('Data/clean_meds.rda')
# load(paste0('../Data/Debug/use_data', 1, '_7.rda'))
load('Data/care_time_df.rda')
care_time_df = cbind(care_time_df, 'green')
care_time_df[care_time_df[,3] == "Intensive Care", 4] = 'red'

level_of_care = read.csv('../Data_clean/Data/_raw_data_new/jw_patient_level_of_care.csv')

EIDs = unique(use_data[,'EID'])

# for(i in EIDs){  
#   temp = use_data[use_data[,'EID']==as.numeric(i), 'b_true']
#   b_temp = matrix( temp, ncol = 1)
#   if(2 %in% temp) {
# 	use_data[use_data[,'EID']==as.numeric(i), 'clinic_rule'] = 1
#   } else {
# 	use_data[use_data[,'EID']==as.numeric(i), 'clinic_rule'] = 0
#   }
# }

# New patients ---------------------------------------------------------------
pdf_title = NULL
if(simulation) {
    pdf_title = paste0('Plots/chart_plot_', trialNum, '_it', itNum, '_seed',toString(args[1]), '_sim.pdf')
} else {
    pdf_title = paste0('Plots/chart_plot_', trialNum, '_it', itNum, '_seed',toString(args[1]), '.pdf')
}
pdf(pdf_title)
# mar=c(b,l,t,r) oma=c(b,l,t,r) 
if(simulation) panels = c(5, 1)  else  panels = c(5,1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    print(which(EIDs == i))
    indices_i = (use_data[,'EID']==i)
    n_i = sum(indices_i)
    t_grid_bar = seq( 0, n_i, by=5)[-1]
    if(simulation) {
        t_grid = seq( 0, n_i, by=5)[-1]   
        rbc_times = which(use_data[indices_i, 'RBC_ordered'] != 0)
        rbc_admin_times = which(diff(use_data[indices_i, 'n_RBC_admin']) > 0) + 1
    } else {
        t_grid = round(use_data[indices_i, 'time'] / 60, digits = 3)
        rbc_times = use_data[which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0), 'time']/60
        rbc_admin_times = use_data[which(use_data[use_data[,'EID']==i, 'RBC_admin'] != 0), 'time']/60
    }
    
    if(simulation){
        b_i = use_data[ indices_i,'b_true']
        to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]	
    }
    # HEART RATE --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'hr']))==0){
        plot.new()
    } else{ 	
        if(simulation) {
            plot(use_data[indices_i, 'hr'], main=paste0('heart rate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule'])), 
                 xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
            abline( v=to_s1, col='dodgerblue', lwd=2)
            abline( v=to_s2, col='firebrick1', lwd=2)
            abline( v=to_s3, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
        } else {
            hr_upper = colQuantiles( mcmc_out_temp$hr_chain[, indices_i, drop=F], probs=.975)
            hr_lower = colQuantiles( mcmc_out_temp$hr_chain[, indices_i, drop=F], probs=.025)
            plotCI( x = t_grid, y=colMeans(mcmc_out_temp$hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
                    main=paste0('heart rate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule'])),
                    xlab='time', ylab=NA, xaxt='n', col.main='green',
                    col.axis='green', pch=20, cex=1, sfrac=.0025)   
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
        }
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    }

    # MAP --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'map']))==0){
        plot.new()
    } else{ 
        if(simulation) {
            plot(use_data[indices_i, 'map'], main=paste0('map: ', i), 
                 xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)   
            abline( v=to_s1, col='dodgerblue', lwd=2)
            abline( v=to_s2, col='firebrick1', lwd=2)
            abline( v=to_s3, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
        } else {
            bp_upper = colQuantiles( mcmc_out_temp$bp_chain[, indices_i, drop=F], probs=.975)
            bp_lower = colQuantiles( mcmc_out_temp$bp_chain[, indices_i, drop=F], probs=.025)
            plotCI(x = t_grid, y = colMeans(mcmc_out_temp$bp_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
                   main=paste0('mean arterial pressure: ', i), xlab='time', ylab=NA, xaxt='n', 
                   col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
        }
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    }
    
    # HEMO --------------------------------------------------------------
    # if(sum(!is.na(use_data[indices_i, 'hemo']))==0){
    #     plot.new()
    # } else{ 
        if(simulation) {
            plot(use_data[indices_i, 'hemo'], main=paste0('hemo: ', i), 
                 xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
            abline( v=to_s1, col='dodgerblue', lwd=2)
            abline( v=to_s2, col='firebrick1', lwd=2)
            abline( v=to_s3, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
        } else {
            hc_upper = colQuantiles( mcmc_out_temp$hc_chain[, indices_i, drop=F], probs=.975)
            hc_lower = colQuantiles( mcmc_out_temp$hc_chain[, indices_i, drop=F], probs=.025)
            plotCI(x = t_grid, y = colMeans(mcmc_out_temp$hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
                   main=paste0('hemoglobin concentration: ', i), xlab='time', ylab=NA, xaxt='n', 
                   col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
        }
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    # }
    
    # LACTATE --------------------------------------------------------------
    # if(sum(!is.na(use_data[indices_i, 'lactate']))==0){
    #     plot.new()
    # } else{ 
        if(simulation) {
            plot(use_data[indices_i, 'lactate'], main=paste0('lactate: ', i), 
                 xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
            abline( v=to_s1, col='dodgerblue', lwd=2)
            abline( v=to_s2, col='firebrick1', lwd=2)
            abline( v=to_s3, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
        } else {
            la_upper = colQuantiles( mcmc_out_temp$la_chain[, indices_i, drop=F], probs=.975)
            la_lower = colQuantiles( mcmc_out_temp$la_chain[, indices_i, drop=F], probs=.025)
            plotCI(x = t_grid, y = colMeans(mcmc_out_temp$la_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
                   main=paste0('lactate levels: ', i), xlab='time', ylab=NA, xaxt='n', 
                   col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025)
            grid( nx=20, NULL, col='white')
            axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
        }
        
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    # }
	
    # BAR PLOTS --------------------------------------------------------------
	barplot( rbind( colMeans(mcmc_out_temp$B_chain[, indices_i] == 1),
					colMeans(mcmc_out_temp$B_chain[, indices_i] == 2),
					colMeans(mcmc_out_temp$B_chain[, indices_i] == 3)), 
					col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
					xlab='time', xaxt='n', space=0, 
					col.main='green', border=NA) 
	grid( nx=NA, NULL, col='white')
	legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
			legend=c( 'Baseline', 'Bleed', 'Recov(B)'), pch=15, pt.cex=1.5, 
					col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
	legend( 'topleft', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
			legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
					col=c( 'darkorchid2', 'deeppink'))				
	axis( side=1, at=t_grid_bar, col.axis='green', labels=t_grid_bar / 4)
	axis( side=2, at=0:1, col.axis='green')
	abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
	abline(v = rbc_admin_times, col = 'grey', lwd = 1)
}
dev.off()
