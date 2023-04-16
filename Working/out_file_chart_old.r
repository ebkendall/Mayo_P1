library(matrixStats)
library(plotrix)
args <- commandArgs(TRUE)
set.seed(args[1])

trialNum = 12 # CHANGE EVERY TIME ******************
itNum = 5

Dir = 'Model_out/'

load(paste0(Dir,'mcmc_out_interm_',toString(args[1]),'_', trialNum,'it',itNum,'.rda'))
mcmc_out_temp$B_chain = mcmc_out_temp$B_chain[300:1000, ]
mcmc_out_temp$hr_chain = mcmc_out_temp$hr_chain[300:1000, ]
mcmc_out_temp$bp_chain = mcmc_out_temp$bp_chain[300:1000, ]
mcmc_out_temp$hc_chain = mcmc_out_temp$hc_chain[300:1000, ]
mcmc_out_temp$la_chain = mcmc_out_temp$la_chain[300:1000, ]

simulation=F
load("Data/data_format_FULL_48hr_update_RBC_sub.rda")
pace_id = c(18075, 108825, 110750, 125025, 173750, 260100, 304700, 307225, 310100,
            382450, 429375, 516150, 533075, 666750, 677225, 732525, 763050, 767500, 
            769025, 777175, 794900, 799125, 819225)
data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
data_format[which(data_format[,"EID"] == 100950)[15:18], 'map'] = 51
data_format[which(data_format[,"EID"] == 100950)[18], 'hemo'] = 6
use_data = data_format

# load(paste0('../Data/Debug/use_data', 1, '_7.rda'))

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
pdf(paste0('Plots/chart_plot_', trialNum, '_it', itNum, '_seed',toString(args[1]), '.pdf'))
# mar=c(b,l,t,r) oma=c(b,l,t,r) 
if(simulation) panels = c(5, 1)  else  panels = c(5,1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (use_data[,'EID']==i)
	n_i = sum(indices_i)
	t_grid = seq( 0, n_i, by=5)[-1]
	rbc_times = which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0)
	rbc_admin_times = which(use_data[use_data[,'EID']==i, 'RBC_admin'] != 0)
	
	# true_state = use_data[ indices_i, 'b_true']
	b_i_check = use_data[ indices_i,'clinic_rule']
	if(simulation){
		b_i = use_data[ indices_i,'b_true']
		to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
		to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
		to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]	
	}
	# HEART RATE --------------------------------------------------------------
	if(sum(!is.na(mcmc_out_temp$hr_chain[, indices_i]))==0){
		plot.new()
	} else{ 
		hr_upper = colQuantiles( mcmc_out_temp$hr_chain[, indices_i, drop=F], probs=.975)
		hr_lower = colQuantiles( mcmc_out_temp$hr_chain[, indices_i, drop=F], probs=.025)
	
		plotCI( x=colMeans(mcmc_out_temp$hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
				main=paste0('Bleeding Clinic Rule: ', b_i_check[1], '\n',
				'heart rate: ', i), xlab='time', ylab=NA, xaxt='n', col.main='green',
						col.axis='green', pch=20, cex=1, sfrac=.0025)
		grid( nx=20, NULL, col='white')
		axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
		abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
		abline(v = rbc_admin_times, col = 'deeppink', lwd = 1)
	}
	if(simulation){
		abline( v=to_s1, col='dodgerblue', lwd=2)
		abline( v=to_s2, col='firebrick1', lwd=2)
		abline( v=to_s3, col='yellow2', lwd=2)
		col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
		abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
	}
	
	# MAP --------------------------------------------------------------
	if(sum(!is.na(mcmc_out_temp$bp_chain[, indices_i]))==0){
		plot.new()
	} else{ 
		bp_upper = colQuantiles( mcmc_out_temp$bp_chain[, indices_i, drop=F], probs=.975)
		bp_lower = colQuantiles( mcmc_out_temp$bp_chain[, indices_i, drop=F], probs=.025)
		plotCI( x=colMeans(mcmc_out_temp$bp_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
				main=paste0('mean arterial pressure: ', i), xlab='time', ylab=NA, xaxt='n', 
						col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025)
		grid( nx=20, NULL, col='white')
		axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
		abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
		abline(v = rbc_admin_times, col = 'deeppink', lwd = 1)
	}
	if(simulation){
		abline( v=to_s1, col='dodgerblue', lwd=2)
		abline( v=to_s2, col='firebrick1', lwd=2)
		abline( v=to_s3, col='yellow2', lwd=2)
		col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
		abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
	}

	# HEMO --------------------------------------------------------------
	if(sum(!is.na(mcmc_out_temp$hc_chain[, indices_i]))==0){
		plot.new()
	} else{ 
		hc_upper = colQuantiles( mcmc_out_temp$hc_chain[, indices_i, drop=F], probs=.975)
		hc_lower = colQuantiles( mcmc_out_temp$hc_chain[, indices_i, drop=F], probs=.025)
		plotCI( x=colMeans(mcmc_out_temp$hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
				main=paste0('hemoglobin concentration: ', i), xlab='time', ylab=NA, xaxt='n', 
						col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025)
		grid( nx=20, NULL, col='white')
		axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
		abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
		abline(v = rbc_admin_times, col = 'deeppink', lwd = 1)
	}
	if(simulation){
		abline( v=to_s1, col='dodgerblue', lwd=2)
		abline( v=to_s2, col='firebrick1', lwd=2)
		abline( v=to_s3, col='yellow2', lwd=2)
		col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
		abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
	}

	# LACTATE --------------------------------------------------------------
	if(sum(!is.na(mcmc_out_temp$la_chain[, indices_i]))==0){
		plot.new()
	} else{ 
		la_upper = colQuantiles( mcmc_out_temp$la_chain[, indices_i, drop=F], probs=.975)
		la_lower = colQuantiles( mcmc_out_temp$la_chain[, indices_i, drop=F], probs=.025)
		plotCI( x=colMeans(mcmc_out_temp$la_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
				main=paste0('lactate levels: ', i), xlab='time', ylab=NA, xaxt='n', 
						col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025)
		grid( nx=20, NULL, col='white')
		axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
		abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
		abline(v = rbc_admin_times, col = 'deeppink', lwd = 1)
	}
	if(simulation){
		abline( v=to_s1, col='dodgerblue', lwd=2)
		abline( v=to_s2, col='firebrick1', lwd=2)
		abline( v=to_s3, col='yellow2', lwd=2)
		col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
		abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
	}
	
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
	axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
	axis( side=2, at=0:1, col.axis='green')
	abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
	abline(v = rbc_admin_times, col = 'deeppink', lwd = 1)
	# if(simulation){
	# 	abline( v=to_s1, col='cyan', lwd=2)
	# 	abline( v=to_s2, col='brown4', lwd=2)
	# 	abline( v=to_s3, col='darkgoldenrod3', lwd=2)
	# 	abline( v=to_s4, col='azure2', lwd=2)
	# 	abline( v=to_s5, col='bisque', lwd=2)
	# 	col_choice = c('cyan', 'brown4', 'darkgoldenrod3', 'azure2', 'bisque')
	# 	abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
	# }
	
	# # Adding a new cumulative posterior probability plot
	# indices_i_new = which(indices_i == T)
	# cumulative_post_prob = matrix(nrow = 2, ncol = length(indices_i_new))
	# for(w in 1:length(indices_i_new)) {
	# 	start_index = indices_i_new[1]
	# 	end_index = indices_i_new[w] 
	# 	if(w - 24 > 0) start_index = indices_i_new[w - 24]
	# 	y_or_n_2 = apply(mcmc_out_temp$B_chain[, start_index:end_index, drop=F],
	# 						1, function(x) (2 %in% x))
	# 	prob_2 = mean(y_or_n_2)
		
	# 	cumulative_post_prob[,w] = c(prob_2, 1-prob_2)
	# }
	# barplot( cumulative_post_prob, 
	# 				col=c('firebrick1', 'black'), 
	# 				main=paste0('cumulative prob.'), xlab='time', xaxt='n', space=0, 
	# 				col.main='green', border=NA) 
	# grid( nx=NA, NULL, col='white')
	# legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
	# 		legend=c( 'bleeding', 'LIMBO'), pch=15, pt.cex=1.5, 
	# 				col=c('firebrick1', 'black'))
	# axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
	# axis( side=2, at=0:1, col.axis='green')
	# abline(v = rbc_times, col = 'darkorchid2', lwd = 2)
	# if(simulation){
	# 	abline( v=to_s1, col='cyan', lwd=3)
	# 	abline( v=to_s2, col='brown4', lwd=3)
	# 	abline( v=to_s3, col='darkgoldenrod3', lwd=3)
	# 	col_choice = c('cyan', 'brown4', 'darkgoldenrod3')
	# 	abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
	# }
}
dev.off()
