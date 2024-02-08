library(matrixStats)
library(plotrix)

data_num = 6
simulation = T
all_seeds = F

if(simulation) {
    trialNum = 2
    itNum = 5
} else {
    trialNum = 7
    itNum = 5
}


# Load the model output -------------------------------------------------------
B_chain   = NULL
Hr_chain  = NULL
Map_chain = NULL
Hc_chain  = NULL
La_chain  = NULL

Dir = 'Model_out/'

if(all_seeds) {
    seed_list = 1:3
} else {
    seed_list = 1
}

for(seed_num in 1:length(seed_list)) {
    print(seed_num)
    if(simulation) {
        load(paste0(Dir,'mcmc_out_interm_', seed_list[seed_num],'_', trialNum,'it',itNum,'_sim.rda'))   
    } else {
        load(paste0(Dir,'mcmc_out_interm_', seed_list[seed_num],'_', trialNum,'it',itNum,'.rda'))
    }
    
    if(seed_num == 1) {
        B_chain   = mcmc_out_temp$B_chain
        Hr_chain  = mcmc_out_temp$hr_chain
        Map_chain = mcmc_out_temp$bp_chain
        Hc_chain  = mcmc_out_temp$hc_chain
        La_chain  = mcmc_out_temp$la_chain
    } else {
        B_chain   = rbind(B_chain  , mcmc_out_temp$B_chain)
        Hr_chain  = rbind(Hr_chain , mcmc_out_temp$hr_chain)
        Map_chain = rbind(Map_chain, mcmc_out_temp$bp_chain)
        Hc_chain  = rbind(Hc_chain , mcmc_out_temp$hc_chain)
        La_chain  = rbind(La_chain , mcmc_out_temp$la_chain)
    }
}
#  ----------------------------------------------------------------------------

if(simulation) {
    load(paste0('Data/use_data1_', data_num, '.rda'))
} else {
    load('Data/data_format_new3.rda')
    use_data = data_format   
}

EIDs = unique(use_data[,'EID'])
load('Data/med_select_FINAL3.rda')

# ------------------------------------------------------------------------------
# Patient Charts ---------------------------------------------------------------
# ------------------------------------------------------------------------------
pdf_title = NULL
if(all_seeds) {
    if(simulation) {
        pdf_title = paste0('Plots/chart_plot_', trialNum, '_it', itNum, '_sim.pdf')
    } else {
        pdf_title = paste0('Plots/chart_plot_', trialNum, '_it', itNum, '.pdf')
    }   
} else {
    if(simulation) {
        pdf_title = paste0('Plots/chart_plot_', trialNum, '_it', itNum, '_seed',seed_list, '_sim.pdf')
    } else {
        pdf_title = paste0('Plots/chart_plot_', trialNum, '_it', itNum, '_seed',seed_list, '.pdf')
    }
}
pdf(pdf_title)
par(mfrow=c(5,1), mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    print(which(EIDs == i))
    indices_i = (use_data[,'EID']==i)
    n_i = sum(indices_i)
    
    med_i = med_select_FINAL[med_select_FINAL$key == i, ]
    
    t_grid = round(use_data[indices_i, 'time'] / 60, digits = 3)
    t_grid_bar = 1:length(t_grid)
    rbc_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0)
    if(simulation) {
        rbc_admin = c(head(use_data[use_data[,'EID']==i, "n_RBC_admin"], 1),
                      diff(use_data[use_data[,'EID']==i, "n_RBC_admin"]))
        rbc_admin_times_bar = which(rbc_admin != 0)
    } else {
        rbc_admin_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_admin'] != 0)   
    }
    rbc_times = t_grid[rbc_times_bar]
    rbc_admin_times = t_grid[rbc_admin_times_bar]
    
    if(simulation) {
        # Put this on the correct scale as the t_grid
        b_i = use_data[ indices_i,'b_true']
        to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    } else {
        
    }

    pb = barplot(rbind(colMeans(B_chain[, indices_i] == 1),
					   colMeans(B_chain[, indices_i] == 2),
					   colMeans(B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
    
    # HEART RATE --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'hr']))==0){
        plot.new()
    } else{ 	
        if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('heart rate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('heart rate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule']))
        }
        if(simulation) {
            plot(x = pb, y = use_data[indices_i, 'hr'], main=title_name, xlab='time', ylab=NA, 
                 xaxt='n', col.main='green',
                 col.axis='green', pch=20, cex=1,
                 xlim = range(pb) + c(-0.5,0.5))
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
            
            abline( v=to_s1-0.5, col='dodgerblue', lwd=2)
            abline( v=to_s2-0.5, col='firebrick1', lwd=2)
            abline( v=to_s3-0.5, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 2)
        } else {
            hr_upper = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.975)
            hr_lower = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.025)
            plotCI( x = pb, y=colMeans(Hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
                    main=title_name,
                    xlab='time', ylab=NA, xaxt='n', col.main='green',
                    col.axis='green', pch=20, cex=1, sfrac=.0025,
                    xlim = range(pb) + c(-0.5,0.5))   
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
        }
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    }

    # MAP --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'map']))==0){
        plot.new()
    } else{ 
        if(simulation) {
            plot(x = pb, y = use_data[indices_i, 'map'], main=paste0('map: ', i), 
                 xlab='time', ylab=NA, xaxt='n', 
                 col.main='green', col.axis='green', pch=20, cex=1, 
                 xlim = range(pb) + c(-0.5,0.5))
            
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
            
            abline( v=to_s1-0.5, col='dodgerblue', lwd=2)
            abline( v=to_s2-0.5, col='firebrick1', lwd=2)
            abline( v=to_s3-0.5, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 2)
        } else {
            bp_upper = colQuantiles( Map_chain[, indices_i, drop=F], probs=.975)
            bp_lower = colQuantiles( Map_chain[, indices_i, drop=F], probs=.025)
            plotCI(x = pb, y = colMeans(Map_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
                   main=paste0('mean arterial pressure: ', i), xlab='time', ylab=NA, xaxt='n', 
                   col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025,
                   xlim = range(pb) + c(-0.5,0.5))
            
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
        }
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    }
    
    # HEMO --------------------------------------------------------------
    # if(sum(!is.na(use_data[indices_i, 'hemo']))==0){
    #     plot.new()
    # } else{ 
        if(simulation) {
            plot(x = pb, y = use_data[indices_i, 'hemo'], main=paste0('hemo: ', i), 
                 xlab='time', ylab=NA, xaxt='n', 
                 col.main='green', col.axis='green', pch=20, cex=1,
                 xlim = range(pb) + c(-0.5,0.5))
            
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
            
            abline( v=to_s1-0.5, col='dodgerblue', lwd=2)
            abline( v=to_s2-0.5, col='firebrick1', lwd=2)
            abline( v=to_s3-0.5, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 2)
        } else {
            hc_upper = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.975)
            hc_lower = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.025)
            plotCI(x = pb, y = colMeans(Hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
                   main=paste0('hemoglobin concentration: ', i), xlab='time', ylab=NA, xaxt='n', 
                   col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025,
                   xlim = range(pb) + c(-0.5,0.5))
            
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
        }
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    # }

    # LACTATE --------------------------------------------------------------
    # if(sum(!is.na(use_data[indices_i, 'lactate']))==0){
    #     plot.new()
    # } else{ 
        if(simulation) {
            plot(x = pb, y = use_data[indices_i, 'lactate'], main=paste0('lactate: ', i), 
                 xlab='time', ylab=NA, xaxt='n', 
                 col.main='green', col.axis='green', pch=20, cex=1,
                 xlim = range(pb) + c(-0.5,0.5))
            
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
            
            abline( v=to_s1-0.5, col='dodgerblue', lwd=2)
            abline( v=to_s2-0.5, col='firebrick1', lwd=2)
            abline( v=to_s3-0.5, col='yellow2', lwd=2)
            col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
            abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 2)
        } else {
            la_upper = colQuantiles( La_chain[, indices_i, drop=F], probs=.975)
            la_lower = colQuantiles( La_chain[, indices_i, drop=F], probs=.025)
            plotCI(x = pb, y = colMeans(La_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
                   main=paste0('lactate levels: ', i), xlab='time', ylab=NA, xaxt='n', 
                   col.main='green', col.axis='green', pch=20, cex=1, sfrac=.0025,
                   xlim = range(pb) + c(-0.5,0.5))
            
            grid( nx=20, NULL, col='white')
            axis( side=1, at=pb, col.axis='green', labels=t_grid)
        }
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    # }
	
    # BAR PLOTS --------------------------------------------------------------
	barplot(rbind(colMeans(B_chain[, indices_i] == 1),
				  colMeans(B_chain[, indices_i] == 2),
				  colMeans(B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main='green', border=NA,
            xlim=range(pb) + c(-0.5,0.5)) 
	grid( nx=20, NULL, col='white')
	legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
			legend=c( 'Baseline', 'Bleed', 'Recov(B)'), pch=15, pt.cex=1.5, 
					col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
	legend( 'topleft', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
			legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
					col=c( 'darkorchid1', 'aquamarine'))				
	axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
	axis( side=2, at=0:1, col.axis='green')
	
	
	abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
	abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
}
dev.off()


# ------------------------------------------------------------------------------ 
# Model evaluation plots -------------------------------------------------------
# ------------------------------------------------------------------------------
pdf_title = NULL
if(all_seeds) {
    if(simulation) {
        pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, '_sim.pdf')
    } else {
        pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, '.pdf')
    }   
} else {
    if(simulation) {
        pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, '_seed',seed_list, '_sim.pdf')
    } else {
        pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, '_seed',seed_list, '.pdf')
    }
}
pdf(pdf_title)
par(mfrow=c(5,1), mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    print(which(EIDs == i))
    indices_i = (use_data[,'EID']==i)
    n_i = sum(indices_i)
    
    med_i = med_select_FINAL[med_select_FINAL$key == i, ]
    
    t_grid = round(use_data[indices_i, 'time'] / 60, digits = 3)
    t_grid_bar = 1:length(t_grid)
    rbc_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0)
    if(simulation) {
        rbc_admin = c(head(use_data[use_data[,'EID']==i, "n_RBC_admin"], 1),
                      diff(use_data[use_data[,'EID']==i, "n_RBC_admin"]))
        rbc_admin_times_bar = which(rbc_admin != 0)
    } else {
        rbc_admin_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_admin'] != 0)   
    }
    rbc_times = t_grid[rbc_times_bar]
    rbc_admin_times = t_grid[rbc_admin_times_bar]
    
    if(simulation) {
        # Put this on the correct scale as the t_grid
        b_i = use_data[ indices_i,'b_true']
        to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    } 
    
    pb = barplot(rbind(colMeans(B_chain[, indices_i] == 1),
                       colMeans(B_chain[, indices_i] == 2),
                       colMeans(B_chain[, indices_i] == 3)), 
                 col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
                 xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
    
    # Heart Rate and MAP double plot -----------------------------------------
    if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
        title_name = paste0('heart rate & MAP: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule']),
                            ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
    } else {
        title_name = paste0('heart rate & MAP: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule']))
    }
    
    hr_upper = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.975)
    hr_lower = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.025)
    bp_upper = colQuantiles( Map_chain[, indices_i, drop=F], probs=.975)
    bp_lower = colQuantiles( Map_chain[, indices_i, drop=F], probs=.025)
    
    hr_map_ylim = c(min(hr_lower, bp_lower), max(hr_upper, bp_upper))
    
    plotCI( x = pb, y=colMeans(Hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', col.main='green',
            col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
            xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim) 
    plotCI( x = pb, y=colMeans(Map_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
            col = 'orange',
            xlim = range(pb) + c(-0.5,0.5), add = T) 
    legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
    grid( nx=20, NULL, col='white')
    axis( side=1, at=pb, col.axis='green', labels=t_grid)
    if(simulation) {
        abline( v=to_s1-0.5, col='dodgerblue', lwd=2)
        abline( v=to_s2-0.5, col='firebrick1', lwd=2)
        abline( v=to_s3-0.5, col='yellow2', lwd=2)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 2)
    } 
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    
    
    # Hemoglobin and Lactate double plot -------------------------------------
    if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
        title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule']),
                            ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
    } else {
        title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule']))
    }
    
    hc_upper = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.975)
    hc_lower = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.025)
    la_upper = colQuantiles( La_chain[, indices_i, drop=F], probs=.975)
    la_lower = colQuantiles( La_chain[, indices_i, drop=F], probs=.025)
    
    hr_map_ylim = c(min(hc_lower, la_lower), max(hc_upper, la_upper))
   
    plotCI(x = pb, y = colMeans(Hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', col.main='green',
            col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
            xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim) 
    plotCI( x = pb, y=colMeans(La_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
            col = 'orange',
            xlim = range(pb) + c(-0.5,0.5), add = T) 
    legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'hemo', 'lactate'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
    grid( nx=20, NULL, col='white')
    axis( side=1, at=pb, col.axis='green', labels=t_grid)
    if(simulation) {
        abline( v=to_s1-0.5, col='dodgerblue', lwd=2)
        abline( v=to_s2-0.5, col='firebrick1', lwd=2)
        abline( v=to_s3-0.5, col='yellow2', lwd=2)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 2)
    } 
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    
    # BAR PLOTS --------------------------------------------------------------
    barplot(rbind(colMeans(B_chain[, indices_i] == 1),
                  colMeans(B_chain[, indices_i] == 2),
                  colMeans(B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
            xlab='time', space=0, col.main='green', border=NA,
            xlim=range(pb) + c(-0.5,0.5)) 
    grid( nx=20, NULL, col='white')
    legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'Baseline', 'Bleed', 'Recov(B)'), pch=15, pt.cex=1.5, 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
    legend( 'topleft', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
            col=c( 'darkorchid1', 'aquamarine'))				
    axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
    axis( side=2, at=0:1, col.axis='green')
    
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    
    # Cumulative PLOTS ---------------------------------------------------------
    cumulative_post_prob = matrix(nrow = 2, ncol = n_i)
    ind = 1
    win_length = 1
    c = 0.18
    
    indices_i_new = which(indices_i == T)
    for(w in 1:length(indices_i_new)) {
        start_index = indices_i_new[1]
        end_index = indices_i_new[w] 
        if(w - win_length > 0) start_index = indices_i_new[w - win_length]
        
        y_or_n_2 = apply(B_chain[, start_index:end_index, drop=F],
                         1, function(x) (2 %in% x))
        prob_2 = mean(y_or_n_2)
        
        cumulative_post_prob[, ind] = c(prob_2, 1-prob_2)
        ind = ind + 1
    }
    
    barplot( cumulative_post_prob,
             col=c('darkred', 'black'),
             main=paste0('cumulative prob.'), xlab='time', space=0, 
             col.main='green', border=NA,
             xlim=range(pb) + c(-0.5,0.5))
    grid( nx=20, NULL, col='white')
    legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'bleeding', 'LIMBO'), pch=15, pt.cex=1.5,
            col=c('darkred', 'black'))
    axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
    axis( side=2, at=0:1, col.axis='green')
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    abline(h = c, col = 'yellow', lwd = 2)
    
    if(simulation) {
        abline( v=to_s1-0.5, col='dodgerblue', lwd=1)
        abline( v=to_s2-0.5, col='firebrick1', lwd=1)
        abline( v=to_s3-0.5, col='yellow2', lwd=1)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 1)
    }
    
    # State verification  ------------------------------------------------------
    bleed_or_no = as.numeric(cumulative_post_prob[1,] > c)
    plot(x=pb, y=bleed_or_no, type = 's', lwd = 2, main = 'Identification of bleeding',
         xlab='time', ylab = ' ', col.main='green', col.lab = 'green',
         xlim = range(pb) + c(-0.5,0.5),
         xaxt='n', yaxt='n', ylim = c(-1,2), col = 'white')
    axis( side=1, at=pb, col.axis='green', labels=t_grid)
    axis( side=2, at=0:1, col.axis='green', labels = c("S1/S3", "S2"),
          cex.axis=1)
    
    if(simulation) {
        correct_choice = bleed_or_no
        for(b in 1:length(b_i)) {
            if(bleed_or_no[b] == 1) {
                if(b_i[b] == 2) {
                    correct_choice[b] = 'green'
                } else {
                    correct_choice[b] = 'red'
                }
            } else {
                if(b_i[b] == 2) {
                    correct_choice[b] = 'red'
                } else {
                    correct_choice[b] = 'green'
                }
            }
        }
        
        abline( v=to_s1-0.5, col='dodgerblue', lwd=1)
        abline( v=to_s2-0.5, col='firebrick1', lwd=1)
        abline( v=to_s3-0.5, col='yellow2', lwd=1)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1-0.5, col = col_choice[b_i[1]], lwd = 1)
        
        points(x=pb, y=bleed_or_no, col = correct_choice, pch=19)   
    }
    
}
dev.off()

