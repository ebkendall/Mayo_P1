# Model evaluation to determine bleeding threshold and optimal time window.

# Adding a new cumulative posterior probability plot
indices_i_new = which(indices_i == T)
cumulative_post_prob = matrix(nrow = 2, ncol = length(indices_i_new))
for(w in 1:length(indices_i_new)) {
	start_index = indices_i_new[1]
	end_index = indices_i_new[w] 
	if(w - 24 > 0) start_index = indices_i_new[w - 24]
	y_or_n_2 = apply(mcmc_out_temp$B_chain[, start_index:end_index, drop=F],
						1, function(x) (2 %in% x))
	prob_2 = mean(y_or_n_2)

	cumulative_post_prob[,w] = c(prob_2, 1-prob_2)
}
barplot( cumulative_post_prob, 
				col=c('firebrick1', 'black'), 
				main=paste0('cumulative prob.'), xlab='time', xaxt='n', space=0, 
				col.main='green', border=NA) 
grid( nx=NA, NULL, col='white')
legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
		legend=c( 'bleeding', 'LIMBO'), pch=15, pt.cex=1.5, 
				col=c('firebrick1', 'black'))
axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
axis( side=2, at=0:1, col.axis='green')
abline(v = rbc_times, col = 'darkorchid2', lwd = 2)
if(simulation){
	abline( v=to_s1, col='cyan', lwd=3)
	abline( v=to_s2, col='brown4', lwd=3)
	abline( v=to_s3, col='darkgoldenrod3', lwd=3)
	col_choice = c('cyan', 'brown4', 'darkgoldenrod3')
	abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
}