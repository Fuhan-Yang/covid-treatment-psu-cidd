

low.col = '#FFFF7F'
hi.col = '#D81111'

bk.s = 10
im.dur = 365

pes = TRUE

load(paste('vac_treat_Summaries_sim', im.dur, '_20250301.Rdata', sep = ''))

file.name = paste('Heatmap_Pessimistic_', im.dur, '_by', bk.s, sep = '')

neu = all_param_combs_list$pes

neu$treat = round(neu$treat, 2)
neu$vac = round(neu$vac, 2)

block.size = bk.s / 100

dat = neu[which(neu$treat %in% round(seq(0.05, 0.95, by = block.size), 2) & 
	neu$vac %in% round(seq(0.25, 0.95, by = block.size), 2)), ]
dat$year_death = dat$cumu_death / 8 / 1000

treats = sort(unique(dat$treat))
vacs = sort(unique(dat$vac))

cols = colorRamp(c(low.col, hi.col))

pdf(paste(file.name, '.pdf', sep = ''), width = 9.5, height = 6.5)
par(mar = c(3, 3, 1, 0))
layout(cbind(matrix(1, nrow = 8, ncol = 8), matrix(2, nrow = 8, ncol = 2)))


plot(0, 0, col = 'white', xlim = c(ifelse(bk.s == 10, 0, 0.025), ifelse(bk.s == 10, 0.97, 0.96)), 
	ylim = c(min(vacs) - (block.size/2), max(vacs) + (block.size/2)), 
	xlab = ' ', ylab = ' ', 
	bty = 'n', xaxt = 'n', yaxt = 'n')
title(xlab = 'Treatment Coverage (%)', ylab = 'Vaccine Coverage (%)', cex.lab = 1.7, line = 1.5)
axis(1, at = treats, labels = 100 * treats, 
	line = -1.7, cex.axis = 1.4)
axis(2, at = vacs, labels = 100 * vacs, 
	line = -2, cex.axis = 1.4)
col.vals = NULL


for (i in 1:nrow(dat)) {
	
	x = dat$treat[i] + (0.5 * c(-block.size, block.size))
	x = c(x, rev(x))
	y = dat$vac[i] + (0.5 * c(-block.size, block.size))
	y = rep(y, c(2, 2))
	
	box.col = as.numeric(cols(dat$year_death[i] / ceiling(max(dat$year_death)))) / 255
	box.col = as.numeric(cols(min(dat$year_death[i] / 200, 1))) / 255
	
	polygon(x, y, col = rgb(box.col[1], box.col[2], box.col[3], 1), 
		border = rgb(box.col[1], box.col[2], box.col[3], 1))
	t.val = ifelse(dat$year_death[i] < 9.95, 
		format(round(dat$year_death[i], digits = 1), nsmall = 1), round(dat$year_death[i], digits = 0))
	text(dat$treat[i], dat$vac[i], t.val, cex = 1.2+(5*block.size))

	
	col.vals = c(col.vals, dat$year_death[i] / ceiling(max(dat$year_death)))
	
	
}
polygon(c(min(treats)-(block.size/2), max(treats)+(block.size/2), max(treats)+(block.size/2), min(treats)-(block.size/2)), 
	c(max(vacs)+(block.size/2), max(vacs)+(block.size/2), min(vacs)-(block.size/2), min(vacs)-(block.size/2)), 
	border = 'black')


if (bk.s == 5) {
	polygon(c(0.15-(block.size/2), 0.15+(block.size/2), 0.15+(block.size/2), 0.15-(block.size/2)), 
		c(0.50+(block.size/2), 0.50+(block.size/2), 0.50-(block.size/2), 0.50-(block.size/2)), 
		border = 'gray30')
} else if (bk.s == 10) {
	polygon(c(0.15-(block.size/2), 0.15+(block.size/2), 0.15+(block.size/2), 0.15-(block.size/2)), 
		c(0.45+(block.size/2), 0.45+(block.size/2), 0.45-(block.size/2), 0.45-(block.size/2)), 
		border = 'gray30')
}



plot.new()
par(mar = c(3, 0.1, 2, 0.1))
num.cols = 13
lgd = rep(NA, num.cols)
lgd[seq(1, num.cols, len = 5)] = round(seq(round(ceiling(max(dat$year_death)), -1), 0, len = 5))
lgd[seq(1, num.cols, len = 5)] = round(seq(200, 0, len = 5))
if (pes == TRUE) { lgd[1] = paste(lgd[1], '+') }
legend(x = 0, y = 0.69, legend = lgd, 
	fill = colorRampPalette(colors = c('white', 'white'))(num.cols), 
	border = NA, y.intersp = 0.75, cex = 2, text.font = NULL, bty = 'n')
legend(x = -0.1, y = 0.7, legend = rep(NA, length(lgd)), text.col = 'white', 
	fill = colorRampPalette(colors = c(hi.col, low.col))(num.cols), 
	border = NA, y.intersp = 0.5, cex = 3, text.font = NULL, bty = 'n')
text(0.38, 0.79, 'Annual\ndeaths,\n2025-2033\n(thousands)', cex = 1.9)


dev.off()








load('treat_cov_baseline_by_state.Rdata')
load('vac_cov_baseline_by_state.Rdata')

cov.state = merge(pct_vacc_by_state[ , c('Location', 'opt_pct_vacc')], treat_cov_by_state_df[ , c('state', 'cov')], 
	by.x = 'Location', by.y = 'state')
names(cov.state) = c('state', 'vac', 'treat')

source('State_pos_shift2.R')

cov.state$vac_pos = (cov.state$vac) + 0
cov.state$treat_pos = (cov.state$treat/100) + 0

cov.state$vac_pos_txt = (cov.state$vac) + (vac.shift2/100)
cov.state$treat_pos_txt = (cov.state$treat/100) + (treat.shift2/100)



pdf(paste(file.name, '_state.pdf', sep = ''), width = 9, height = 6.5)
par(mar = c(3, 3, 1, 0))
layout(cbind(matrix(1, nrow = 8, ncol = 8), matrix(2, nrow = 8, ncol = 2)))

if (bk.s == 5) {
	xl = c(0.0322, 0.42)
	yl = c(0.235, 0.82)
} else if (bk.s == 10) {
	xl = c(0.02, 0.97)
	yl = c(0.215, 0.99)
}

plot(0, 0, col = 'white', xlim = xl, # c(min(treats) - (block.size/2), max(treats) + (block.size/2)), 
	ylim = yl, # c(min(vacs) - (block.size/2), max(vacs) + (block.size/2)), 
	xlab = ' ', ylab = ' ', 
	bty = 'n', xaxt = 'n', yaxt = 'n')
title(xlab = 'Treatment Coverage (%)', ylab = 'Vaccine Coverage (%)', cex.lab = 1.4, line = 1.7)
axis(1, at = treats, labels = 100 * treats, 
	line = -1, cex.axis = 1.25)
axis(2, at = vacs, labels = 100 * vacs, 
	line = -1, cex.axis = 1.25)
col.vals = NULL

for (i in 1:nrow(dat)) {
	
	x = dat$treat[i] + (0.5 * c(-block.size, block.size))
	x = c(x, rev(x))
	y = dat$vac[i] + (0.5 * c(-block.size, block.size))
	y = rep(y, c(2, 2))
	
	box.col = as.numeric(cols(dat$year_death[i] / ceiling(max(dat$year_death)))) / 255
	box.col = as.numeric(cols(min(dat$year_death[i] / 200, 1))) / 255
	
	polygon(x, y, col = rgb(box.col[1], box.col[2], box.col[3], 1), 
		border = rgb(box.col[1], box.col[2], box.col[3], 1))

	col.vals = c(col.vals, dat$year_death[i] / ceiling(max(dat$year_death)))
	
	
}
	
for (qq in 1:nrow(cov.state)) {
	if (add.lines[qq] == 1) {
		lines(c(cov.state$treat_pos[qq], cov.state$treat_pos[qq]+(0.75*treat.shift2[qq]/100)), 
			c(cov.state$vac_pos[qq], cov.state$vac_pos[qq]+(0.75*vac.shift2[qq]/100)), 
			col = 'grey50')
	}
}	

source('add_contour.R')


points(cov.state$treat_pos, cov.state$vac_pos, pch = pchs, col = all.col)
text(cov.state$treat_pos_txt, cov.state$vac_pos_txt, cov.state$state, 
	col = all.col, cex = 1.2)
polygon(c(min(treats)-(block.size/2), max(treats)+(block.size/2), max(treats)+(block.size/2), min(treats)-(block.size/2)), 
	c(max(vacs)+(block.size/2), max(vacs)+(block.size/2), min(vacs)-(block.size/2), min(vacs)-(block.size/2)), 
	border = 'black')



plot.new()
par(mar = c(3, 0.1, 2, 0.1))
num.cols = 13
lgd = rep(NA, num.cols)
lgd[seq(1, num.cols, len = 5)] = round(seq(round(ceiling(max(dat$year_death)), -1), 0, len = 5))
lgd[seq(1, num.cols, len = 5)] = round(seq(200, 0, len = 5))
if (pes == TRUE) { lgd[1] = paste(lgd[1], '+') }
legend(x = 0, y = 0.69, legend = lgd, 
	fill = colorRampPalette(colors = c('white', 'white'))(num.cols), 
	border = NA, y.intersp = 0.75, cex = 2, text.font = NULL, bty = 'n')
legend(x = -0.1, y = 0.7, legend = rep(NA, length(lgd)), text.col = 'white', 
	fill = colorRampPalette(colors = c(hi.col, low.col))(num.cols), 
	border = NA, y.intersp = 0.5, cex = 3, text.font = NULL, bty = 'n')
text(0.38, 0.79, 'Annual\ndeaths,\n2025-2033\n(thousands)', cex = 1.9)


dev.off()








