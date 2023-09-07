

d.mat = matrix(NA, nrow = length(vacs), ncol = length(treats))
rownames(d.mat) = sort(vacs, decreasing = TRUE)
colnames(d.mat) = treats

for (i in 1:length(treats)) {
	for (j in 1:length(vacs)) {
		
		d.mat[j, i] = 
		dat$year_death[which(dat$treat == treats[i] & dat$vac == vacs[length(vacs)-j+1])]
		
	}
}


rotate = function(x) t(apply(x, 2, rev))


a = contour(x = treats, y = vacs, z = rotate(d.mat), 
	lwd = 2, col = 'grey60', labcex = 1, 
	levels = c(10, 50), axes = FALSE, frame.plot = FALSE, add = TRUE, 
	method = 'edge')
	
a = contourLines(z = rotate(d.mat), levels = c(10, 50))
	
	