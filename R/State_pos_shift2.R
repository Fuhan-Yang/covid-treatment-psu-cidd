

shift2 = c(
	-0.8, -0.9, #AK
	1, 0, #AL
	0, -1, #AR
	0, 1, #AZ
	1, 0, #CA
	-1, 0, #CO
	0, -1, #CT
	0, 1, #DE
	0.4, 2.8, #FL
	1, 0, #GA
	-1, 0, #HI
	1, -0.2, #IA
	-1, 0, #ID
	0.9, 0, #IL
	0, -1, #IN
	0, 1, #KS
	-1, 0, #KY
	0.3, 1, #LA
	1, 0, #MA
	0, -1, #MD
	-1, 0, #ME
	1, 0, #MI
	1, -0.3, #MN
	0.5, -1, #MO
	-0.3, -1, #MS
	-1, 0.2, #MT
	3, 1.1, #NC
	-1, 0, #ND
	-1, 0, #NE
	0, -1, #NH
	0, -1, #NJ
	1, 0, #NM
	1.8, 2, #NV
	1, 0, #NY
	1, 0, #OH
	0.7, -0.8, #OK
	0, 1.05, #OR
	1, 0, #PA
	-1, 0, #RI
	1, 0, #SC
	-1, 0, #SD
	-0.5, 0.85, #TN
	0.2, -1, #TX
	-1, 0, #UT
	1, 0, #VA
	1, 0, #VT
	-0.3, 1, #WA
	-1, 0, #WI
	1.1, 0, #WV
	-1, 0 #WY
)


vac.shift2 = shift2[seq(2, length(shift2), by = 2)] * 1.15
treat.shift2 = shift2[seq(1, length(shift2), by = 2)] * 0.9

add.lines = rep(0, 50)
add.lines[which(cov.state$state %in% 
	c('FL', 'NC', 'NV')
	)] = 1


pchs = rep(20, 50)
