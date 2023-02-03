#### After 2022-08-18, use treatment successful rate to model treatment efficacy
# 'treat-eff-x0' is no longer used #

#### 2022-10-27 version #####
### only susceptible get vaccinated, beta is calibrated based on that ###
### waning immunity duration = 365 ###
### update the RI vaccinated population, beta is calibrated ###

generate_baseline_params_RI <- function(odepath="cpp-v9-rebound//",
                                        treat_beginday = 732,
                                        treat_endday = 1065,
                                        beta_base_period_start,  
                                        beta_base_period,
                                        sim_time,
                                        beta_transition_period,
                                        beta_multiplier = 1,
                                        seas_amp,
                                        seas_dur,
                                        seas_start,
                                        beta_adopt = F,
                                        adopted_beta = NULL){

  ### mcmc initialization (don't change)
  # end.day = day index of the last date with available data with 2020-01-01 being day 1 
  ### starts from 2020-02-27 to 2022-11-30
  end.day = 1065
  num.days = 1008
  
  ## The introduction day of the first case ## 
  introday = end.day - num.days + 1
  
  # source("../../data.process.R")
  # dp <- data.process(ri.data)
  # str(dp)
  
  ### Create spline expansions
  
  ## start spline at 2020-03-01 ##
  start_daynum = 61
  
  # Cubic spline with one basis function every 7 days starting from Mar 1 2020
  bspl <- create.bspline.basis(c(start_daynum, end.day), nbasis=round(num.days/7))
  
  # create basis evaluation function matrices
  Z <- eval.basis(bspl, start_daynum:end.day)
  
  
  ### create I-spline basis functions, evaluated every day from start_daynum:end.day
  Z.rr <- iSpline(x = start_daynum:end.day, knots = c(92, 122, 183, 288, 336, 381, 410), degree = 2)
  Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
  Z.rr <- Z.rr[,c(1,3:8)]
  
  
  
  ###############################################################################
  ### 2. Starting values for parameters
  ###############################################################################
  
  # starting values for reporting rate
  rr.start <- c(0.0424058267858393,0.458615599545913,0.210846692864851,0.0146983392464192,0.249688232988845,0.00937456422422415,0.00506977941951041)
  rr.daily <- Z.rr %*% rr.start
  # plot(rr.daily)
  
  
  # starting values for odesim parameters and other parameters
  prms.start <- c(10.5774733781657,
                  0.698716865548575,
                  0.660921841557831,0.957784878956134,
                  0.52845627678712,0.256526655243724,0.306050459199176,
                  156.276171748767,344.423409006052,
                  
                  0.0122239664155569,0.0209363491765736,0.0304561389882463,0.0482261102338912,0.0780703865458093,0.147326127082464,0.285027801413735,0.314214187084188,
                  0.0122239664155569*0.46,
                  0.0209363491765736*1.01,
                  0.0304561389882463*0.79,
                  0.0482261102338912*0.44,
                  0.0780703865458093*0.49,
                  0.147326127082464*0.38,
                  0.285027801413735*0.63,
                  0.314214187084188*0.77,
                  721,
                  
                  2.55513057168369,3.05020653841806,3.02253004409028,2.58876335497109,1.84026766217073,1.15495881700266,1.10255613331181,2.25217171491254,
                  2.91812626656877,2.32151696890213,1.61945634763761,1.07940122510788,0.863041731949425,0.494255363822617,0.366734867193275,0.519274018826153,
                  3.73843817900114,2.98992836747547,1.30173265877784,0.653444145160029,0.478416814779162,0.269999187994631,0.17699128750395,0.130498436380522,
                  158.021957094757,389.924577671391,
                  0,rep(0.137,8),0,rep(0,8),
                  treat_beginday,treat_endday)
  
  names(prms.start) <- c("mean-time-vent", 
                         "tv-dev-len-hospstay",
                         "prob-icu-vent", "dev-ventdeath-mid",
                         "tv-dev-icu-frac_1","tv-dev-icu-frac_2", "tv-dev-icu-frac_3",
                         "tv-dev-icu-frac-endday_1", "tv-dev-icu-frac-endday_2",
                         
                         "tv-hosp-frac-10_1","tv-hosp-frac-20_1","tv-hosp-frac-30_1","tv-hosp-frac-40_1","tv-hosp-frac-50_1","tv-hosp-frac-60_1","tv-hosp-frac-70_1", "tv-hosp-frac-80_1",
                         "tv-hosp-frac-10_2","tv-hosp-frac-20_2","tv-hosp-frac-30_2","tv-hosp-frac-40_2","tv-hosp-frac-50_2","tv-hosp-frac-60_2","tv-hosp-frac-70_2", "tv-hosp-frac-80_2",
                         "tv-hosp-frac-endday_1",
                         
                         "tv-contact-rate-10_1", "tv-contact-rate-20_1", "tv-contact-rate-30_1", "tv-contact-rate-40_1", "tv-contact-rate-50_1", "tv-contact-rate-60_1", "tv-contact-rate-70_1", "tv-contact-rate-80_1",
                         "tv-contact-rate-10_2", "tv-contact-rate-20_2", "tv-contact-rate-30_2", "tv-contact-rate-40_2", "tv-contact-rate-50_2", "tv-contact-rate-60_2", "tv-contact-rate-70_2", "tv-contact-rate-80_2",
                         "tv-contact-rate-10_3", "tv-contact-rate-20_3", "tv-contact-rate-30_3", "tv-contact-rate-40_3", "tv-contact-rate-50_3", "tv-contact-rate-60_3", "tv-contact-rate-70_3", "tv-contact-rate-80_3",
                         "tv-contact-rate-endday_1", "tv-contact-rate-endday_2",
                         
                         "tv-treat-cov-00_1","tv-treat-cov-10_1","tv-treat-cov-20_1","tv-treat-cov-30_1","tv-treat-cov-40_1","tv-treat-cov-50_1","tv-treat-cov-60_1","tv-treat-cov-70_1","tv-treat-cov-80_1",
                         "tv-treat-cov-00_2","tv-treat-cov-10_2","tv-treat-cov-20_2","tv-treat-cov-30_2","tv-treat-cov-40_2","tv-treat-cov-50_2","tv-treat-cov-60_2","tv-treat-cov-70_2","tv-treat-cov-80_2",
                         "treat-beginday",
                         "tv-treat-cov-endday")
  
  ## names of parameters that are NOT odesim parameters, but are still to be estimated
  non.odesim.params <- NULL
  ## inputs to odesim that are to be FIXED, and not estimated
  
  ## fixed hosp-frac: arithmetic mean of combined posterior from RI 0911, CT 0916, MA 0920 runs with different hosp-fracs
  const.prms.start <- c("", 2,
                        4,
                        "5\t6",
                        721,
                        "6\t3.42",
                        721,
                        "0\t0",
                        "0\t0",
                        "0\t0",
                        "0\t0",
                        "0\t0",
                        "0\t0",
                        # "0\t0",
                        # "0\t0",
                        # "0\t0",
                        # "0.0130658513676206\t0",
                        # "0.0421950867804589\t0",
                        # "0.22715636067416\t0",
                        # "0.0130658513676206\t0.001437244",
                        # "0.0421950867804589\t0.006751214",
                        # "0.22715636067416\t0.08177629",
                        "0.0130658513676206\t0.0013",
                        "0.0421950867804589\t0.006106533",
                        "0.22715636067416\t0.07396739",
                        721,
                        "347	354	361	368	375	382	389	396	403	410	417	424	431	438	445	452	459	466	473	480	487	494	501	508	515	522	529	536	543	550	557	564	571	578	585	592	599	606	613	620	627	634	641	648	655	662	669	676	683	690	697	704	711	718	725	732	739	746	753	760	767	774	781	788	795	802	809	816	823	830	837	844	851	858	865	872	879	886	893	900	907	914	921	928	935	942	949	956	963	970	977	984	991	998	1005	1012	1019	1026	1033	1040	1047	1054	1061	1068",
                        "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	630	5202	4078	1387	686	1251	1357	780	900	814	854	648	443	301	345	292	268	251	174	133	151	100	127	93	92	97	87	82	84	71	70	81	59	64	76	218	279	309	398	424	375	235	402	343	436	356	311	258	312	289	287	228	217	124	241",
                        "0	0	0	20	43	48	148	143	184	208	174	155	208	238	174	446	867	654	1074	1788	2486	6671	2934	2674	7190	8514	4328	2256	1694	1376	1417	1305	1144	1277	1388	1576	1581	1536	1485	1319	1161	952	748	655	599	480	394	361	368	586	3022	2478	1042	616	888	1081	746	844	764	864	758	556	386	340	336	261	283	236	157	135	124	119	148	120	108	112	108	87	113	80	120	74	68	76	83	137	120	186	296	162	98	176	125	160	161	123	128	173	154	138	107	93	83	116",
                        "0	0	8	496	868	554	1462	1232	1457	1461	1249	1156	1564	1526	1080	2855	4245	2314	3457	5140	5225	12851	5790	4100	2474	3033	2672	1760	1342	1002	1021	959	912	895	1028	1163	1332	1166	1117	1221	1317	1333	859	823	684	591	564	525	546	434	551	535	538	386	408	474	497	474	503	540	512	453	288	240	195	146	149	143	113	99	93	91	110	74	82	85	82	64	65	41	64	60	70	51	64	103	216	152	252	112	51	190	163	242	280	219	229	264	252	214	170	147	104	162",
                        "0	9	15	922	1298	564	2044	1586	1677	1473	1647	1405	2167	1950	1686	4124	5084	3053	4008	5412	5316	12125	5569	3689	2185	2602	2330	1507	1203	930	944	780	671	670	789	932	1189	1103	1081	1336	1423	1378	928	790	724	640	654	490	558	382	503	445	420	334	357	388	433	372	340	439	399	294	229	202	152	117	114	100	83	80	59	61	67	84	47	76	77	73	37	45	36	39	44	38	41	103	216	152	252	112	51	190	163	242	280	219	229	264	252	214	170	147	104	162",
                        "0	5	7	723	1098	564	2179	1633	1753	1537	1830	1779	2376	1955	2037	5043	5869	3473	4034	5722	7022	10824	5338	3199	2014	2212	1949	1324	938	759	751	578	566	568	642	769	783	865	739	874	845	898	594	543	491	389	359	339	364	282	367	349	309	238	264	283	282	268	256	309	270	206	124	114	100	76	97	74	70	53	52	37	43	62	55	52	39	43	44	32	29	26	23	30	24	103	216	152	252	112	51	190	163	242	280	219	229	264	252	214	170	147	104	162",
                        "0	5	16	852	1404	799	2913	2171	2497	1995	2402	2400	3366	3567	3251	7614	10958	5599	7026	12672	10136	9298	4631	3148	1852	1939	1744	1138	801	707	613	544	537	555	609	719	805	719	690	760	717	678	459	434	401	360	312	306	352	256	351	295	285	209	269	268	248	239	255	223	225	167	112	87	64	68	50	96	108	105	90	74	97	94	100	88	62	63	47	45	35	39	52	62	36	103	216	152	252	112	51	190	163	242	280	219	229	264	252	214	170	147	104	162",
                        "0	0	9	586	1009	717	2290	1756	2270	2056	2461	2744	4460	15877	13862	14420	11247	8306	8345	7222	3341	3523	2099	1634	993	1048	850	660	441	388	361	331	317	295	349	388	451	383	318	352	372	333	276	262	220	237	198	215	247	183	207	205	172	122	131	137	146	129	114	118	120	84	58	54	42	33	24	111	179	141	111	127	111	108	125	131	59	62	50	39	23	20	64	58	49	88	137	105	173	79	33	155	155	219	244	218	227	245	234	181	144	118	75	133",
                        "0	0	0	75	144	314	927	717	1210	1575	2423	3382	6074	19110	14991	9551	3956	1042	1069	1232	868	911	652	534	310	320	290	221	153	116	120	105	112	118	123	149	163	120	116	95	93	117	103	77	111	161	133	135	113	67	81	72	80	48	50	56	54	36	68	48	43	33	24	25	10	19	21	91	171	177	109	119	105	108	94	98	52	58	40	39	21	27	47	53	34	74	56	58	95	46	15	120	146	197	207	219	224	227	215	147	118	90	46	104",
                        "0	0	0	5	26	601	1214	852	1453	1478	2890	3483	3792	9648	6332	2545	1307	389	318	407	273	358	266	213	138	150	112	89	61	51	43	45	40	36	50	68	50	67	60	37	37	46	62	58	51	128	85	64	86	41	53	47	41	23	28	21	32	14	21	18	22	14	7	9	10	10	12	47	90	100	71	102	75	72	52	70	34	35	18	28	16	14	26	39	29	74	56	58	95	46	15	120	146	197	207	219	224	227	215	147	118	90	46	104",
                        730, # 10000,
                        ### from the rebounding paper, all ages assume the same ###
                        rep(0.9647,9),
                        c(   0.0122239664155569*0.46,
                             0.0209363491765736*1.01,
                             0.0304561389882463*0.79,
                             0.0482261102338912*0.44,
                             0.0780703865458093*0.49,
                             0.147326127082464*0.38,
                             0.285027801413735*0.63,
                             0.314214187084188*0.77)*0.12,
                        0.95)
  names(const.prms.start) <- c("symp-frac-davies", "steps-per-day",
                               "time-symp-to-hosp",
                               "tv-infect-period",
                               "tv-infect-period-endday",
                               "tv-incub-period",
                               "tv-incub-period-endday",
                               "tv-death-prob-home-00",
                               "tv-death-prob-home-10",
                               "tv-death-prob-home-20",
                               "tv-death-prob-home-30",
                               "tv-death-prob-home-40",
                               "tv-death-prob-home-50",
                               "tv-death-prob-home-60",
                               "tv-death-prob-home-70",
                               "tv-death-prob-home-80",
                               "tv-death-prob-home-endday",
                               "tv-vaccinees-endday",
                               "tv-vaccinees-00",
                               "tv-vaccinees-10",
                               "tv-vaccinees-20",
                               "tv-vaccinees-30",
                               "tv-vaccinees-40",
                               "tv-vaccinees-50",
                               "tv-vaccinees-60",
                               "tv-vaccinees-70",
                               "tv-vaccinees-80", 
                               "len-immunity",
                               "treat-success-00", 
                               "treat-success-10",
                               "treat-success-20",
                               "treat-success-30",
                               "treat-success-40",
                               "treat-success-50",
                               "treat-success-60",
                               "treat-success-70",
                               "treat-success-80",
                               # "treat-hosp-frac-00",
                               "treat-hosp-frac-10",
                               "treat-hosp-frac-20",
                               "treat-hosp-frac-30",
                               "treat-hosp-frac-40",
                               "treat-hosp-frac-50",
                               "treat-hosp-frac-60",
                               "treat-hosp-frac-70",
                               "treat-hosp-frac-80",
                               "tv-vac-efficacy")
  # 'tv-vac-efficacy-endday'
  
  # starting values for beta (length should == ncol(Z))
  
  beta.strt <- c(1.720765,0.6847355,0.2626981,0.1098193,0.1397066,0.0728003,
                 0.03618704,0.04517943,0.03937023,0.04179407,0.03797227,
                 0.01944581,0.02654224,0.0269852,0.02519945,0.07283215,
                 0.09123634,0.2011001,0.1871146,0.2267943,0.2221589,
                 0.2365618,0.1128427,0.0001972004,0.1545295,0.05364528,
                 0.08979767,0.09273189,0.1384901,0.2966932,0.1304654,
                 0.1186107,0.1357858,0.1464862,0.206332,0.1727486,0.1849924,
                 0.1540724,0.1405076,0.1779844,0.1377176,0.1761172,0.09987833,
                 0.1316974,0.1194551,0.1831427,0.2408771,0.1944747,0.1234956,
                 0.08781628,0.134314,0.2102527,0.237092,0.09222687,0.1457107,
                 0.2513297,0.1943738,0.4588658,0.2448166,0.1902813,0.1757216,
                 0.1784988,0.0002163204,0.0796581,0.09339351,0.337311,0.5882754,
                 0.5214095,0.5684911,0.7051709,0.8290295,1.387568,0.5067564,
                 0.4928869,0.4929972,0.4376321,0.5318749,0.3197462,0.4990711,
                 0.4585915,0.500883,0.1967081,0.6853107,0.5484323,0.5124942,
                 0.4867156,0.1986309,0.2586745,1.597719,1.089512,0.6585058,
                 0.5311778,0.6340663,1.403016,1.123881,1.525742,1.811137,
                 1.380535,3.008206,0.9147705,2.948115,1.549177,1.670181,
                 1.471259,1.689512,1.272016,1.668222,1.514672,1.576633,1.135684,
                 1.665079,1.281641,1.329348,1.254889,0.9466516,1.541413,
                 0.6542388,1.59326,0.4666356,2.051524,1.069608,0.9544311,
                 1.240686,1.648283,0.8999421,1.329196,1.454918,0.6616658,
                 1.611116,1.097197,1.280505,0.8618507,2.111266,1.144453,1.209791,
                 1.300466,1.166139,0.8433937,1.794371,1.27201,1.248367,1.281557,
                 1.277921,1.267835)
  
  beta = generate_beta2(start_beta = Z %*% beta.strt,
                        base_period_start = beta_base_period_start,
                        base_period = beta_base_period,
                        extend_period = sim_time,
                        transition_period = beta_transition_period,
                        multiplier = beta_multiplier,
                        amp = seas_amp, 
                        seas_dur = seas_dur, 
                        seas_start = seas_start)
  
  ### if baseline beta is adopted from other sources ###
  if(beta_adopt) {
    
    beta = generate_beta2(start_beta = adopted_beta,
                          base_period_start = beta_base_period_start,
                          base_period = beta_base_period,
                          extend_period = sim_time,
                          transition_period = beta_transition_period,
                          multiplier = beta_multiplier,
                          amp = seas_amp, 
                          seas_dur = seas_dur, 
                          seas_start = seas_start)
  }
               
  baseline_params <- list(beta = beta,
                          params = prms.start,
                          const.params = const.prms.start,
                          non.odesim.params=NULL,
                          introday = introday,
                          tf= length(beta) + 60,
                          data_endday = end.day,
                          odepath=odepath,
                          loc="RI", 
                          symp = NULL,
                          mean.report.rate = mean(rr.daily))
  return(baseline_params)
}
