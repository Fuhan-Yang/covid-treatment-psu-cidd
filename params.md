| param | description | current value | new value | ref| notes |
| :-: | :-: | :-: | :-: | :-: |  :-: | 
| hard-coded | (a)symptomatic probability  | 0.71, 0.79, 0.73, 0.67, 0.60, 0.51, 0.37, 0.31, 0.31 (Davies et al) | 0.467 (0-18); 0.321 (19-59); 0.197 (60+) | Sah 2021 | Meta analysis of papers published through Apr 2021, proportions of asymptomatic cases |
| hard-coded | incubation period | 6 | omicron: 3.42 | Yu et al 2022 JAMA Open | Meta-analysis on the incubation periods of multiple variants |
| hard-coded | infectious period |5 |6 | Boucau et al 2022, NEJM (https://pubmed.ncbi.nlm.nih.gov/35767428/) ||
| `death-prob-home-60` | probability that a patient aged 60-69 dies of COVID-19 without being hospitalized | RI: 0.013 | 0.013*0.22|Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `death-prob-home-70` | probability that a patient aged 70-79 dies of COVID-19 without being hospitalized | RI:0.04 | 0.04*0.26 | Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `death-prob-home-80` | probability that a patient aged 80 and above dies of COVID-19 without being hospitalized | RI: 0.227 | 0.227*0.46  | Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `time-symp-to-hosp` | average length from symptomatic to hospitalized | RI: 4| 4 |Reichert et al 2022, Lancet Microbe supp. (https://pubmed.ncbi.nlm.nih.gov/36057266/) | | 
| `tv-hosp-frac-10` | hospitalization probablity for symptomatic cases aged 0-9 and 10-19  | RI: 0.012 | 0.012*0.83  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-20` | hospitalization probablity for symptomatic cases aged 20-29  | RI: 0.021|0.021*0.55  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-30` | hospitalization probablity for symptomatic cases aged 30-39  | RI: 0.03|0.03*0.44  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-40` | hospitalization probablity for symptomatic cases aged 40-49  | RI: 0.048|0.048*0.33  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-50` | hospitalization probablity for symptomatic cases aged 50-59  | RI: 0.078|0.078*0.26  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-60` | hospitalization probablity for symptomatic cases aged 60-69  | RI: 0.147|0.147*0.25  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-70` | hospitalization probablity for symptomatic cases aged 70-79  | RI: 0.285|0.285*0.36  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-80` | hospitalization probablity for symptomatic cases aged 80 and above  | RI: 0.314|0.314*0.47  |Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
| `tv-hosp-frac-x0` | hospitalization probablity for symptomatic cases | RI | delta/BA.1: 3.41 | Wolter et al Nat Comm | cohort study in South Africa comparing hosp odd ratio among Delta BA.1 BA.4 BA.5 for people infected with COVID-19 |
| hard-coded | average fraction of patients progressing from hospital standard floor care to ICU | RI: | delta: 0.104; omicron:0.05 |Bouzid et al 2022 Annuals of Internal Medicine | retrospective cohort study on patients in ED in Paris |
|  |  | RI: | delta: 20.9%; omicron: 7% | Goethem et al 2022 viruses | retrospective cohort study on hospitalized patients in Belgium; median age: 78(67-86)  |
|  |  | RI: | delta vs omicron aHR: 0.5 | Lewnard et al nature medicine 2022 | 222,688 omicron+ patients vs23,305 delta+ patients confirmed by PCR, recruited during Dec 2021 to Mar 2022 in Southern California | 
|  |  | RI: | delta vs omicron: 2 | Strasser et al JAMA NO | retrospective cohort study on covid-positive patients in New England, vaccinated and unvaccinated; age: 40-50 |
|  |  | RI: | delta vs SGTF-omicron: 1:0.3 | Wolten et al Lancet | cohort study on pcr-positive covid patients in South Africa; the ratio is aOR of severe diseases among hospitalized patients |
|  |  | RI: | pre-omicron: 4.3% vs omicron: 1% | Abdullah et al Int Journal Infec Diseases | hospitalized patients since Nov 14 2021(might included delta) are compared to the hospitalized patients before omicron
| hard-coded | average length of hospital standard floor care stay | 10.7 (Lewnard et al) |Delta: 5; Omicron:1  | Bouzid et al 2022 Annuals of Internal Medicine | 
| | | 10.7 (Lewnard et al) |Delta: 10; Omicron: 8  | Goethem et al 2022 viruses | retrospective cohort study on hospitalized patients in Belgium; median age: 78(67-86)  |
|  |  | 10.7 (Lewnard et al) | pre-omicron: 8.8 vs omicron: 4 | Abdullah et al Int Journal Infec Diseases | hospitalized patients since Nov 14 2021(might included delta) are compared to the hospitalized patients before omicron
| hard-coded | probability of death among non-ICU hospitalized patients |  | delta: 17.6%; omicron: 11.6% | Goethem et al 2022 viruses | retrospective cohort study on hospitalized patients in Belgium; median age: 78(67-86) |
|  |  |  | pre-omicron: 21.3% vs omicron: 4.5% | Abdullah et al Int Journal Infec Diseases | hospitalized patients since Nov 14 2021(might included delta) are compared to the hospitalized patients before omicron
|  |  |  | pre-omicron/omicron aHR: 0.68 | Krutikov et al Lancet Health Longev | Residents who were positive in long-term care facilities in England from Sep 1 2021 to Feb 1 2022 | 
| hard-coded | average length of ICU stay | | Delta: 9; Omicron:5  | Goethem et al 2022 viruses | retrospective cohort study on hospitalized patients in Belgium; median age: 78(67-86)  |
| `prob-icu-vent` | probability of progressing from ICU to mechanical ventilation | RI: 0.66 | delta:0.37; omicron: 0.425 | Bouzid et al 2022 Annuals of Internal Medicine | retrospective cohort study on patients in ED in Paris  |
|  |  | RI: 0.66 | delta:0.38 omicron: 0.26 | Goethem et al 2022 viruses | retrospective cohort study on hospitalized patients in Belgium; median age: 78(67-86) |
|  |  | RI: 0.66 | delta adjusted-OR vs BA2: 4.36; omicron adjusted OR vs BA2: 3.55 (delta/omicron: 1.23) | Strasser et al JAMA NO | retrospective cohort study on covid-positive patients in New England, vaccinated and unvaccinated; age: 40-50 |
|  |  | RI: | delta vs omicron aHR: 0.36 | Lewnard et al nature medicine 2022 | 222,688 omicron+ patients vs23,305 delta+ patients confirmed by PCR, recruited during Dec 2021 to Mar 2022 in Southern California |
| `mean-time-vent` | mean number of days a surviving patient spends on a ventilator | RI: 10.58 | assume the same |  |  |
|  |  | RI: | delta vs SGTF-omicron: 1:0.3 | Wolten et al Lancet | cohort study on pcr-positive covid patients in South Africa; the ratio is aOR of severe diseases among hospitalized patients |
| hard-coded | probability of dying given the patient is on ventilator | 0.03125, 0.05119, 0.15, 0.15, 0.400*`dev-ventdeath-mid`, 0.460*`dev-ventdeath-mid`, 0.585*`dev-ventdeath-mid`, 0.7, 0.9 (Lewnard et al; Graselli; Yang) | assume the same |  |  |
| hard-coded | death probablity among Covid-19 positive patients | | prob ratio: delta/omicron: No observed death below 30; 0.28(30-39); 0.25(40-49); 0.16(50-59); 0.22(60-69);0.26(70-79);0.46(80+)| Nyberg et al Lancet 2022  | cohort study in England comparing hosp risk between Delta and Omicron(B.1.1) |
