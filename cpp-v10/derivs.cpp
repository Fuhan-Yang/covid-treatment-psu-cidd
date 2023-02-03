#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "derivs.h"
#include "essentials.h"
#include "prms.h"

extern prms* ppc;
extern bool G_B_USESOCIALCONTACTMATRIX;
extern double G_C_COMM[NUMAC][NUMAC];
extern double G_C_HOSP[NUMAC][NUMAC];
extern double G_C_ICU[NUMAC][NUMAC];
extern double G_C_VENT[NUMAC][NUMAC];


void derivs( double t, double *y, double *dydt)
{

    // clinical parameter changes during omicron 

    if(t >= 721) 
    {
    // death prob in hosp only apply for 70,80
     ppc->v_prob_HA4_D[7] = 0.05 / 2.0 * 0.21; // Lewnard et al nature medicine death aHR of omicron vs delta : 0.21
     ppc->v_prob_HA4_D[8] = 0.05 * 0.21 ; 

    // // mean time stay in hospital 
    ppc->v[ i_len_medicalfloor_hospital_stay ] = 10.7 * 1.24; // Lewnard et al
    
    // // prob to ICU 
    ppc->v_prob_HA_CA[0] = 0.304*0.5;    // Lewnard et al nature medicine: delta vs omicron aHR : 1:0.5
    ppc->v_prob_HA_CA[1] = 0.293*0.5;     // NOTE the estimated Lewnard et al (medRxiv, Apr 16) probabilities are used here
    ppc->v_prob_HA_CA[2] = 0.2825*0.5;    //      averaged over male/female equally
    ppc->v_prob_HA_CA[3] = 0.301*0.5;
    ppc->v_prob_HA_CA[4] = 0.463*0.5;
    ppc->v_prob_HA_CA[5] = 0.4245*0.5;
    ppc->v_prob_HA_CA[6] = 0.460*0.5;
    ppc->v_prob_HA_CA[7] = 0.4835*0.5;
    ppc->v_prob_HA_CA[8] = 0.416*0.5;

    // // death prob in ICU = prob not be ventilated * death prob on ventilator
    
    // // prob from ICU to ventilator 
    ppc->v_prob_CA_V[0] = 0.66*0.36; // Lewnard et al nature medicine: delta vs omicron aHR : 1: 0.36
    ppc->v_prob_CA_V[1] = 0.66*0.36;
    ppc->v_prob_CA_V[2] = 0.66*0.36;
    ppc->v_prob_CA_V[3] = 0.66*0.36;
    ppc->v_prob_CA_V[4] = 0.66*0.36;
    ppc->v_prob_CA_V[5] = 0.66*0.36;
    ppc->v_prob_CA_V[6] = 0.66*0.36;
    ppc->v_prob_CA_V[7] = 0.66*0.36;
    ppc->v_prob_CA_V[8] = 0.66*0.36;

    }
    int i,ac,byac,onac;
    double popsize=ppc->v[i_N];
    for (i = STARTD; i < STARTDHOSP+NUMAC; i++)
    {
        popsize -= y[i];
    }
    

    //printf("\n\n\t ************* %1.3f \t %1.3f \t %1.3f \n\n", ppc->v[0],ppc->v[1],ppc->v[2]); fflush(stdout);

    
    // ### 0 ### compute the force of infection 
    //
    double foi=0.0;
    double foi_ofcomm[NUMAC]; // force of infection of the community
    double foi_ofhosp[NUMAC]; // force of infection coming from hospitalized patients
    for(onac=0;onac<NUMAC;onac++) foi_ofcomm[onac]=0.0;
    for(onac=0;onac<NUMAC;onac++) foi_ofhosp[onac]=0.0;
    
    
    if( !G_B_USESOCIALCONTACTMATRIX )
    {
    
        for(i=NUME-2; i<NUME; i++) // loop through the last two stages of exposed individuals 
        {
            for(ac=0; ac<NUMAC; ac++)
            {
                foi += ppc->v[i_beta] * ppc->v[i_phi_incub] * y[STARTE + i*NUMAC + ac];   
            }
        }
        for(i=0; i<NUMAC*NUMA; i++) // loop through all asymptomatic individuals
        {
            foi += ppc->v[i_beta] * ppc->v[i_phi_asymp] * y[STARTA + i];   // STARTI is the starting index of all of the I-classes 
        }
        for(i=0; i<NUMAC*NUMI; i++) // you want to loop across all NUMI stages and all NUMAC ages of infected individuals
        {
            foi += ppc->v[i_beta] * y[STARTI + i];      // STARTI is the starting index of all of the I-classes 
        }                                               // no need to multiply by a relative infectiousness parameter because this is 1.0 here (the reference case)
        for(i=0; i<NUMAC*NUMHA; i++) // loop through all hospitalized individuals
        {
            foi += ppc->v[i_beta_hosp] * ppc->v[i_phi_hosp] * y[STARTHA + i];   
        }
        for(i=0; i<NUMAC; i++) // loop through all ICU individuals
        {
            foi += ppc->v[i_beta_icu] * ppc->v[i_phi_icu] * y[STARTCA + i];   
        }
        for(i=0; i<NUMV*NUMAC; i++) // loop through all ventilated individuals
        {
            foi += ppc->v[i_beta_vent] * ppc->v[i_phi_vent] * y[STARTV + i];   
        }
        for(i=0; i<NUMAC; i++) // loop through all CR individuals
        {
            foi += ppc->v[i_beta_icu] * ppc->v[i_phi_icu] * y[STARTCR + i];   
        }
        for(i=0; i<NUMAC; i++) // loop through all HR individuals
        {
            // foi += ppc->v[i_beta_hosp] * ppc->v[i_phi_hosp] * y[STARTHR + i];    // NOTE this is set to zero-contribution to the FOI right now because        
        }                                                                           // these individuals are on day 10-15 of their infection, and Wolfel et al (Nature, 2020) suggest that
                                                                                    // positivity by virus-culture should be low by this time 
    }

    double beta00;
    double beta10;
    double beta20;
    double beta30;
    double beta40; // this is the true beta (transmission and contact rate) experienced by 40-49 year-olds
    double beta50; // this is the true beta (transmission and contact rate) experienced by 50-59 year-olds
    double beta60; // this is the true beta (transmission and contact rate) experienced by 60-69 year-olds
    double beta70; // this is the true beta (transmission and contact rate) experienced by 70-79 year-olds
    double beta80; // this is the true beta (transmission and contact rate) experienced by 80-89 year-olds

    
    
    // 
    // in the lines below, if the current beta is too low, then reset the beta00 to beta80 values to a minimum
    // acceptable level as these (older) age groups will not be able to have very little contact
    //
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_00 ] )  beta00 = ppc->v[ i_min_relbeta_00 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta00 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_10 ] )  beta10 = ppc->v[ i_min_relbeta_10 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta10 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_20 ] )  beta20 = ppc->v[ i_min_relbeta_20 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta20 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_30 ] )  beta30 = ppc->v[ i_min_relbeta_30 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta30 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_40 ] )  beta40 = ppc->v[ i_min_relbeta_40 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta40 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_50 ] )  beta50 = ppc->v[ i_min_relbeta_50 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta50 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_60 ] )  beta60 = ppc->v[ i_min_relbeta_60 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta60 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_70 ] )  beta70 = ppc->v[ i_min_relbeta_70 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta70 = ppc->v[i_beta];
    if( (ppc->v[i_beta] / ppc->v[ i_prelockdown_beta ]) < ppc->v[ i_min_relbeta_80 ] )  beta80 = ppc->v[ i_min_relbeta_80 ] * ppc->v[ i_prelockdown_beta ];
    else                                                                                beta80 = ppc->v[i_beta];



    

    
    if( G_B_USESOCIALCONTACTMATRIX )
    {
        for(onac=0; onac<NUMAC; onac++) // force of infection ON a particular age class
        {
            
            // ### 1 ### loop through all the last two classed of exposed individuals, and add their FOI
            for(i=NUME-2; i<NUME; i++) // loop through the last two stages of exposed individuals 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class
                {
                    //if(onac<4)  foi_ofcomm[onac] += ppc->v[i_beta] * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==0) foi_ofcomm[onac] +=         beta00 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==1) foi_ofcomm[onac] +=         beta10 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==2) foi_ofcomm[onac] +=         beta20 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==3) foi_ofcomm[onac] +=         beta30 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==4) foi_ofcomm[onac] +=         beta40 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==5) foi_ofcomm[onac] +=         beta50 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==6) foi_ofcomm[onac] +=         beta60 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];
                    if(onac==7) foi_ofcomm[onac] +=         beta70 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];   
                    if(onac==8) foi_ofcomm[onac] +=         beta80 * ppc->v[i_phi_incub] * G_C_COMM[onac][byac] * y[STARTE + i*NUMAC + byac];   
                    
                }
            }
            
            // ### 2 ### loop through all asymptomatic individuals, and add their FOI
            for(i=0; i<NUMA; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class
                {
                    //if(onac<4)  foi_ofcomm[onac] += ppc->v[i_beta] * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];   
                    if(onac==0) foi_ofcomm[onac] +=         beta00 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];
                    if(onac==1) foi_ofcomm[onac] +=         beta10 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];
                    if(onac==2) foi_ofcomm[onac] +=         beta20 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];
                    if(onac==3) foi_ofcomm[onac] +=         beta30 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];
                    if(onac==4) foi_ofcomm[onac] +=         beta40 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];
                    if(onac==5) foi_ofcomm[onac] +=         beta50 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];   
                    if(onac==6) foi_ofcomm[onac] +=         beta60 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];   
                    if(onac==7) foi_ofcomm[onac] +=         beta70 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];   
                    if(onac==8) foi_ofcomm[onac] +=         beta80 * ppc->v[i_phi_asymp] * G_C_COMM[onac][byac] * y[STARTA + i*NUMAC + byac];   
                }
            }
            
            // ### 3 ### loop through all infected+symptomatic individuals, and add their FOI
            for(i=0; i<NUMI; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class - there is no phi-parameter because it is 1.0 here (this is the reference class for infectivity)
                {
                    //if(onac<4)  foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] * ppc->v[i_beta] * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];
                    if(onac==0) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta00 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];
                    if(onac==1) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta10 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];
                    if(onac==2) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta20 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];
                    if(onac==3) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta30 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];
                    if(onac==4) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta40 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];
                    if(onac==5) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta50 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];   
                    if(onac==6) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta60 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];   
                    if(onac==7) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta70 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];   
                    if(onac==8) foi_ofcomm[onac] += ppc->v[i_selfisolation_factor] *         beta80 * G_C_COMM[onac][byac] * y[STARTI + i*NUMAC + byac];   
                }
            }

            // ### 4 ### loop through all hospitalized individuals, and add their FOI
            for(i=0; i<NUMHA; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class 
                {
                    foi_ofhosp[onac] += ppc->v[i_beta_hosp] * ppc->v[i_phi_hosp] * G_C_HOSP[onac][byac] * y[STARTHA + i*NUMAC + byac];   
                }
            }
            for(i=0; i<1; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class 
                {
                    foi_ofhosp[onac] += ppc->v[i_beta_hosp] * ppc->v[i_phi_hosp_recovering] * G_C_HOSP[onac][byac] * y[STARTHR + byac];   
                }
            }

            // ### 5 ### loop through all ICU individuals, and add their FOI
            for(i=0; i<1; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class 
                {
                    foi_ofhosp[onac] += ppc->v[i_beta_icu] * ppc->v[i_phi_icu] * G_C_ICU[onac][byac] * y[STARTCA + byac];   
                }
            }
            for(i=0; i<1; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class 
                {
                    foi_ofhosp[onac] += ppc->v[i_beta_icu] * ppc->v[i_phi_icu] * G_C_ICU[onac][byac] * y[STARTCR + byac];   
                }
            }
            
            // ### 6 ### loop through all ICU individuals on ventilators, and add their FOI
            for(i=0; i<NUMV; i++) 
            {
                for(byac=0; byac<NUMAC; byac++) // infection BY a particular age class
                {
                    foi_ofhosp[onac] += ppc->v[i_beta_vent] * ppc->v[i_phi_vent] * G_C_VENT[onac][byac] * y[STARTV + i*NUMAC + byac];   
                }
            }
            
        }
        
    }    

    
    //NOTE the variables below are transitions rate (tr) variables between different compartment types; e.g. trv is the transition rate among V-classes
    //
    // this is the transition rate among the E-classes 
    double tre = ((double)NUME) / ppc->v[i_len_incub_period];
    
    // this is the transition rate among the I-classes
    double tri1 = (((double)NUMI)/2.0) / ppc->v[ i_len_symptomatic_infectious_period_phase_1 ];
    double tri2 = (((double)NUMI)/2.0) / ppc->v[ i_len_symptomatic_infectious_period_phase_2 ];
    
    // this is the transition rate among the HA-classes
    double trha = ((double)NUMHA) / ppc->v[ i_len_medicalfloor_hospital_stay ];
    
    
    // this is the transition rate among the A-classes  this is hardcoded as 5.333 days for now
    double tra = 0.75;
    
    // this is the transition rate out of the 'acute ICU' stage meaning your mean time here is two days
    double trca = 0.5;
    
    // this is the transition rate among the V-classes ... each stage is about 1.8 days
    // meaning you have about 10.8 days on a ventilator in the 6 V classes
    // double trv = 1.0 / 1.8;
    double trv = 1.0 / ( (ppc->v[i_mean_time_vent]/6.0) );
    
    double trhr = 0.4;          // NOTE BethG says this should be about 2-3 days.
    double trcr = 1.0/2.6;      // NOTE Bhatraju paper has this at 2.6 days; no other source
    
    
    double from_I2 = 0.0;   // amount from I2 (helper in case of treatment)
    
        
    // ### 1 ### from index=0 to index=NUMAC-1 you have the S-classes, one for each age class
    //
    // first the S-classes
    if( !G_B_USESOCIALCONTACTMATRIX )
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[ac] = - foi * y[ac] / popsize; 
        }
    }
    
    if( G_B_USESOCIALCONTACTMATRIX )
    {
        for(onac=0;onac<NUMAC;onac++)
        {
            dydt[onac] = - ppc->v_rel_susc[onac] * ( foi_ofcomm[onac] + foi_ofhosp[onac] ) * y[onac] / popsize; 
        }
    }

    // immunity waning (R, RHOSP, VAX back to S)
    if (t >= 524 && ppc->v[i_len_immunity] > 0){
        for (ac = 0; ac < NUMAC; ac++){
            dydt[ac] += 1.0/ppc->v[i_len_immunity] * (y[STARTR + ac] + y[STARTRHOSP + ac] + y[STARTVAX + ac]);
        }
    }
    
    // death rate for S class
    /* for (ac = 0; ac < NUMAC; ac++)
    {
        dydt[ac] -= y[ac] * ppc->v_deathrate[ac];
    } */
    
    
    
    // ### 2 ### from index=NUMAC to index=NUMAC+NUME*NUMAC you have the E-classes, one for each age class
    // there should be 54 E-classes (54 = NUME*NUMAC)
    if( !G_B_USESOCIALCONTACTMATRIX )
    {
        for(i=0;i<NUME;i++)
        {
            for(ac=0;ac<NUMAC;ac++)
            {
                if(i==0) // meaning this is the E_1 class and the S-classes have to flow into it
                {
                    dydt[STARTE + i*NUMAC + ac] = foi * y[ac] / popsize - tre * y[STARTE + i*NUMAC + ac];
                }
                else // these are the classes E_2 and higher
                {
                    dydt[STARTE + i*NUMAC + ac] = tre * y[STARTE + (i-1)*NUMAC + ac] - tre * y[STARTE + i*NUMAC + ac];
                }
                
            }
        }
    }

    if( G_B_USESOCIALCONTACTMATRIX )
    {
        for(i=0;i<NUME;i++)
        {
            for(ac=0;ac<NUMAC;ac++)
            {
                if(i==0) // meaning this is the E_1 class and the S-classes have to flow into it
                {
                    dydt[STARTE + i*NUMAC + ac] = ppc->v_rel_susc[ac] * ( foi_ofcomm[ac] + foi_ofhosp[ac] ) * y[ac] / popsize - tre * y[STARTE + i*NUMAC + ac];
                }
                else // these are the classes E_2 and higher
                {
                    dydt[STARTE + i*NUMAC + ac] = tre * y[STARTE + (i-1)*NUMAC + ac] - tre * y[STARTE + i*NUMAC + ac];
                }
                
            }
        }
    }
    
    
    // ### 3 ###    ASYMPTOMATIC, INFECTED, SOMEWHAT INFECTIOUS CLASSES (the A-class)
    //              there should be NUMAC*NUMA = 36 of these 
    for(i=0;i<NUMA;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the A_1 class and a fraction of the E_final class has to flow into it
            {
                dydt[STARTA + i*NUMAC + ac] = ppc->v_prob_E_A[ac]*tre*y[STARTE + (NUME-1)*NUMAC + ac] - tra * y[STARTA + i*NUMAC + ac]; 
            }
            else // meaning this is the A_2 class or higher
            {
                dydt[STARTA + i*NUMAC + ac] = tra*y[STARTA + (i-1)*NUMAC + ac] - tra * y[STARTA + i*NUMAC + ac];
            }
            
        }
    }
    
    
    
    // ### 4 ###    INFECTED, INFECTIOUS, AND SYMPTOMATIC CLASSES (THE I-CLASS)
    // there are 36 I-classes (36 = NUMI*NUMAC)
    for(i=0;i<NUMI;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the I_1 class and a fraction of the E_final class has to flow into it
            {
                dydt[STARTI + i*NUMAC + ac] = (1.0-ppc->v_prob_E_A[ac])*tre*y[STARTE + (NUME-1)*NUMAC + ac] - tri1 * y[STARTI + i*NUMAC + ac]; 
            }
            else if(i==1) // meaning this is the I_2 class 
            {
                dydt[STARTI + i*NUMAC + ac  ] = tri1 * y[STARTI + (i-1)*NUMAC + ac] - tri1 * y[STARTI + i*NUMAC + ac];
            }
            else if(i==2) // meaning this is the I_3 class .... NOTE you have to be **non-hospitalized** to flow into this class
            {
                if (t >= ppc->v[i_treatment_beginday]){
                    // dydt[STARTI + i*NUMAC + ac] = (1.0-ppc->v_prob_I2_H_treatment[ac])*tri1 * y[STARTI + (i-1)*NUMAC + ac] - tri2 * y[STARTI + i*NUMAC + ac];

                    // outflow after tri2 waiting time
                    dydt[STARTI + i*NUMAC + ac] = - tri2 * y[STARTI + i*NUMAC + ac];

                    // from I2
                    from_I2 = tri1 * y[STARTI + (i-1)*NUMAC + ac] ;

                    // inflow from non-treatment group
                    // dydt[STARTI + i*NUMAC + ac] += (1.0 - ppc->v_treatment_coverage[ac]) * (1.0 - ppc->v_prob_I2_H[ac]) * from_I2;
                    dydt[STARTI + i*NUMAC + ac] += (1.0 - ppc->tvp_treatment_coverage->get_current_data(ac) ) * (1.0 - ppc->v_prob_I2_H[ac]) * from_I2;
                    
                    // inflow from treatment group (treatment failure and not hospitalized)
                    // dydt[STARTI + i*NUMAC + ac] += ppc->v_treatment_coverage[ac] * (1.0 - ppc->v_treatment_success[ac]) * (1.0 - ppc->v_prob_I2_H_treatment[ac]) * from_I2;
                    dydt[STARTI + i*NUMAC + ac] += ppc->tvp_treatment_coverage->get_current_data(ac) * (1.0 - ppc->v_treatment_success[ac]) * (1.0 - ppc->v_prob_I2_H_treatment[ac]) * from_I2;

                } else {
                    dydt[STARTI + i*NUMAC + ac] = (1.0-ppc->v_prob_I2_H[ac])*tri1 * y[STARTI + (i-1)*NUMAC + ac] - tri2 * y[STARTI + i*NUMAC + ac];
                }
            }
            else if(i==3) // meaning this is the I_4 class 
            {
                dydt[STARTI + i*NUMAC + ac] = tri2 * y[STARTI + (i-1)*NUMAC + ac] - tri2 * y[STARTI + i*NUMAC + ac];
            }
            else
            {
                assert(false);   
            }
                        
        }
    }

    
    
    // ### 5 ###    HOSPITALIZED INDIVIDUALS, ACUTE-PHASE, MEANING THAT THEY RECENTLY ARRIVED AT THE HOSPITAL AND THEY ARE AT RISK FOR WORSENING CONDITIONS
    //              there should be NUMAC*NUMHA = 36 of these as well 
    for(i=0;i<NUMHA;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the HA_1 class and a fraction of the I_2 class flows into here
            {
                if (t >= ppc->v[i_treatment_beginday]){
                    // dydt[STARTHA + i*NUMAC + ac] = ppc->v_prob_I2_H_treatment[ac]*tri1*y[STARTI + NUMAC + ac] - trha * y[STARTHA + i*NUMAC + ac]; 

                    // outflow after trha waiting time
                    dydt[STARTHA + i*NUMAC + ac] = - trha * y[STARTHA + i*NUMAC + ac];

                    // from I2
                    from_I2 = tri1 * y[STARTI + NUMAC + ac] ;

                    // inflow from non-treatment group
                    // dydt[STARTHA + i*NUMAC + ac] += (1.0 - ppc->v_treatment_coverage[ac]) * ppc->v_prob_I2_H[ac] * from_I2;
                    dydt[STARTHA + i*NUMAC + ac] += (1.0 - ppc->tvp_treatment_coverage->get_current_data(ac)) * ppc->v_prob_I2_H[ac] * from_I2;

                    // inflow from treatment group (treatment failures and hospitalized)
                    // dydt[STARTHA + i*NUMAC + ac] += ppc->v_treatment_coverage[ac] * (1.0 - ppc->v_treatment_success[ac]) * ppc->v_prob_I2_H_treatment[ac] * from_I2;
                    dydt[STARTHA + i*NUMAC + ac] += ppc->tvp_treatment_coverage->get_current_data(ac) * (1.0 - ppc->v_treatment_success[ac]) * ppc->v_prob_I2_H_treatment[ac] * from_I2;

                } else {
                    dydt[STARTHA + i*NUMAC + ac] = ppc->v_prob_I2_H[ac]*tri1*y[STARTI + NUMAC + ac] - trha * y[STARTHA + i*NUMAC + ac]; 
                }
            }
            else if(i==1) // meaning this is the HA_2 class; you only go into HA2 if you did not progress to the ICU
            {
                dydt[STARTHA + i*NUMAC + ac] = (1.0 - ppc->v_prob_HA_CA[ac]) * trha * y[STARTHA + (i-1)*NUMAC + ac] - trha * y[STARTHA + i*NUMAC + ac];                
            }
            else // meaning this is the HA_3 class or higher
            {
                dydt[STARTHA + i*NUMAC + ac] = trha * y[STARTHA + (i-1)*NUMAC + ac] - trha * y[STARTHA + i*NUMAC + ac];
            }
            
        }
    }
    
    
    
    // ### 6 ###    CRITICAL-CARE/ICU INDIVIDUALS, ACUTE-PHASE, MEANING THAT THEY RECENTLY ARRIVED IN THE ICU AND THEY ARE AT RISK FOR WORSENING CONDITIONS
    //              THERE IS ONLY ONE "STAGE" OF ACUTE ICU FOR NOW
    //
    for(i=0;i<1;i++) // number of ICU stages is 1, so no loop here (just a placeholder in case this is changed)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[STARTCA + ac] = 0.0;
            
            // the only stage ("stg") of the HA-classes that allows progression to the ICU is the first one, stg=0
            int stg=0;
            dydt[STARTCA + ac] = ppc->v_prob_HA_CA[ac] * trha * y[STARTHA + stg*NUMAC + ac] - trca * y[STARTCA + ac];
        }
    }
    
    

    
    // ### 7 ###    INDIVIDUALS ON VENTILATORS
    //              there should be NUMAC*NUMV = 54 of these classes
    
    // this is the index of the 'special' V-class (currently, the fourth class) from which a patient can die
    int VK = 3;
    
    for(i=0;i<NUMV;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the V_1 class and a fraction of the CA class flows into here
            {
                dydt[STARTV + i*NUMAC + ac] = ppc->v_prob_CA_V[ac]*trca*y[STARTCA + ac] - trv * y[STARTV + i*NUMAC + ac]; 
            }
            else if(i==VK+1) // meaning, this is the fifth class and patients who progress to death will not enter this class
            {
                dydt[STARTV + i*NUMAC + ac] = trv * (1.0-ppc->v_prob_V_D[ac]) * y[STARTV + (i-1)*NUMAC + ac] - trv * y[STARTV + i*NUMAC + ac];
            }
            else // meaning this is the V_2 class or higher
            {
                dydt[STARTV + i*NUMAC + ac] = trv * y[STARTV + (i-1)*NUMAC + ac] - trv * y[STARTV + i*NUMAC + ac];
            }
            
        }
    }

    
    // ### 8 ###    CRITICAL-CARE/ICU INDIVIDUALS, RECOVERING-PHASE, MEANING THAT THEY JUST GOT OFF A VENTILATOR
    //              THERE IS ONLY ONE "STAGE" OF RECOVERING ICU FOR NOW
    for(i=0;i<1;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[STARTCR + ac] = trv * y[STARTV + (NUMV-1)*NUMAC + ac] - trcr * y[STARTCR + ac];            
        }
    }

    
    // ### 9 ###    HOSPITALIZED INDIVIDUALS, RECOVERING-PHASE, MEANING THAT THEY JUST GOT HERE FROM THE ICU
    //              THERE IS ONLY ONE "STAGE" OF "HR" FOR NOW
    for(i=0;i<1;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[STARTHR + ac] = (1.0-ppc->v_prob_CR_D[ac]) * trcr * y[STARTCR + ac]  + (1.0-ppc->v_prob_CA_V[ac]-ppc->v_prob_CA_D[ac]) * trca * y[STARTCA + ac] - trhr * y[STARTHR + ac]; 
        }
    }
    
    
    
    // ### 10 ###    DEAD INDIVIDUALS (D) WHO DIED AT HOME 
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        dydt[STARTD + ac] = ppc->v_prob_I4_D[ac] * tri2 * y[STARTI + (NUMI-1)*NUMAC + ac]; 
    }

    
    
    // ### 11 ###    DEAD INDIVIDUALS (DHOSP) WHO DIED IN THE HOSPITAL - 
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        // set deriv equal to zero
        dydt[STARTDHOSP + ac] = 0.0;
        
        // add asymptomatics recovering from A4 - WE ASSUME THAT THE DEATH RATE FORM THE ASYMPTOMATIC CLASS IS ZERO
        // dydt[STARTD + ac] += tra * y[STARTA + (NUMA-1)*NUMAC + ac]; 
        
        // add infecteds coming from I4 
        // dydt[STARTD + ac] += ppc->v_prob_I4_D[ac] * tri2 * y[STARTI + (NUMI-1)*NUMAC + ac]; 

        // add HA-individuals coming from HA4 
        dydt[STARTDHOSP + ac] += ppc->v_prob_HA4_D[ac] * trha * y[STARTHA + (NUMHA-1)*NUMAC + ac]; 
        
        // add HR-individuals coming from HR; same death rate as HA4
        dydt[STARTDHOSP + ac] += ppc->v_prob_HA4_D[ac] * trhr * y[STARTHR + ac]; 

        // add CA-individuals coming from CA 
        dydt[STARTDHOSP + ac] += ppc->v_prob_CA_D[ac] * trca * y[STARTCA + ac]; 
        
        // add in V-individuals who will die (just V4 for now)
        dydt[STARTDHOSP + ac] += trv * ppc->v_prob_V_D[ac] * y[STARTV + VK*NUMAC + ac];

        // add CR-individuals coming from CR 
        dydt[STARTDHOSP + ac] += ppc->v_prob_CR_D[ac] * trcr * y[STARTCR + ac]; 
    }

    
    
    
    // ### 12 ###    RECOVERED INDIVIDUALS (R) WHO RECOVERED AT HOME
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        // set deriv equal to zero
        dydt[STARTR + ac] = 0.0;
        
        // add asymptomatics recovering from A4
        dydt[STARTR + ac] += tra * y[STARTA + (NUMA-1)*NUMAC + ac]; 
        
        // add infecteds coming from I4
        dydt[STARTR + ac] += (1.0-ppc->v_prob_I4_D[ac]) * tri2 * y[STARTI + (NUMI-1)*NUMAC + ac]; 

        // add HA-individuals coming from HA4
        //dydt[STARTR + ac] += (1.0 - ppc->v_prob_HA4_D[ac] - ppc->v_prob_HA_CA[ac]) * trha * y[STARTHA + (NUMHA-1)*NUMAC + ac];
        // add HR-individuals coming from HR
        //dydt[STARTR + ac] += (1.0-ppc->v_prob_HA4_D[ac]) * trhr * y[STARTHR + ac];

        // add infecteds from I2 after successful treatment
        if (t >= ppc->v[i_treatment_beginday]){
            // dydt[STARTR + ac] += ppc->v_treatment_coverage[ac] * ppc->v_treatment_success[ac] * tri1 * y[STARTI + NUMAC + ac] ;
            dydt[STARTR + ac] += ppc->tvp_treatment_coverage->get_current_data(ac) * ppc->v_treatment_success[ac] * tri1 * y[STARTI + NUMAC + ac] ; 
        }

        // immunity waning (R back to S)
        if (t >= 524 && ppc->v[i_len_immunity] > 0){
            dydt[STARTR + ac] -= 1.0/ppc->v[i_len_immunity] * y[STARTR + ac];
        }
        
    }

    // ### 13 ###    RECOVERED INDIVIDUALS (RHOSP) WHO RECOVERED AT THE HOSPITAL AND WERE DISCHARGED
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        // set deriv equal to zero
        dydt[STARTRHOSP + ac] = 0.0;
        
        // add asymptomatics recovering from A4
        //dydt[STARTRHOSP + ac] += tra * y[STARTA + (NUMA-1)*NUMAC + ac]; 
        
        // add infecteds coming from I4
        // dydt[STARTRHOSP + ac] += (1.0-ppc->v_prob_I4_D[ac]) * tri2 * y[STARTI + (NUMI-1)*NUMAC + ac]; 

        // add HA-individuals coming from HA4
        dydt[STARTRHOSP + ac] += (1.0 - ppc->v_prob_HA4_D[ac]) * trha * y[STARTHA + (NUMHA-1)*NUMAC + ac];
        // add HR-individuals coming from HR
        dydt[STARTRHOSP + ac] += (1.0-ppc->v_prob_HA4_D[ac]) * trhr * y[STARTHR + ac];
        
        // immunity waning (RHOSP back to S)
        if (t >= 524 && ppc->v[i_len_immunity] > 0){
            dydt[STARTRHOSP + ac] -= 1.0/ppc->v[i_len_immunity] * y[STARTRHOSP + ac];
        }
    }

    // ### 14 ###    VACCINATED CLASS
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        // set deriv equal to zero
        dydt[STARTVAX + ac] = 0.0;
        
        // immunity waning (V back to S)
        if (t >= 524 && ppc->v[i_len_immunity] > 0){
            dydt[STARTVAX + ac] -= 1.0/ppc->v[i_len_immunity] * y[STARTVAX + ac];
        }
    }
    
    
    
    // ### 15 ###    CUMULATIVE SYMPTOMATIC INCIDENCE CLASSES
    //               the 9 J-classes
    for(ac=0;ac<NUMAC;ac++)
    {
        dydt[STARTJ+ac] = (1.0-ppc->v_prob_E_A[ac]) * tre * y[STARTE + (NUME-1)*NUMAC + ac];
    }
    
    
    // ### 16 ###    CUMULATIVE HOSPITALIZATION INCIDENCE CLASSES
    //               finally the 9 K-classes
    for(ac=0;ac<NUMAC;ac++)
    {
        if (t >= ppc->v[i_treatment_beginday]){
            // dydt[STARTK+ac] = ( ppc->v_prob_I2_H_treatment[ac] ) * tri1 * y[STARTI + NUMAC + ac] ;   

            // from I2
            from_I2 = tri1 * y[STARTI + NUMAC + ac] ;

            // inflow from non-treatment group
            // dydt[STARTK+ac] = (1.0 - ppc->v_treatment_coverage[ac]) * ppc->v_prob_I2_H[ac] * from_I2;
            dydt[STARTK+ac] = (1.0 - ppc->tvp_treatment_coverage->get_current_data(ac)) * ppc->v_prob_I2_H[ac] * from_I2;

            // inflow from treatment group (treatment failures and hospitalized)
            // dydt[STARTK+ac] += ppc->v_treatment_coverage[ac] * (1.0 - ppc->v_treatment_success[ac]) * ppc->v_prob_I2_H_treatment[ac] * from_I2;
            dydt[STARTK+ac] += ppc->tvp_treatment_coverage->get_current_data(ac) * (1.0 - ppc->v_treatment_success[ac]) * ppc->v_prob_I2_H_treatment[ac] * from_I2;

        } else {
            dydt[STARTK+ac] = ( ppc->v_prob_I2_H[ac] ) * tri1 * y[STARTI + NUMAC + ac] ;
        }
    }
    
    
    
}

