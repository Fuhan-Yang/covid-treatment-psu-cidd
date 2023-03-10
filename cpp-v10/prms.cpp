//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "essentials.h"
#include "assert.h"
#include "prms.h"
// #include "time_varying_prms_age.h"
// #include "time_varying_prms_nonage.h"


// constructor
prms::prms()
{
    ///// remember to delete these pointers in destructor !!
    tvp_contact_rate = new time_varying_prms_age<double>();
    tvp_contact_coeff = new time_varying_prms_age<double>();
    tvp_hosp_frac =  new time_varying_prms_age<double>();
    tvp_vac_num =  new time_varying_prms_age<double>();
    tvp_vac_efficacy = new time_varying_prms_nonage<double>();
    tvp_dev_icu_frac = new time_varying_prms_nonage<double>();
    tvp_dev_len_hospstay = new time_varying_prms_nonage<double>();

    tvp_incubation_period = new time_varying_prms_nonage<double>();
    tvp_infectious_period = new time_varying_prms_nonage<double>();

    tvp_death_home = new time_varying_prms_age<double>();

    // printf("%d\n", tvp_contact_rate->good_integrity() );
    // tvp_contact_rate->advance_indices();
    // printf("%d\n", tvp_contact_rate->get_day_index() );
    // printf("%d\n", tvp_contact_rate->get_data_index() );

    v.insert( v.begin(), num_params, 0.0 );
    v[i_treatment_beginday] = -999;
    v[i_len_immunity] = -999; // immunity waning is turned off by default
    assert( v.size()==num_params );
    
    // v_treatment_coverage.insert( v_treatment_coverage.begin(), NUMAC, 0.0 );
    tvp_treatment_coverage = new time_varying_prms_age<double>();
    v_treatment_success.insert( v_treatment_success.begin(), NUMAC, 1.0 );

    v_deathrate.insert( v_deathrate.begin(), NUMAC, 0.0 );

    v_pop_ac.insert( v_pop_ac.begin(), NUMAC, 0.0 );
    v_cumul_vac.insert( v_cumul_vac.begin(), NUMAC, 0.0 );

    v_rel_susc.insert( v_rel_susc.begin(), NUMAC, 0.0 );
    /* v_mixing_level.insert( v_mixing_level.begin(), NUMAC, 0.0 );
    v_mixing_level_postld.insert( v_mixing_level_postld.begin(), NUMAC, 0.0 ); */
    v_prob_E_A.insert( v_prob_E_A.begin(), NUMAC, 0.0 );
    v_prob_I2_H.insert(  v_prob_I2_H.begin(),  NUMAC, 0.0 );
    v_prob_I2_H_treatment.insert(  v_prob_I2_H_treatment.begin(),  NUMAC, 0.0 );
    v_prob_I4_D.insert(  v_prob_I4_D.begin(),  NUMAC, 0.0 );
    v_prob_HA4_D.insert( v_prob_HA4_D.begin(), NUMAC, 0.0 );
    v_prob_HA_CA.insert( v_prob_HA_CA.begin(), NUMAC, 0.0 );
    v_prob_V_D.insert( v_prob_V_D.begin(), NUMAC, 0.0 );
    v_prob_CA_D.insert( v_prob_CA_D.begin(), NUMAC, 0.0 );
    v_prob_CR_D.insert( v_prob_CR_D.begin(), NUMAC, 0.0 );
    v_prob_CA_V.insert( v_prob_CA_V.begin(), NUMAC, 0.0 );
    
    // assert( v_treatment_coverage.size() == NUMAC);
    assert( v_treatment_success.size() == NUMAC);
    assert( v_deathrate.size() == NUMAC);
    assert( v_pop_ac.size()==NUMAC );  
    assert( v_cumul_vac.size()==NUMAC );  
    assert( v_rel_susc.size()==NUMAC );  
    /* assert( v_mixing_level.size()==NUMAC );  
    assert( v_mixing_level_postld.size()==NUMAC );   */
    assert( v_prob_E_A.size()==NUMAC );
    assert( v_prob_I2_H.size() ==NUMAC );
    assert( v_prob_I2_H_treatment.size() ==NUMAC );
    assert( v_prob_I4_D.size() ==NUMAC );
    assert( v_prob_HA4_D.size()==NUMAC );
    assert( v_prob_HA_CA.size()==NUMAC );
    assert( v_prob_V_D.size()==NUMAC );
    assert( v_prob_CA_D.size()==NUMAC );
    assert( v_prob_CR_D.size()==NUMAC );
    assert( v_prob_CA_V.size()==NUMAC );
    
    v_betas.clear();
    v_betatimes.clear();
    assert( v_betas.size() == 0 );
    assert( v_betatimes.size() == 0 );
    
    earlymarch_highhosp_period = false;
    earlymarch_highhosp_factor = 1.0;
    earlymarch_highhosp_endday = -1.0;
    
    index_current_beta=-1;
    
}

// destructor
prms::~prms()
{
    delete tvp_contact_rate;
    delete tvp_contact_coeff;
    delete tvp_hosp_frac;
    delete tvp_vac_num;
    delete tvp_vac_efficacy;
    delete tvp_dev_icu_frac;
    delete tvp_dev_len_hospstay;
    delete tvp_incubation_period;
    delete tvp_infectious_period;
    delete tvp_death_home;
    delete tvp_treatment_coverage;
}


void prms::assign_new_beta( void )
{
    assert( v_betas.size() == v_betatimes.size() );
    assert( v_betas.size() > 0 );
    
    // if it's invalid, just set it to zero
    // this is what happens at initialization as well
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        index_current_beta = 0;
    }
    else // if it's valid, increment it if you can
    {
        if( index_current_beta != v_betas.size() - 1 ) // if it's NOT the last element, then increment it
        {
            index_current_beta++;
        }
    }
    
    // assign the possibly new beta value to the main parameter vector "v"
    v[i_beta] = v_betas[ index_current_beta ];
}


double prms::get_new_update_time( void )
{
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        return 1000000.0;
    }
    else // if the index is within range, return the next update time
    {
        if( index_current_beta != v_betas.size() - 1 ) // if the index does NOT point to the last element, then return the next element
        {
            return v_betatimes[ index_current_beta+1 ];
        }
        else // if it is the last element in the array, just return 1000000.0, as in the next update time is a million days from now
        {
            return 1000000.0;
        }
    }
    
    
}

void prms::apply_earlymarch_hosprates( void )
{
    for(int ac=0; ac<NUMAC; ac++)
    {
        v_prob_I2_H[ac] = 1.0 - ( (1.0-v_prob_I2_H[ac]) / earlymarch_highhosp_factor );
    }
    earlymarch_highhosp_period = true;
}

void prms::end_earlymarch_hosprates( void )
{
    for(int ac=0; ac<NUMAC; ac++)
    {
        v_prob_I2_H[ac] = 1.0 - ( (1.0-v_prob_I2_H[ac]) * earlymarch_highhosp_factor );
    }
    earlymarch_highhosp_period = false;
}


// ensure 0 <= v_prob_I2_H_treatment <= 1,  0 <= v_treatment_success <= 1, 0 <= v_treatment_coverage <= 1
// new prob_I2_H_treatment = old (current) prob_I2_H_treatment / (1 - success fraction)
void prms::calculate_prob_I2_H_treatment( void ){
    // = hospitalization fraction given treatment / (1 - success fraction)

    for (int ac = 0; ac < NUMAC; ac++)
    {
        // ensure 0 <= v_prob_I2_H_treatment <= 1
        if ( v_prob_I2_H_treatment[ac] < 0.0){
            v_prob_I2_H_treatment[ac] = 0.0;
        } else if ( v_prob_I2_H_treatment[ac] > 1.0){
            v_prob_I2_H_treatment[ac] = 1.0;
        } else {
            // do nothing, v_prob_I2_H_treatment = itself
        }

        // new prob_I2_H_treatment = old (current) prob_I2_H_treatment / (1 - success fraction)
        if ( v_treatment_success[ac] >= 1.0 )
        {
            v_treatment_success[ac] = 1.0;
            v_prob_I2_H_treatment[ac] = 0.0;
        } 
        else if ( v_treatment_success[ac] < 0.0 )
        {
            v_treatment_success[ac] = 0.0;
            // v_prob_I2_H_treatment = itself
        } 
        else {
            // v_treatment_success = itself
            v_prob_I2_H_treatment[ac] = v_prob_I2_H_treatment[ac] / (1.0 - v_treatment_success[ac]);
        }

        // ensure 0 <= treatment coverage <= 1
        /* if ( v_treatment_coverage[ac] < 0.0){
            v_treatment_coverage[ac] = 0.0;
        } else if ( v_treatment_coverage[ac] > 1.0){
            v_treatment_coverage[ac] = 1.0;
        } else {
            // do nothing, v_treatment_coverage = itself
        } */
    }
    
}