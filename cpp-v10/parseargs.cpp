#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "essentials.h"
#include "prms.h"
#include "assert.h"



// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double* yic;  
extern prms* ppc;
extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;
extern FILE* OutFile;
extern double G_CLO_TF;
extern int G_CLO_STEPS_PER_DAY;
extern string G_CLO_LOCATION;
extern double G_CLO_P_HOSP_TO_ICU;
extern bool G_B_DIAGNOSTIC_MODE;
extern bool G_B_CHECKPOP_MODE;
extern bool G_B_USESOCIALCONTACTMATRIX;
extern bool G_B_REMOVE_99COLUMNS;
extern bool G_B_BINARY_OUTPUT;
extern bool G_B_DEATHRATE_OUTPUT;
extern bool G_B_S_TO_R_OUTPUT;

extern double G_CLO_SYMP_FRAC;          //SHOULD BE DEPRECATED SOON
extern double G_CLO_SYMP_FRAC_EQUAL;
extern double G_B_SYMP_FRAC_DAVIES_20200616;
extern double G_B_SYMP_FRAC_SIMPLEAVERAGE;



extern double G_CLO_HOSPFRAC_YOUNG_DEV; //DEPRECATED
extern double G_CLO_HOSPFRAC_OLD_DEV;   //DEPRECATED
extern double G_CLO_HOSPFRAC_MID_DEV;   //DEPRECATED

/* extern double G_CLO_HOSPFRAC_10;
extern double G_CLO_HOSPFRAC_20;
extern double G_CLO_HOSPFRAC_30;
extern double G_CLO_HOSPFRAC_40;
extern double G_CLO_HOSPFRAC_50;
extern double G_CLO_HOSPFRAC_60;
extern double G_CLO_HOSPFRAC_70;
extern double G_CLO_HOSPFRAC_80; */

/* extern double G_CLO_MIXINGLEVEL[];
extern double G_CLO_MIXINGLEVEL_POSTLD[];
extern int G_CLO_FIRSTLOCKDOWN_ENDDAY;
extern bool G_CLO_POSTLD_MIXING_SET; */


// extern double G_CLO_DEATHPROB_HOME_60;
extern double G_CLO_DEATHPROB_HOME_70;
extern double G_CLO_DEATHPROB_HOME_80;

extern double G_CLO_DEATHPROB_POSTVENT;

extern double G_CLO_PROB_NONICU_DEATH_80;

extern double G_CLO_MIN_MIXINGLEVEL_00;
extern double G_CLO_MIN_MIXINGLEVEL_10;
extern double G_CLO_MIN_MIXINGLEVEL_20;
extern double G_CLO_MIN_MIXINGLEVEL_30;
extern double G_CLO_MIN_MIXINGLEVEL_40;
extern double G_CLO_MIN_MIXINGLEVEL_50;
extern double G_CLO_MIN_MIXINGLEVEL_60;
extern double G_CLO_MIN_MIXINGLEVEL_70;
extern double G_CLO_MIN_MIXINGLEVEL_80;

extern double G_CLO_RELSUSC_0;
extern double G_CLO_RELSUSC_10;
extern double G_CLO_RELSUSC_20;
extern double G_CLO_RELSUSC_30;
extern double G_CLO_RELSUSC_40;
extern double G_CLO_RELSUSC_50;
extern double G_CLO_RELSUSC_60;
extern double G_CLO_RELSUSC_70;
extern double G_CLO_RELSUSC_80;

extern double G_CLO_TIME_SYMPTOHOSP;
extern double G_CLO_SELFISOLATION_FACTOR;

extern double G_CLO_EARLYMARCH_HOSPRATE;     // scaling factor showing how much more likely hospitalization was in early March than later in the epidemic
extern double G_CLO_EARLYMARCH_ENDDAY;       // the day that the high-hosp rate ended (probably because the patient load became too high)




// extern double G_CLO_ICUFRAC_DEV;
// extern double G_CLO_ICUFRAC_DEV_SECONDPHASE;
// extern int G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY;
extern double G_CLO_PROB_ICU_TO_VENT;
extern double G_CLO_VENTDEATH_MID_DEV;
extern double G_CLO_VENTDEATH_70_DEV;
extern double G_CLO_VENTDEATH_80_DEV;
extern double G_CLO_MEANTIME_ON_VENT_SURV;
extern double G_CLO_RELATIVE_BETA_HOSP;

// extern double G_CLO_DEV_LEN_HOSPSTAY;

extern double G_C_COMM[NUMAC][NUMAC];

//  END  ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void PrintUsageModes();
void SetLocationData( string loc );

bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}

void PrintUsageModes()
{
    printf("\n\tUSAGE: ./odesim  outfilename   [-options] \n\n");
}


//
// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=2;
    
    // this counter is here to check if we have set any of the G_CLO_MIXINGLEVEL_POSTLD parameters
    // this number has to be either 0 or 8; anything else leads to an error
    int num_postld_contact_rates_set = 0;

    if( argc<2 )
    { 
        PrintUsageModes(); 
	exit(-1);
    }
	
    if( argv[1][0]=='n' && argv[1][1]=='o' && argv[1][2]=='n' && argv[1][3]=='e' && argv[1][4]==0 )
    {
        //fprintf(stderr, "\n\n\tnot printing to outfile\n\n");
    }
    else 
    {
        if( argv[1][0]=='-'){
            fprintf(stderr, "You probably meant to use '%s' as an option, not the output file.\n", argv[1]);
            exit(-1);
        }
        OutFile = fopen( argv[1], "w" );
    }

    string str;

    i=start;
    
    // read in options
    while(i<argc)
    {
	str = argv[i];
        
        // ### 1 ### IF BLOCK FOR BETA
        if( str == "-beta" )
        {
            ppc->v_betas.clear();
            ppc->v_betatimes.clear();
            i++;
            bool bFirstBetaPushedBack = false;
            
            //BEGIN LOOPING THROUGH THE LIST OF BETAS
            while(i<argc)
            {
                string s( argv[i] );
                if( isFloat(s) ) // if the current string is a floating point number, write it into the v_betas array
                {
                    // if the command line argument is <0, just set it back to zero
                    double d = atof( argv[i] );
                    if( d < 0.0 ) d = 0.0;
                    
                    ppc->v_betas.push_back( d );
                    
                    // if it's the first beta value being read in, then push it back twice, as it will be used as 
                    // the Jan 1 to March 1 beta (when there were almost no cases) as well as the beta for the 
                    // first time period we are invesigating
                    if(!bFirstBetaPushedBack) {ppc->v_betas.push_back( d );bFirstBetaPushedBack=true;}
                    
                    // increment and move on in this sub-loop
                    i++;
                }
                else
                {
                    // if the current string is NOT a float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }
                 
            } 
            //END OF LOOPING THROUGH THE LIST OF BETAS
             
            // MAKE SURE AT LEAST ONE BETA VALUE WAS READ IN 
            if( ppc->v_betas.size() == 0 )
            {
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-beta' option must be followed by at least one floating point number.\n\n");
                exit(-1);
            }
            else
            {
                double total_time_period = G_CLO_TF - 60.0;
                double time_step = total_time_period / ( (double) ppc->v_betas.size() - 1 );
                
                ppc->v_betatimes.push_back(0.0);
                
                for(int jj=1; jj<ppc->v_betas.size(); jj++)
                {
                    ppc->v_betatimes.push_back( 60.0 + time_step*((double)(jj-1)) );
                }
                
            }
                         
        }
        //
        // ### - ### END IF BLOCK FOR BETA
        //
        else if( str == "-binary-output" )  {           G_B_BINARY_OUTPUT = true;                   }   
        else if( str == "-checkpop" )
        {
            G_B_CHECKPOP_MODE = true;
        }
        /*else if( str == "-dev-hosp-young" )
        {
            G_CLO_HOSPFRAC_YOUNG_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_YOUNG_DEV > 3.0 ) G_CLO_HOSPFRAC_YOUNG_DEV = 3.0;
            if( G_CLO_HOSPFRAC_YOUNG_DEV < 0.0 ) G_CLO_HOSPFRAC_YOUNG_DEV = 0.0;
        }
        else if( str == "-dev-hosp-mid" )
        {
            G_CLO_HOSPFRAC_MID_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_MID_DEV > 3.0 ) G_CLO_HOSPFRAC_MID_DEV = 3.0;
            if( G_CLO_HOSPFRAC_MID_DEV < 0.0 ) G_CLO_HOSPFRAC_MID_DEV = 0.0;
        }
        else if( str == "-dev-hosp-old" )
        {
            G_CLO_HOSPFRAC_OLD_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_OLD_DEV > 1.3 ) G_CLO_HOSPFRAC_OLD_DEV = 1.3;
            if( G_CLO_HOSPFRAC_OLD_DEV < 0.0 ) G_CLO_HOSPFRAC_OLD_DEV = 0.0;
        }*/
        /* else if( str == "-contact-rate-10" )        {     G_CLO_MIXINGLEVEL[1] = atof( argv[++i] );        }
        else if( str == "-contact-rate-20" )        {     G_CLO_MIXINGLEVEL[2] = atof( argv[++i] );        }
        else if( str == "-contact-rate-30" )        {     G_CLO_MIXINGLEVEL[3] = atof( argv[++i] );        }
        else if( str == "-contact-rate-40" )        {     G_CLO_MIXINGLEVEL[4] = atof( argv[++i] );        }
        else if( str == "-contact-rate-50" )        {     G_CLO_MIXINGLEVEL[5] = atof( argv[++i] );        }
        else if( str == "-contact-rate-60" )        {     G_CLO_MIXINGLEVEL[6] = atof( argv[++i] );        }
        else if( str == "-contact-rate-70" )        {     G_CLO_MIXINGLEVEL[7] = atof( argv[++i] );        }
        else if( str == "-contact-rate-80" )        {     G_CLO_MIXINGLEVEL[8] = atof( argv[++i] );        }
        else if( str == "-contact-rate-postld-10" ) {     G_CLO_MIXINGLEVEL_POSTLD[1] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-20" ) {     G_CLO_MIXINGLEVEL_POSTLD[2] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-30" ) {     G_CLO_MIXINGLEVEL_POSTLD[3] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-40" ) {     G_CLO_MIXINGLEVEL_POSTLD[4] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-50" ) {     G_CLO_MIXINGLEVEL_POSTLD[5] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-60" ) {     G_CLO_MIXINGLEVEL_POSTLD[6] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-70" ) {     G_CLO_MIXINGLEVEL_POSTLD[7] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        }
        else if( str == "-contact-rate-postld-80" ) {     G_CLO_MIXINGLEVEL_POSTLD[8] = atof( argv[++i] ); num_postld_contact_rates_set++; G_CLO_POSTLD_MIXING_SET=true;        } */
        else if( str == "-tv-contact-rate-10" )   {           
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-contact-rate-10
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_contact_rate->append_data(0, 1.0); // push default value 1.0 to the end of v_data[0] (age group 0-9) vector of tvp_contact_rate 
                    ppc->tvp_contact_rate->append_data(1, atof( argv[i] ) ); // push read value to the end of v_data[1] (age group 10-19) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-contact-rate-10 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(1) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-10' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-20" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(2, atof( argv[i] ) ); // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-20 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(2) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-20' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-30" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(3, atof( argv[i] ) ); // push the read value to the end of v_data[3] (age group 30-39) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-30 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(3) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-30' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-40" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(4, atof( argv[i] ) ); // push the read value to the end of v_data[4] (age group 40-49) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-40 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(4) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-40' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-50" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(5, atof( argv[i] ) ); // push the read value to the end of v_data[5] (age group 50-59) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-50 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(5) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-50' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-60" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(6, atof( argv[i] ) ); // push the read value to the end of v_data[6] (age group 60-69) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-60 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(6) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-60' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-70" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(7, atof( argv[i] ) ); // push the read value to the end of v_data[7] (age group 70-79) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-70 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(7) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-70' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-rate-80" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_rate->append_data(8, atof( argv[i] ) ); // push the read value to the end of v_data[8] (age group 80+) vector of tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-rate-40 VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_data(8) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-80' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        // else if( str == "-firstlockdown-endday" )   {     G_CLO_FIRSTLOCKDOWN_ENDDAY = atoi( argv[++i] );       }
        else if( str == "-tv-contact-rate-beginday" ){        
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-contact-rate-beginday
            while(i<argc)
            {
                // int d = -1;
                // if the current string is convertible to non-negative integer 
                if( isdigit( argv[i][0] ) ){
                    // d = atoi( argv[i] );
                    ppc->tvp_contact_rate->append_begin_day( atoi(argv[i]) ); // push to the end of v_begin_days vector in tvp_contact_rate 
                    i++; // increment and move on in this sub-loop
                    // printf("%d\t", d);
                } else {
                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->tvp_contact_rate->empty_begin_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-beginday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        //// time-varying contact-rate
        else if( str == "-tv-contact-rate-endday" ){        
            i++; // increase index to get value after the str
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative integer
                    ppc->tvp_contact_rate->append_end_day( atoi(argv[i]) ); 
                    i++; 
                } else {
                    i--;   
                    break;
                }   
            }  //END OF WHILE LOOP
            if( ppc->tvp_contact_rate->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-rate-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }


        else if( str == "-tv-contact-coeff-00" )   {           
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-contact-coeff-10
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_contact_coeff->append_data(0, atof( argv[i] ) ); // push read value to the end of v_data[0] (age group 0-9) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-contact-coeff-00 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(0) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-00' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-10" )   {           
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-contact-coeff-10
            while(i<argc){
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_contact_coeff->append_data(1, atof( argv[i] ) ); // push read value to the end of v_data[1] (age group 10-19) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                }
                else {
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-contact-coeff-10 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(1) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-10' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-20" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(2, atof( argv[i] ) ); // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-20 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(2) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-20' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-30" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(3, atof( argv[i] ) ); // push the read value to the end of v_data[3] (age group 30-39) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-30 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(3) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-30' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-40" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(4, atof( argv[i] ) ); // push the read value to the end of v_data[4] (age group 40-49) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-40 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(4) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-40' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-50" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(5, atof( argv[i] ) ); // push the read value to the end of v_data[5] (age group 50-59) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-50 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(5) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-50' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-60" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(6, atof( argv[i] ) ); // push the read value to the end of v_data[6] (age group 60-69) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-60 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(6) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-60' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-70" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(7, atof( argv[i] ) ); // push the read value to the end of v_data[7] (age group 70-79) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-70 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(7) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-70' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-contact-coeff-80" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_contact_coeff->append_data(8, atof( argv[i] ) ); // push the read value to the end of v_data[8] (age group 80+) vector of tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-contact-coeff-40 VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_data(8) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-80' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        // else if( str == "-firstlockdown-endday" )   {     G_CLO_FIRSTLOCKDOWN_ENDDAY = atoi( argv[++i] );       }
        else if( str == "-tv-contact-coeff-beginday" ){        
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-contact-coeff-beginday
            while(i<argc)
            {
                // int d = -1;
                // if the current string is convertible to non-negative integer 
                if( isdigit( argv[i][0] ) ){
                    // d = atoi( argv[i] );
                    ppc->tvp_contact_coeff->append_begin_day( atoi(argv[i]) ); // push to the end of v_begin_days vector in tvp_contact_coeff 
                    i++; // increment and move on in this sub-loop
                    // printf("%d\t", d);
                } else {
                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->tvp_contact_coeff->empty_begin_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-beginday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        //// time-varying contact-coeff
        else if( str == "-tv-contact-coeff-endday" ){        
            i++; // increase index to get value after the str
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative integer
                    ppc->tvp_contact_coeff->append_end_day( atoi(argv[i]) ); 
                    i++; 
                } else {
                    i--;   
                    break;
                }   
            }  //END OF WHILE LOOP
            if( ppc->tvp_contact_coeff->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-contact-coeff-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        /* else if( str == "-dev-icu-frac" )
        {
            G_CLO_ICUFRAC_DEV = atof( argv[++i] );
            if( G_CLO_ICUFRAC_DEV > 2.0 ) G_CLO_ICUFRAC_DEV = 2.0;
            if( G_CLO_ICUFRAC_DEV < 0.0 ) G_CLO_ICUFRAC_DEV = 0.0;
        }
        else if( str == "-dev-icu-frac-phase2" )
        {
            G_CLO_ICUFRAC_DEV_SECONDPHASE = atof( argv[++i] );
            if( G_CLO_ICUFRAC_DEV_SECONDPHASE > 2.0 ) G_CLO_ICUFRAC_DEV_SECONDPHASE = 2.0;
            if( G_CLO_ICUFRAC_DEV_SECONDPHASE < 0.0 ) G_CLO_ICUFRAC_DEV_SECONDPHASE = 0.0;
        }
        else if( str == "-dev-icu-frac-phase2beginday" )   {     G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY = atoi( argv[++i] );       } */
        else if( str == "-tv-dev-icu-frac")
        {
            i++;
            double dif_v;
            //BEGIN LOOPING THROUGH THE LIST OF tv-dev-icu-frac
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    dif_v = atof( argv[i]);
                    dif_v = (dif_v < 0.0) ? 0.0 : dif_v;
                    dif_v = (dif_v > 2.0) ? 2.0 : dif_v;
                    ppc->tvp_dev_icu_frac->append_data( dif_v ); // push read value to the end of v_data
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-dev-icu-frac VALUE WAS READ IN 
            if( ppc->tvp_dev_icu_frac->empty_data() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-dev-icu-frac' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-dev-icu-frac-beginday" )
        {
            i++;
            while(i<argc)
            {
                // if the current string is convertible to non-negative integer 
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_dev_icu_frac->append_begin_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else {
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->tvp_dev_icu_frac->empty_begin_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-dev-icu-frac-beginday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-dev-icu-frac-endday" )
        {
            i++;
            while(i<argc)
            {
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_dev_icu_frac->append_end_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else {
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-endday VALUE WAS READ IN 
            if( ppc->tvp_dev_icu_frac->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-dev-icu-frac-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        /* else if( str == "-dev-len-hospstay" )
        {
            G_CLO_DEV_LEN_HOSPSTAY = atof( argv[++i] );
            if( G_CLO_DEV_LEN_HOSPSTAY > 2.5 ) G_CLO_DEV_LEN_HOSPSTAY = 2.5;
            if( G_CLO_DEV_LEN_HOSPSTAY < 0.0 ) G_CLO_DEV_LEN_HOSPSTAY = 0.0;
        } */
        else if( str == "-tv-dev-len-hospstay")
        {
            i++;
            double dif_v;
            //BEGIN LOOPING THROUGH THE LIST OF tv-dev-len-hospstay
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    dif_v = atof( argv[i]);
                    dif_v = (dif_v < 0.0) ? 0.0 : dif_v;
                    dif_v = (dif_v > 2.5) ? 2.5 : dif_v;
                    ppc->tvp_dev_len_hospstay->append_data( dif_v ); // push read value to the end of v_data
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-dev-icu-frac VALUE WAS READ IN 
            if( ppc->tvp_dev_len_hospstay->empty_data() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-dev-len-hospstay' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-dev-len-hospstay-beginday" )
        {
            i++;
            while(i<argc)
            {
                // if the current string is convertible to non-negative integer 
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_dev_len_hospstay->append_begin_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else {
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->tvp_dev_len_hospstay->empty_begin_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-dev-len-hospstay-beginday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-dev-len-hospstay-endday" )
        {
            i++;
            while(i<argc)
            {
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_dev_len_hospstay->append_end_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else {
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-endday VALUE WAS READ IN 
            if( ppc->tvp_dev_len_hospstay->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-dev-len-hospstay-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        else if( str == "-tv-incub-period")
        {
            i++;
            double incub_v;
            //BEGIN LOOPING THROUGH THE LIST OF tv-dev-len-hospstay
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    incub_v = atof( argv[i]);

                    if (incub_v >=0 && incub_v <= 10)
                    {
                        ppc->tvp_incubation_period->append_data( incub_v ); // push read value to the end of v_data
                    } else {
                        fprintf(stderr,"\n\n\t'-tv-incub-period' option must be within [0, 10].\n\n");
                    }
                    
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-incub-period VALUE WAS READ IN 
            if( ppc->tvp_incubation_period->empty_data() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-incub-period' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-incub-period-endday" )
        {
            i++;
            while(i<argc)
            {
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_incubation_period->append_end_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else {
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-incub-period-endday VALUE WAS READ IN 
            if( ppc->tvp_incubation_period->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-incub-period-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        else if( str == "-tv-infect-period")
        {
            i++;
            double infect_v;
            //BEGIN LOOPING THROUGH THE LIST OF tv-dev-len-hospstay
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    infect_v = atof( argv[i]);

                    if (infect_v >0 ){
                        ppc->tvp_infectious_period->append_data( infect_v ); // push read value to the end of v_data
                    } else {
                        fprintf(stderr,"\n\n\t'-tv-infect-period' value must be positive.\n\n");
                    }
                    
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-infect-period VALUE WAS READ IN 
            if( ppc->tvp_infectious_period->empty_data() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-infect-period' option must be followed by at least one positive number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-infect-period-endday" )
        {
            i++;
            while(i<argc)
            {
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_infectious_period->append_end_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else {
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-incub-period-endday VALUE WAS READ IN 
            if( ppc->tvp_infectious_period->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-infect-period-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        else if( str == "-dev-ventdeath-mid" )
        {
            G_CLO_VENTDEATH_MID_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_MID_DEV > 1.7 ) G_CLO_VENTDEATH_MID_DEV = 1.7;
            if( G_CLO_VENTDEATH_MID_DEV < 0.0 ) G_CLO_VENTDEATH_MID_DEV = 0.0;
        }
        else if( str == "-dev-ventdeath-70" )
        {
            G_CLO_VENTDEATH_70_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_70_DEV > 1.4 ) G_CLO_VENTDEATH_70_DEV = 1.4;
            if( G_CLO_VENTDEATH_70_DEV < 0.8 ) G_CLO_VENTDEATH_70_DEV = 0.8;
        }
        else if( str == "-dev-ventdeath-80" )
        {
            G_CLO_VENTDEATH_80_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_80_DEV > 1.1 ) G_CLO_VENTDEATH_80_DEV = 1.1;
            if( G_CLO_VENTDEATH_80_DEV < 0.7 ) G_CLO_VENTDEATH_80_DEV = 0.7;
        }
        // else if( str == "-death-prob-home-60" )     {     G_CLO_DEATHPROB_HOME_60 = atof( argv[++i] );          }
        // else if( str == "-death-prob-home-70" )     {     G_CLO_DEATHPROB_HOME_70 = atof( argv[++i] );          }
        // else if( str == "-death-prob-home-80" )     {     G_CLO_DEATHPROB_HOME_80 = atof( argv[++i] );          } 
        
        else if( str == "-tv-death-prob-home-00" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_data(0, atof( argv[i] ) );   i++;
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(0) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-00' option must be followed by at least one non-negative number.\n\n"); exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-10" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_death_home->append_data(1, atof( argv[i] ) );    i++;
                } else { i--;   break; }   
            } //END OF WHILE LOOP
            if( ppc->tvp_death_home->empty_data(1) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-10' option must be followed by at least one non-negative number.\n\n"); exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-20" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_data(2, atof( argv[i] ) );    i++;
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(2) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-20' option must be followed by at least one non-negative number.\n\n");    exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-30" )   {           
            i++;
            while(i<argc){
                if( isdigit( 
                    argv[i][0] ) ){ ppc->tvp_death_home->append_data(3, atof( argv[i] ) );    i++;
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(3) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-30' option must be followed by at least one non-negative number.\n\n");   exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-40" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_death_home->append_data(4, atof( argv[i] ) );    i++;
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(4) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-40' option must be followed by at least one non-negative number.\n\n");    exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-50" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_data(5, atof( argv[i] ) );   i++; 
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(5) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-50' option must be followed by at least one non-negative number.\n\n");    exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-60" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_data(6, atof( argv[i] ) );   i++; 
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(6) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-60' option must be followed by at least one non-negative number.\n\n");    exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-70" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_data(7, atof( argv[i] ) );   i++;
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(7) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-70' option must be followed by at least one non-negative number.\n\n");   exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-80" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_data(8, atof( argv[i] ) );  i++;
                } else { i--; break; }   
            } 
            if( ppc->tvp_death_home->empty_data(8) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-80' option must be followed by at least one non-negative number.\n\n");    exit(-1);
            }
        }
        else if( str == "-tv-death-prob-home-endday" ){        
            i++; 
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ 
                    ppc->tvp_death_home->append_end_day( atoi(argv[i]) );   i++; 
                } else { i--;   break; }   
            }  //END OF WHILE LOOP
            if( ppc->tvp_death_home->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-death-prob-home-endday' option must be followed by at least one non-negative integer.\n\n"); exit(-1);
            }
        }
        

        else if( str == "-death-prob-nonicu-80" )   {     G_CLO_PROB_NONICU_DEATH_80 = atof( argv[++i] );       }
        else if( str == "-death-prob-postvent" )    {     G_CLO_DEATHPROB_POSTVENT   = atof( argv[++i] );       }
        
        else if( str == "-death-rate-output" )      {     G_B_DEATHRATE_OUTPUT = true;                          }
        else if( str == "-s-to-r-output" )          {     G_B_S_TO_R_OUTPUT = true;                             }
        
        else if( str == "-diag" )                   {     G_B_DIAGNOSTIC_MODE = true;                           }
        else if( str == "-earlymarch-hosp-factor" ) {     G_CLO_EARLYMARCH_HOSPRATE = atof( argv[++i] );        }
        
        //// time-varying hosp-frac
        else if( str == "-tv-hosp-frac-beginday" ){        
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-hosp-frac-beginday
            while(i<argc)
            {
                // int d = -1;
                // if the current string is convertible to non-negative integer 
                if( isdigit( argv[i][0] ) ){
                    // d = atoi( argv[i] );
                    ppc->tvp_hosp_frac->append_begin_day( atoi(argv[i]) ); // push to the end of v_begin_days vector in tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                    // printf("%d\t", d);
                } else {
                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_begin_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-beginday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-endday" ){        
            i++; // increase index to get value after the str
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative integer
                    ppc->tvp_hosp_frac->append_end_day( atoi(argv[i]) ); 
                    i++; 
                } else {
                    i--;   
                    break;
                }   
            }  //END OF WHILE LOOP
            if( ppc->tvp_hosp_frac->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        else if( str == "-introday" )       {           G_CLO_INTRODUCTION_TIME = atof( argv[++i] );        }
        else if( str == "-loc" )            {           G_CLO_LOCATION = argv[++i];                         }

        
        else if( str == "-tv-hosp-frac-10" )   {           
            i++;
            double hf_v;
            //BEGIN LOOPING THROUGH THE LIST OF tv-hosp-frac-10
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(0, hf_v ); // push default value 1.0 to the end of v_data[0] (age group 0-9) vector of tvp_hosp_frac 
                    ppc->tvp_hosp_frac->append_data(1, hf_v ); // push read value to the end of v_data[1] (age group 10-19) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-hosp-frac-10 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(1) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-10' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-20" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(2, hf_v ); // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-20 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(2) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-20' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-30" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(3, hf_v ); // push the read value to the end of v_data[3] (age group 30-39) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-30 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(3) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-30' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-40" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(4, hf_v ); // push the read value to the end of v_data[4] (age group 40-49) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-40 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(4) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-40' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-50" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(5, hf_v ); // push the read value to the end of v_data[5] (age group 50-59) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-50 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(5) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-50' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-60" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(6, hf_v ); // push the read value to the end of v_data[6] (age group 60-69) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-60 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(6) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-60' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-70" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(7, hf_v ); // push the read value to the end of v_data[7] (age group 70-79) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-70 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(7) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-70' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-hosp-frac-80" )   {           
            i++;
            double hf_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    hf_v = atof( argv[i] );
                    hf_v = (hf_v < 0.0) ? 0.0 : hf_v;
                    hf_v = (hf_v > 1.0) ? 1.0 : hf_v;
                    ppc->tvp_hosp_frac->append_data(8, hf_v ); // push the read value to the end of v_data[8] (age group 80+) vector of tvp_hosp_frac 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-hosp-frac-40 VALUE WAS READ IN 
            if( ppc->tvp_hosp_frac->empty_data(8) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-hosp-frac-80' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }

        // vaccinations
        else if( str == "-tv-vaccinees-00" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(0, atof( argv[i] ) ); // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-20 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(0) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-00' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-10" )   {           
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-vaccinees-10
            while(i<argc)
            {
                // if the current string is convertible to non-negative float
                if( isdigit( argv[i][0] ) ){
                    // ppc->tvp_vac_num->append_data(0, 0.0); // push default value 0.0 to the end of v_data[0] (age group 0-9) vector of tvp_vac_num 
                    ppc->tvp_vac_num->append_data(1, atof( argv[i] ) ); // push read value to the end of v_data[1] (age group 10-19) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                }
                else {
                    // if the current string is NOT a convertible to float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }   
            } //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE tv-vaccinees-10 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(1) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-10' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-20" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(2, atof( argv[i] ) ); // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-20 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(2) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-20' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-30" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(3, atof( argv[i] ) ); // push the read value to the end of v_data[3] (age group 30-39) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-30 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(3) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-30' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-40" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(4, atof( argv[i] ) ); // push the read value to the end of v_data[4] (age group 40-49) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-40 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(4) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-40' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-50" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(5, atof( argv[i] ) ); // push the read value to the end of v_data[5] (age group 50-59) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-50 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(5) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-50' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-60" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(6, atof( argv[i] ) ); // push the read value to the end of v_data[6] (age group 60-69) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-60 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(6) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-60' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-70" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(7, atof( argv[i] ) ); // push the read value to the end of v_data[7] (age group 70-79) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-70 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(7) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-70' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-80" )   {           
            i++;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    ppc->tvp_vac_num->append_data(8, atof( argv[i] ) ); // push the read value to the end of v_data[8] (age group 80+) vector of tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-40 VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_data(8) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-80' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-beginday" ){        
            i++;
            //BEGIN LOOPING THROUGH THE LIST OF tv-vaccinees-beginday
            while(i<argc)
            {
                // int d = -1;
                // if the current string is convertible to non-negative integer 
                if( isdigit( argv[i][0] ) ){
                    // d = atoi( argv[i] );
                    ppc->tvp_vac_num->append_begin_day( atoi(argv[i]) ); // push to the end of v_begin_days vector in tvp_vac_num 
                    i++; // increment and move on in this sub-loop
                    // printf("%d\t", d);
                } else {
                    // if the current string is NOT an integer, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;  
                    break;
                }   
            }  //END OF WHILE LOOP
            // MAKE SURE AT LEAST ONE vac1-phase1-beginday VALUE WAS READ IN 
            if( ppc->tvp_vac_num->empty_begin_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-beginday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-vaccinees-endday" ){        
            i++; // increase index to get value after the str
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative integer
                    ppc->tvp_vac_num->append_end_day( atoi(argv[i]) ); 
                    i++; 
                } else {
                    i--;   
                    break;
                }   
            }  //END OF WHILE LOOP
            if( ppc->tvp_vac_num->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vaccinees-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        // vaccine efficacy
        else if( str == "-tv-vac-efficacy")
        {
            i++;
            double eff_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){
                    eff_v = atof( argv[i]);
                    eff_v = (eff_v < 0.0) ? 0.0 : eff_v;
                    eff_v = (eff_v > 1.0) ? 1.0 : eff_v;
                    ppc->tvp_vac_efficacy->append_data( eff_v ); // push read value to the end of v_data
                    i++; // increment and move on in this sub-loop
                }
                else { i--; break; }   
            } 
            if( ppc->tvp_vac_efficacy->empty_data() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vac-efficacy' option must be followed by at least one number in the range [0, 1].\n\n"); exit(-1);
            }
        }
        else if( str == "-tv-vac-efficacy-endday" )
        {
            i++;
            while(i<argc)
            {
                if( isdigit( argv[i][0] ) ){
                    ppc->tvp_vac_efficacy->append_end_day( atoi(argv[i]) ); // push to the end of v_begin_days vector
                    i++; // increment and move on in this sub-loop
                } else { i--;  break; }   
            }  
            if( ppc->tvp_vac_efficacy->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-vac-efficacy-endday' option must be followed by at least one non-negative integer.\n\n");  exit(-1);
            }
        }

        else if( str == "-no-vaccine-wastage" ){
            ppc->vaccine_wastage_flag = false;
        }
        else if( str == "-rr"){
            ppc->v[ i_avg_reporting_rate ] = atof(argv[++i]);
        }
        
        /* else if( str == "-treat-cover-00" ) {       ppc->v_treatment_coverage[0] = atof( argv[++i] );       }
        else if( str == "-treat-cover-10" ) {       ppc->v_treatment_coverage[1] = atof( argv[++i] );       }
        else if( str == "-treat-cover-20" ) {       ppc->v_treatment_coverage[2] = atof( argv[++i] );       }
        else if( str == "-treat-cover-30" ) {       ppc->v_treatment_coverage[3] = atof( argv[++i] );       }
        else if( str == "-treat-cover-40" ) {       ppc->v_treatment_coverage[4] = atof( argv[++i] );       }
        else if( str == "-treat-cover-50" ) {       ppc->v_treatment_coverage[5] = atof( argv[++i] );       }
        else if( str == "-treat-cover-60" ) {       ppc->v_treatment_coverage[6] = atof( argv[++i] );       }
        else if( str == "-treat-cover-70" ) {       ppc->v_treatment_coverage[7] = atof( argv[++i] );       }
        else if( str == "-treat-cover-80" ) {       ppc->v_treatment_coverage[8] = atof( argv[++i] );       } */

        else if( str == "-tv-treat-cov-00" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(0, tc_v ) ; // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            // MAKE SURE AT LEAST ONE tv-vaccinees-20 VALUE WAS READ IN 
            if( ppc->tvp_treatment_coverage->empty_data(0) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-00' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-10" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(1, tc_v ) ; // push the read value to the end of v_data[1] (age group 10-19) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(1) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-10' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-20" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(2, tc_v ) ; // push the read value to the end of v_data[2] (age group 20-29) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(2) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-20' option must be followed by at least one non-negative number.\n\n"); exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-30" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(3, tc_v ) ;// push the read value to the end of v_data[3] (age group 30-39) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(3) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-30' option must be followed by at least one non-negative number.\n\n"); exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-40" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(4, tc_v ) ; // push the read value to the end of v_data[4] (age group 40-49) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(4) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-40' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-50" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(5, tc_v ) ; // push the read value to the end of v_data[5] (age group 50-59) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(5) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-50' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-60" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(6, tc_v ) ; // push the read value to the end of v_data[6] (age group 60-69) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(6) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-60' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-70" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(7, tc_v ) ; // push the read value to the end of v_data[7] (age group 70-79) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(7) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-70' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-80" )   {           
            i++;
            double tc_v;
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative float
                    tc_v = atof( argv[i]);
                    if (tc_v < 0.0){ tc_v = 0.0;}
                    if (tc_v > 1.0){ tc_v = 1.0;}
                    ppc->tvp_treatment_coverage->append_data(8, tc_v ) ; // push the read value to the end of v_data[8] (age group 80+) vector of tvp_treatment_coverage 
                    i++; // increment and move on in this sub-loop
                } else { i--; break; }   
            } 
            if( ppc->tvp_treatment_coverage->empty_data(8) ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-80' option must be followed by at least one non-negative number.\n\n");
                exit(-1);
            }
        }
        else if( str == "-tv-treat-cov-endday" ){        
            i++; // increase index to get value after the str
            while(i<argc){
                if( isdigit( argv[i][0] ) ){ // if the current string is convertible to non-negative integer
                    ppc->tvp_treatment_coverage->append_end_day( atoi(argv[i]) ); 
                    i++; 
                } else { i--;   break; }   
            }  //END OF WHILE LOOP
            if( ppc->tvp_treatment_coverage->empty_end_days() ){
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-tv-treat-cov-endday' option must be followed by at least one non-negative integer.\n\n");
                exit(-1);
            }
        }

        else if( str == "-treat-success-00" ) {       ppc->v_treatment_success[0] = atof( argv[++i] );       }
        else if( str == "-treat-success-10" ) {       ppc->v_treatment_success[1] = atof( argv[++i] );       }
        else if( str == "-treat-success-20" ) {       ppc->v_treatment_success[2] = atof( argv[++i] );       }
        else if( str == "-treat-success-30" ) {       ppc->v_treatment_success[3] = atof( argv[++i] );       }
        else if( str == "-treat-success-40" ) {       ppc->v_treatment_success[4] = atof( argv[++i] );       }
        else if( str == "-treat-success-50" ) {       ppc->v_treatment_success[5] = atof( argv[++i] );       }
        else if( str == "-treat-success-60" ) {       ppc->v_treatment_success[6] = atof( argv[++i] );       }
        else if( str == "-treat-success-70" ) {       ppc->v_treatment_success[7] = atof( argv[++i] );       }
        else if( str == "-treat-success-80" ) {       ppc->v_treatment_success[8] = atof( argv[++i] );       }

        else if( str == "-treat-hosp-frac-00" ) {       ppc->v_prob_I2_H_treatment[0] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-10" ) {       ppc->v_prob_I2_H_treatment[1] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-20" ) {       ppc->v_prob_I2_H_treatment[2] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-30" ) {       ppc->v_prob_I2_H_treatment[3] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-40" ) {       ppc->v_prob_I2_H_treatment[4] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-50" ) {       ppc->v_prob_I2_H_treatment[5] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-60" ) {       ppc->v_prob_I2_H_treatment[6] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-70" ) {       ppc->v_prob_I2_H_treatment[7] = atof( argv[++i] );       }
        else if( str == "-treat-hosp-frac-80" ) {       ppc->v_prob_I2_H_treatment[8] = atof( argv[++i] );       }

        else if( str == "-treat-beginday" ) {       ppc->v[i_treatment_beginday] = atof( argv[++i] );       }

        else if( str == "-len-immunity" ) {       ppc->v[i_len_immunity] = atof( argv[++i] );       }
        
        else if( str == "-mean-time-vent" ) {           G_CLO_MEANTIME_ON_VENT_SURV = atof( argv[++i] );    }        
        else if( str == "-min-mixinglevel-00" ) {       G_CLO_MIN_MIXINGLEVEL_00 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-10" ) {       G_CLO_MIN_MIXINGLEVEL_10 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-20" ) {       G_CLO_MIN_MIXINGLEVEL_20 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-30" ) {       G_CLO_MIN_MIXINGLEVEL_30 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-40" ) {       G_CLO_MIN_MIXINGLEVEL_40 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-50" ) {       G_CLO_MIN_MIXINGLEVEL_50 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-60" ) {       G_CLO_MIN_MIXINGLEVEL_60 = atof( argv[++i] );       }
        else if( str == "-min-mixinglevel-70" ) {       G_CLO_MIN_MIXINGLEVEL_70 = atof( argv[++i] );       }        
        else if( str == "-min-mixinglevel-80" ) {       G_CLO_MIN_MIXINGLEVEL_80 = atof( argv[++i] );       }
        else if( str == "-printIndices" )
        {
            printf("STARTS %d\n", 0 );     // start index of E
            printf("NUMAC %d\n", NUMAC );       // number age groups
            printf("NUME %d\n", NUME );         // number of Exposed-class
            printf("STARTE %d\n", STARTE );     // start index of E
            printf("NUMA %d\n", NUMA );         // Asymptomatic-class
            printf("STARTA %d\n", STARTA );
            printf("NUMI %d\n", NUMI );         // I (symptomatic)-class
            printf("STARTI %d\n", STARTI );
            printf("NUMHA %d\n", NUMHA );       // Hospotalized Acute-class
            printf("STARTHA %d\n", STARTHA );
            printf("STARTCA %d\n", STARTCA );   // Critical care (ICU) Acute-class
            printf("NUMV %d\n", NUMV );         // mechanical Ventilator-class
            printf("STARTV %d\n", STARTV );
            printf("STARTCR %d\n", STARTCR );   // Critical care (ICU) Recovering-class
            printf("STARTHR %d\n", STARTHR );   // Hospitalized Recovering-class
            printf("STARTD %d\n", STARTD );             // Died-at-home class (age-stratified)
            printf("STARTDHOSP %d\n", STARTDHOSP );     // Died-in-hospital class (age-stratified)
            printf("STARTR %d\n", STARTR );             // Recovered-at-home-class (age-stratified)
            printf("STARTRHOSP %d\n", STARTRHOSP );         // Recovered-in-hospital-and-discharged class (age-stratified)
            printf("STARTVAX %d\n", STARTVAX );         // Vaccinated class
            printf("STARTJ %d\n", STARTJ );   // Cumulative I-class (age-stratified)
            printf("STARTK %d\n", STARTK );   // Cumulative hospitalization incidence-class (age-stratified)

            printf("DIMENSION %d\n", DIMENSION );

            exit(0);
        }
        else if( str == "-prob-icu-vent" )          {           G_CLO_PROB_ICU_TO_VENT     = atof( argv[++i] ); }
        
        else if( str == "-remove99" )               {           G_B_REMOVE_99COLUMNS = true;                    }   
        else if( str == "-version" ) 
        {
            printf("%d\n", 5);
            exit(0);
        }
        else if( str == "-rel-beta-hosp" )
        {
            G_CLO_RELATIVE_BETA_HOSP = atof( argv[++i] );
            if( G_CLO_RELATIVE_BETA_HOSP > 1.0 ) G_CLO_RELATIVE_BETA_HOSP = 1.0;
            if( G_CLO_RELATIVE_BETA_HOSP < 0.0 ) G_CLO_RELATIVE_BETA_HOSP = 0.0;
        }
        else if( str == "-scm" )
        {
            G_B_USESOCIALCONTACTMATRIX = true;
            //G_C_COMM[0][3] = 0.5;
        }
        else if( str == "-self-isolation-factor" )  {            G_CLO_SELFISOLATION_FACTOR = atof( argv[++i] );    }
        else if( str == "-steps-per-day" )          {            G_CLO_STEPS_PER_DAY = atof( argv[++i] );           }
        else if( str == "-susc-0-20" )
        {
            G_CLO_RELSUSC_0  = atof( argv[++i] );
            G_CLO_RELSUSC_10 = G_CLO_RELSUSC_0;
        }
        else if( str == "-susc-60-100" )
        {
            G_CLO_RELSUSC_60 = atof( argv[++i] );
            G_CLO_RELSUSC_70 = G_CLO_RELSUSC_60;
            G_CLO_RELSUSC_80 = G_CLO_RELSUSC_60;            
        }
        else if( str == "-symp-frac" ) // what you're really setting here is the symp fraction for 30-39 year-olds .. should be DEPRECATED
        {
            //G_CLO_SYMP_FRAC = atof( argv[++i] );
            //if( G_CLO_SYMP_FRAC > 0.325 ) G_CLO_SYMP_FRAC = 0.325;
            //if( G_CLO_SYMP_FRAC < 0.0 ) G_CLO_SYMP_FRAC = 0.0;
        }
        else if( str == "-symp-frac-davies" ) // use the Davies et al (Nature Medicine) numbers
        {
            G_B_SYMP_FRAC_DAVIES_20200616 = true;
            G_B_SYMP_FRAC_SIMPLEAVERAGE = false;
        }
        else if( str == "-symp-frac-equal" ) // here set the symp fraction to be equal for all age groups
        {
            G_B_SYMP_FRAC_DAVIES_20200616 = false;
            G_B_SYMP_FRAC_SIMPLEAVERAGE = false;
            G_CLO_SYMP_FRAC_EQUAL = atof( argv[++i] ); 
        }
        else if( str == "-tf" )                     {            G_CLO_TF = atof( argv[++i] );                      }
        else if( str == "-time-symp-to-hosp" )      {            G_CLO_TIME_SYMPTOHOSP = atof( argv[++i] );         }
        // ### FINAL ### IF BLOCK FOR AN UNKNOWN COMMAND-LINE OPTIONS            
        else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
        
            

        
        
        
        
        //END OF MAIN WHILE-LOOP BLOCK; INCREMENT AND MOVE ON
        i++;
    }
    
    // after parsing all of the arguments, do the following
    
    // ### 1 ### ERROR CHECK
    if( num_postld_contact_rates_set==0 || num_postld_contact_rates_set==8 )
    {
        // do nothing
    }
    else
    {
        fprintf(stderr, "\n\tIf you are setting the -contact-rate-postld-nn parameters, you need to set all eight of them (%d were set).\n\n", num_postld_contact_rates_set);
        exit(-1);
    }
    
    // ### 2 ### set the postld-contact data if it needs to be set
    /* if( !G_CLO_POSTLD_MIXING_SET )
    {
        // if the post-ld mixing paramters were not set, just set them to the regular (pre-ld) mixing params   
        for(int ac=1; ac<NUMAC; ac++)
        {
            G_CLO_MIXINGLEVEL_POSTLD[ac] = G_CLO_MIXINGLEVEL[ac];
        }
    } */
    
    // ### 3 ### set the specific age structure, population, and introduction time for this location
    SetLocationData( G_CLO_LOCATION );


    return;
}

void SetLocationData( string loc )
{
    for(int i=0;i<STARTK+NUMAC;i++) yic[i]=0.0; // zero everything out
    
    if( loc=="RI" )
    {
        // 2019 est pop is 1,059,361
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 1059361.0;
        yic[0]  = ppc->v[i_N]   *   .105;   //  0-9
        yic[1]  = ppc->v[i_N]   *   .123;   //  10-19
        yic[2]  = ppc->v[i_N]   *   .140;   //  20-29
        yic[3]  = ppc->v[i_N]   *   .127;   //  30-39
        yic[4]  = ppc->v[i_N]   *   .124;   //  40-49
        yic[5]  = ppc->v[i_N]   *   .135;   //  50-59
        yic[6]  = ppc->v[i_N]   *   .120;   //  60-69
        yic[7]  = ppc->v[i_N]   *   .074;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 55.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;
        //G_CLO_EARLYMARCH_ENDDAY = 84.5; 
        
    }
    else if( loc=="MA" )
    {
        // 2018 est pop is 6,897,212
        

        // data from https://www.census.gov/data/tables/time-series/demo/popest/2010s-national-total.html
        ppc->v[i_N] = 6897212.0;
        yic[0]  = ppc->v[i_N]   *   0.10586466;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.12243686;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.14498786;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.13384234;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.12230812;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.14064248;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.11801015;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.06958116;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 55.0; // 61.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;        
        
    }
    else if( loc=="PA" )
    {
        // 2019 est pop is 12,800,721
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 12800721.0;
        yic[0]  = ppc->v[i_N]   *   0.11160395;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.12229803;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.13156525;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.12581869;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.11809624;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.13878546;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.1270166;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.07657303;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 62.0;
        
        G_CLO_INTRODUCTION_COUNT = 2;        
        
    }
    else if( loc=="FL" )
    {
        // 2019 est pop is 21.5M
        

        // numbers from https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-detail.html
        ppc->v[i_N] = 21477737.0;
        yic[0]  = ppc->v[i_N]   *   0.107371275;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.112138583;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.124684365;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.126472682;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.121251461;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.133090744;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.125931750;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.094975183;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 61.0;
        
        G_CLO_INTRODUCTION_COUNT = 2;        
        
    }
    else if( loc=="CT" )
    {
        // 2020 est pop is 3.565M
        
        // numbers from https://data.census.gov/cedsci/table?q=connecticut&tid=ACSST1Y2019.S0101&hidePreview=false
        ppc->v[i_N] = 3565287.0;
        // fraction: 0.104918341777254, 0.127581313930688, 0.129650151586675, 0.123066109404376, 0.12305713396986, 0.143385932184422, 0.127443036142672, 0.074562300314112, 0.046335680689942
        yic[0]  = 374064;   //  0-9
        yic[1]  = 454864;   //  10-19
        yic[2]  = 462240;   //  20-29
        yic[3]  = 438766;   //  30-39
        yic[4]  = 438734;   //  40-49
        yic[5]  = 511212;   //  50-59
        yic[6]  = 454371;   //  60-69
        yic[7]  = 265836;   //  70-79
        yic[8]  = 165200;   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 61.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;        
        
    }
    else if( loc=="USA" )
    {
        // 332,599,000
        ppc->v[i_N] = 332599000.0;

        yic[0]  = 39888000;   //  0-9
        yic[1]  = 43455000;   //  10-19
        yic[2]  = 44501000;   //  20-29
        yic[3]  = 45194000;   //  30-39
        yic[4]  = 41153000;   //  40-49
        yic[5]  = 43033000;   //  50-59
        yic[6]  = 39034000;   //  60-69
        yic[7]  = 24136000;   //  70-79
        yic[8]  = 12205000;   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 55.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;        
        
    }
    else
    {
        fprintf(stderr, "\n\tUnknown location [%s] entered after -loc on command line.\n\n", loc.c_str() );
        exit(-1);
    }

    for (int ac = 0; ac < NUMAC; ac++){ ppc->v_pop_ac[ac] = yic[ac];    }
    
    
}





