/* This is the code for the EAS single fiber model described in:
 *  
 *   Kipping, D., and Nogueira, W. (under review). "A computational model of 
 *   a single auditory nerve fiber for electric-acoustic stimulation", to appear 
 *   in JARO - Journal of the Association for Research in Otolaryngology
 * 
 * Please cite this paper if you publish any research
 * results obtained with this code or any modified versions of this code.
 * 
 *
 * This implementation of the EAS model is based on code of I. Bruce et al. 
 * and S. Joshi et al. described in:
 * 
 *   Bruce, I.C., Erfani, Y., Zilany, M.S.A. (2018). "A phenomenological 
 *   model of the synapse between the inner hair cell and auditory nerve:
 *   Implications of limited neurotransmitter release sites." Hear. Res. 360, 
 *   40–54. https://doi.org/10.1016/j.heares.2017.12.016
 *
 *   Joshi, S.N., Dau, T., Epp, B. (2017). "A Model of Electrically Stimulated 
 *   Auditory Nerve Fiber Responses with Peripheral and Central Sites of Spike 
 *   Generation." JARO - J. Assoc. Res. Otolaryngol. 18, 323–342. 
 *   https://doi.org/10.1007/s10162-016-0608-2
 *
 *
 * Daniel Kipping, Aug 2021
 *
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 
 * Input:
 *   0: pAc   - sampled at  100 kHz
 *   1: pEl   - sampled at 1000 kHz (1 MHz)      
 *   2: nrep
 *   3: reptime
 *   4: cf
 *   5: cohc
 *   6: cihc
 *   7: noiseType
 *   8: implnt
 *   9: spont
 *  10: tabs
 *  11: trel
 *  12: Cper
 *  13: Ccen
 *  14: sigmaPer
 *  15: sigmaCen
 *  16: coupling
 * 
 * Output:
 *   0: psth_AS
 *   1: psth_ES
 *   2: meanrate
 *   3: varrate
 *   4: Sout
 *   5: trd_vector
 *   6: trel_vector
 *   7: vihc
 *   8: vmem
 *   9: isub
 *  10: isupra
 *  11: inoise
 *  12: ives
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <time.h>

#include "complex.hpp"

#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif

/* This function is the MEX "wrapper", to pass the input and output variables between the .mex* file and Matlab or Octave */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    const double dtAc = 1e-5 ;                  // time resolution of acoustic input in s (BEZ2018)    
    const double dtEl = 1e-6 ;                  // time resolution of electric input / timestep of Joshi2017 model in s (Joshi2017)
    /* note: it is not recommended to run the Joshi model with a time step larger than 1 us */

    void check_input_arguments(double*,int,double*,int,double*,int,double,double,double,double,double*,int,double,double,double,double,double,double,double,double,double,double,double,double*,int) ;

    /* --- Check for proper number of input arguments --- */
    if ((nrhs < 17) || (nrhs > 17)) {
        mexPrintf("nrhs = %d \n",nrhs); 
        mexErrMsgTxt("model_EAS requires 17 input arguments.");
    }
	if ((nlhs <  2) || (nlhs > 13)) {
        mexPrintf("nlhs = %d \n",nlhs); 
        mexErrMsgTxt("model_EAS requires 2 to 13 output argument.");
    }
    
    /* --- process inputs --- */
    
    double *pAcTmp        =  mxGetPr(prhs[ 0]);         // input sound wave in Pa (BEZ2018)
    double *pElTmp        =  mxGetPr(prhs[ 1]);         // input electric pulse train in A (Joshi2017)
    double *nrepTmp       =  mxGetPr(prhs[ 2]);         // number of stimulus repetitions for psth 
    double  reptime       = *mxGetPr(prhs[ 3]);         // time between stimulus repetitions in s 
    double  cf            = *mxGetPr(prhs[ 4]);         // characteristic frequency in Hz (BEZ2018)
    double  cohc          = *mxGetPr(prhs[ 5]);         // OHC impairment (BEZ2018)
    double  cihc          = *mxGetPr(prhs[ 6]);         // IHC impairment (BEZ2018)
    double *speciesTmp    = mxCalloc(1, sizeof(double)); 
            speciesTmp[0] = 1.0  ;                      // fixed cat species (1=cat; 2,3=human) (BEZ2018)      
    double  noiseType     = *mxGetPr(prhs[ 7]);         // fixed or variable fGn: 0=fixed; 1=variable (BEZ2018)
    double  implnt	      = *mxGetPr(prhs[ 8]);         // implementation of the power-law functions: 0=approxiate; 1=actual (BEZ2018)
    double  spont	      = *mxGetPr(prhs[ 9]);         // spontaneous firing rate of the fiber in /s (BEZ2018)
    double  tabsAc        = *mxGetPr(prhs[10]);         // absolute refractory period in s (BEZ2018)
    double  trelAc        = *mxGetPr(prhs[11]);         // baseline relative refractory period in s (BEZ2018)
    double  Cper          = *mxGetPr(prhs[12]);         // membrane capacitance for peripheral neuron (Joshi2017)
    double  Ccen          = *mxGetPr(prhs[13]);         // membrane capacitance for central neuron (Joshi2017)
    double  sigmaCen      = *mxGetPr(prhs[15]);         // noise current amplitude for central neuron (Joshi2017)
    double  sigmaPer      = *mxGetPr(prhs[14]);         // noise current amplitude for peripheral neuron (Joshi2017)
    double *couplingTmp   =  mxGetPr(prhs[16]);         // models are coupled for coupling=1,2 and uncoupled for coupling=0
        
    const int pAcBins     = (int) mxGetN(prhs[0]);      // length of input sound wave (BEZ2018)
    const int pElBins     = (int) mxGetN(prhs[1]);      // length of input electric stimulus (Joshi2017)
    
    const int nrep        = (int) nrepTmp[0];
    const int species     = (int) speciesTmp[0];
    const int coupling    = (int) couplingTmp[0];       

    /* --- check individual input arguments --- */
    check_input_arguments(pAcTmp,pAcBins,pElTmp,pElBins,nrepTmp,nrep,reptime,cf,cohc,cihc,speciesTmp,species,noiseType,implnt,spont,tabsAc,trelAc,Cper,Ccen,sigmaPer,sigmaCen,dtAc,dtEl,couplingTmp,coupling) ;
    mxFree(speciesTmp);

    /* --- create arrays for the return arguments --- */
    int nbinsAc = (int) floor(reptime/dtAc+0.5);        // number of samples for total repetition time, acoustic model
    int nbinsEl = (int) floor(reptime/dtEl+0.5);        // number of samples for total repetition time, electric model
    mwSize outsize[2];                                  // size of output arrays
    
    /* PSTH: */
    outsize[0] = 1; outsize[1] = nbinsAc;
    plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // psth_AS
    plhs[1] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // psth_ES

    /* acoustic model: */
    outsize[0] = 1; outsize[1] = nbinsAc;
    plhs[2] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // meanrate
    plhs[3] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // varrate
    
    outsize[1] = nbinsAc*nrep;
    plhs[4] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // Sout
    plhs[5] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // trd_vector
    plhs[6] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);  // trel_vector
    plhs[7] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);	 // vihc

    /* electric model: */
    outsize[0] = 2;
    plhs[ 8] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL); // vMem
    plhs[ 9] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL); // iSub
    plhs[10] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL); // iSupra
    plhs[11] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL); // iNoise
    plhs[12] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL); // iVes

    /* ---  assign pointers to the outputs --- */
    double *psth_AS	    = mxGetPr(plhs[ 0]);     // psth_AS
    double *psth_ES   	= mxGetPr(plhs[ 1]);     // psth_ES

    double *meanrate    = mxGetPr(plhs[ 2]);     // analytical estimate of the mean firing rate in /s for each time bin (BEZ2018)
    double *varrate     = mxGetPr(plhs[ 3]);     // analytical estimate of the variance in firing rate in /s for each time bin (BEZ2018)
    double *Sout        = mxGetPr(plhs[ 4]);     // synapse output rate in /s for each time bin (before the effects of redocking are considered) (BEZ2018)
    double *trd_vector  = mxGetPr(plhs[ 5]);     // vector of the mean redocking time in s for each time bin (BEZ2018)
    double *trel_vector = mxGetPr(plhs[ 6]);     // vector of the mean relative refractory period in s for each time bin (BEZ2018)
    double *vihc        = mxGetPr(plhs[ 7]);     // inner hair cell relative transmembrane potential in V (BEZ2018)
    
    double *vMem        = mxGetPr(plhs[ 8]);     // relative transmembrane potential in V for each time bin (Joshi2017)
    double *iSub        = mxGetPr(plhs[ 9]);     // subthreshold adaptation current in A for each time bin (Joshi2017)
    double *iSupra      = mxGetPr(plhs[10]);     // suprathreshold adaptation current in A for each time bin (Joshi2017)
    double *iNoise      = mxGetPr(plhs[11]);     // noise current in A for each time bin (Joshi2017)
    double *iVes        = mxGetPr(plhs[12]);     // vesicle release current in A for each time bin (only for coupling=2) (Joshi2017)
 
    /* ============================ */
    /* ======  run the model  ===== */
    /* ============================ */
    {   
        double *pAc = (double*) mxCalloc(nbinsAc, sizeof(double));  // mxCalloc initializes all elements to zero
        double *pEl = (double*) mxCalloc(nbinsEl, sizeof(double));  // mxCalloc initializes all elements to zero
        int     i;
        
        /* function declarations */
        void SingleFiber_EAS(double*,double,int,double*,double,int,int,double,double,int,double,double,double,double,double,double,double,double,double,double,int,     // input
                             double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);                                  // output
        
        /* Put stimuli (length=pAcBins or pElBins) into stimulus waveforms (length=nbinsAc or nbinsEl) */
        for (i=0; i<pAcBins; i++) pAc[i] = pAcTmp[i];
        for (i=0; i<pElBins; i++) pEl[i] = pElTmp[i];

        /* run Synapse and EAS spike generation model for IHC output and electric input */
        SingleFiber_EAS(pAc,dtAc,nbinsAc,pEl,dtEl,nbinsEl,nrep,noiseType,cf,species,cohc,cihc,implnt,spont,tabsAc,trelAc,Cper,Ccen,sigmaPer,sigmaCen,coupling,          // input
                        psth_AS,psth_ES,vihc,Sout,meanrate,varrate,trd_vector,trel_vector,vMem,iSub,iSupra,iNoise,iVes);                                                // output

        mxFree(pAc); mxFree(pEl);
    }

    // mexPrintf("finished EAS model.\n");
}

/* =========================== */
/* ======  Ac+El spiking ===== */
/* =========================== */

void SingleFiber_EAS(double *pAc, double dtAc, int nbinsAc, double *pEl, double dtEl, int nbinsEl, int nrep, double noiseType,                                          // input Ac+El
                     double cf, int species, double cohc, double cihc, double implnt, double spont, double tabsAc, double trelAc,                                       // parameters Ac
                     double Cper, double Ccen, double sigmaPer, double sigmaCen, int coupling,                                                                          // parameters El
                     double *psth_AS, double *psth_ES,                                                                                                                  // output psth
                     double *vihc, double *Sout, double *meanrate, double *varrate, double *trd_vector, double *trel_vector,                                            // output Ac
                     double *vMem, double *iSub, double *iSupra, double *iNoise, double *iVes) {                                                                        // output El

    /* ---------------------------------------------------------------- */
    /* ----- variables and preparation for BEZ2018 acoustic model ----- */
    /* ---------------------------------------------------------------- */

    /*variables for the signal-path, control-path and onward */
    double       *spTime;
    double       *spTimeAc;
    double       *spTimeEl;
    double        MeanISI;
    const double  SignalLength = nbinsAc * nrep * dtAc;
    long          maxSpikes;
    long          maxSpikesAc;
    long          maxSpikesEl;
    
    /* ---  Synaptic Release/Spike Generation Parameters --- */
    const double  tau       = 60.0e-3;                              // time constant for short-term adaptation (in mean redocking time)
    const double  t_rd_rest = 14.0e-3;                              // resting value of the mean redocking time
    const double  t_rd_jump =  0.4e-3;                              // size of jump in mean redocking time when a redocking event occurs
    const double  t_rd_init = t_rd_rest + 0.02e-3*spont-t_rd_jump;  // Initial value of the mean redocking time
    double        trelAc_i;
    const int     nSites = 4;                                       // number of synpatic release sites
    int           i;
    int           nSpikesAc = 0;
    int           nSpikesEl = 0;
    int          *nSpikesAcPtr = &nSpikesAc;
    int          *nSpikesElPtr = &nSpikesEl;
    int           ipst;
    double        I;
    const double  sampFreqSyn = 10e3;                               // Sampling frequency used in the synapse
    double        total_mean_rate = 0.0;
    
    /* --- function declarations --- */
    void   IHCAN(double*,double,int,double,int,double,double,int,double*);
    double Synapse(double*,double,double,int,int,double,double,double,double,double*);
    void   SpikeGenerator_EAS(double*,double,int,double*,double,int,int,double,double,double,double,double,    // input 
                                int,double,double,double,double,double,double,double,double,long,long,int,     // input
                                double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);    // output
        
    /* --- Run the IHC model --- */    
    // mexPrintf(" * running IHC model \n"); mexEvalString("pause(.001);"); 
    IHCAN(pAc,cf,nrep,dtAc,nbinsAc,cohc,cihc,species,vihc);

    /* --- Run the synapse model --- */
    // mexPrintf(" * running Synapse model \n"); mexEvalString("pause(.001);"); 
    I = Synapse(vihc,dtAc,cf,nbinsAc,nrep,spont,noiseType,implnt,sampFreqSyn,Sout);
    
    /* Calculate the overall mean synaptic rate */
    for(i = 0; i<I ; i++) {
        total_mean_rate = total_mean_rate + Sout[i]/I;
    };
    
    /* We register only the spikes at times after zero, the sufficient array size (more than 99.7% of cases) to register spike times  after zero is :/
     * /*MaxN=signalLengthInSec/meanISI+ 3*sqrt(signalLengthInSec/MeanISI)= nSpike_average +3*sqrt(nSpike_average)*/
    MeanISI = 1/total_mean_rate + t_rd_init/nSites + tabsAc + trelAc;
    // maxSpikesAc = ceil ( (long) ( SignalLength/MeanISI + 3* sqrt(SignalLength/MeanISI) ) ); // original value
    maxSpikesAc = ceil ( (long) 2 *  ( SignalLength/MeanISI + 3* sqrt(SignalLength/MeanISI) ) ) + 10 ; // amount ~doubled
    maxSpikesEl = ceil ( (long) 2 * ( ( SignalLength+10e-3)/tabsAc ) );                                                                      // PROVISIONAL - FIND A BETTER ESTIMATION !!!!!!!!
    // maxSpikes = ceil ( (long) ( SignalLength/MeanISI + 3* sqrt(SignalLength/MeanISI) ) + SignalLength/100e-6);               // PROVISIONAL - FIND A BETTER ESTIMATION !!!!!!!!

    spTimeAc = (double*) mxCalloc(maxSpikesAc, sizeof(double));
    spTimeEl = (double*) mxCalloc(maxSpikesEl, sizeof(double));

    do {
        if (nSpikesAc<0) { // Deal with cases where the acoustic spike time array was not long enough
            mxFree(spTimeAc);
            maxSpikesAc += 100 ; // Make the spike time array 100 elements larger
            spTimeAc = (double*) mxCalloc(maxSpikesAc, sizeof(double));
        }
        if (nSpikesEl<0) { // Deal with cases where the spike time array was not long enough
            mxFree(spTimeEl);
            maxSpikesEl += 100 ; // Make the spike time array 100 elements larger
            spTimeEl = (double*) mxCalloc(maxSpikesEl, sizeof(double));
        }
        SpikeGenerator_EAS(Sout,dtAc,nbinsAc,pEl,dtEl,nbinsEl,nrep,noiseType,t_rd_rest,t_rd_init,tau,t_rd_jump,nSites,tabsAc,trelAc,
                                    Cper,Ccen,sigmaPer,sigmaCen,spont,total_mean_rate,maxSpikesAc,maxSpikesEl,coupling,
                                    spTimeAc,spTimeEl,trd_vector,vMem,iSub,iSupra,iNoise,iVes,nSpikesAcPtr,nSpikesElPtr);
        
    } while ( nSpikesAc<0 || nSpikesEl<0 );  // Repeat if spike time array was not long enough

    /* Calculate the analytical estimates of meanrate and varrate and wrapping them up based on no. of repetitions */
    for(i = 0; i<I ; i++) {
        
        ipst = (int) (fmod(i, nbinsAc));
        if ( Sout[i]>0 ) {
            trelAc_i = __min(trelAc*100/Sout[i], trelAc);
            trel_vector[i] = trelAc_i;
            
            meanrate[ipst] = meanrate[ipst] + Sout[i]/(Sout[i]*(tabsAc + trd_vector[i]/nSites + trelAc_i) + 1)/nrep;  // estimated instantaneous mean rate
            varrate[ipst]  = varrate[ipst] + ((11*pow(Sout[i],7)*pow(trd_vector[i],7))/2 + (3*pow(Sout[i],8)*pow(trd_vector[i],8))/16 + 12288*pow(Sout[i],2)*pow(trelAc_i,2) + trd_vector[i]*(22528*pow(Sout[i],3)*pow(trelAc_i,2) + 22528*Sout[i]) + pow(trd_vector[i],6)*(3*pow(Sout[i],8)*pow(trelAc_i,2) + 82*pow(Sout[i],6)) + pow(trd_vector[i],5)*(88*pow(Sout[i],7)*pow(trelAc_i,2) + 664*pow(Sout[i],5)) + pow(trd_vector[i],4)*(976*pow(Sout[i],6)*pow(trelAc_i,2) + 3392*pow(Sout[i],4)) + pow(trd_vector[i],3)*(5376*pow(Sout[i],5)*pow(trelAc_i,2) + 10624*pow(Sout[i],3)) + pow(trd_vector[i],2)*(15616*pow(Sout[i],4)*pow(trelAc_i,2) + 20992*pow(Sout[i],2)) + 12288)/(pow(Sout[i],2)*pow(Sout[i]*trd_vector[i] + 4,4)*(3*pow(Sout[i],2)*pow(trd_vector[i],2) + 40*Sout[i]*trd_vector[i] + 48)*pow(trd_vector[i]/4 + tabsAc + trelAc_i + 1/Sout[i],3))/nrep; // estimated instananeous variance in the discharge rate
        }
        else trel_vector[i] = trelAc;
    };

    /* Generate PSTHs */
    for(i = 0; i < nSpikesAc; i++) { ipst = (int) (fmod(spTimeAc[i], dtAc*nbinsAc) / dtAc); psth_AS[ipst]  = psth_AS[ipst]  + 1; };
    for(i = 0; i < nSpikesEl; i++) { ipst = (int) (fmod(spTimeEl[i], dtAc*nbinsAc) / dtAc); psth_ES[ipst]  = psth_ES[ipst]  + 1; };

} // End of the SingleAN function

void SpikeGenerator_EAS(double *Sout, double dtAc, int nbinsAc, double *pEl, double dtEl, int nbinsEl, int nrep, double noiseType, 
                       double t_rd_rest, double t_rd_init, double tau, double t_rd_jump, int nSites, 
                       double tabsAc, double trelAc, double Cper, double Ccen, double sigmaPer, double sigmaCen, 
                       double spont, double total_mean_rate, long maxSpikesAc, long maxSpikesEl, int coupling,
                       double *spTimeAc, double *spTimeEl, double *trd_vector, 
                       double *vMem, double *iSub, double *iSupra, double *iNoise, double *iVes, int *nSpikesAcPtr, int *nSpikesElPtr) {

    /* ---------------------------------------------------------------- */
    /* ----- variables and preparation for BEZ2018 acoustic model ----- */
    /* ---------------------------------------------------------------- */

    /* ----- variables ----- */
    int     siteNo;                                                                                 // selected release site
    long    k, kInit;                                                                               // current time bin
    double  currRefracPeriod;                                                                       // end time of the current refractory period
    long    spCountAc  = 0;                                                                         // total number of acoustically evoked spikes fired
    long    spCountEl  = 0;                                                                         // total number of electrically evoked spikes fired
    int     t_rd_decay = 1;                                                                         // Logical "true" as to whether to decay the value of currRedockPeriod at the end of the time step 
    int     rd_first   = 0;                                                                         // Logical "false" as to whether to a first redocking event has occurred
    double  prevRedockPeriod;                                                                       // dynamic mean redocking time
    double  currRedockPeriod;                                                                       // dynamic mean redocking time

    double *randNums;                                                                               // random number array
    long    randBufIndex = 0;                                                                       // index indicating the next random number to be used
    long    randBufLen;

    double *oneSiteRedock_ptr    = (double*) mxCalloc(nSites, sizeof(double));
    double *prevReleaseTimes_ptr = (double*) mxCalloc(nSites, sizeof(double));    
    int    *unitRateInterval_ptr = (int*)    mxCalloc(nSites, sizeof(double));
    double *elapsedTime_ptr      = (double*) mxCalloc(nSites, sizeof(double));
    double *currReleaseTimes_ptr = (double*) mxCalloc(nSites, sizeof(double));
    double *Xsum_ptr             = (double*) mxCalloc(nSites, sizeof(double));
    double *preReleaseTimeBinsSorted;   
    double *currRefracPeriod_ptr = &currRefracPeriod;
    int    *t_rd_decay_ptr       = &t_rd_decay;
    int    *rd_first_ptr         = &rd_first;
    double *prevRedockPeriod_ptr = &prevRedockPeriod;
    double *currRedockPeriod_ptr = &currRedockPeriod;
    long   *randBufIndex_ptr     = &randBufIndex;

    int     registerSpike(double*, int, double) ;
    
    /* ----- initialize variables ----- */
    {
        /* --- create random number array: --- */
        /* Estimating Max number of spikes and events (including before zero elements)
        * The sufficient array size (more than 99.7% of cases) to register event times  after zero is:
        * MaxN=signalLengthInSec/meanEvents+ 3*sqrt(signalLengthInSec/MeanEvents) */
        double    MeanInterEvents = (1/total_mean_rate) + (t_rd_init)/nSites;
        long      maxEventsAc     = ceil( (long) (nbinsAc*nrep*dtAc/MeanInterEvents + 3*sqrt(nbinsAc*nrep*dtAc/MeanInterEvents)) ) + nSites;

        /* Max random array Size:   
        *      - nSites elements for oneSiteRedock initialization 
        *      - nSites elements for preRelease_initialGuessTimeBins initialization
        *      - 1 element for Tref initialization
        *      - maxSpikesAc elements for Tref in the loop
        *      - maxEventsAc elements one time for redocking another time for rate intervals
        *      - averageley 2nSites for before zero elements (redock-and unitRate)
        *      - nSites (Max) for Trefs */
        randBufLen  = (long) ceil( 2*nSites + 1 + 2*maxSpikesAc + 2*maxEventsAc + 3*nSites);                // original value for acoustic model only
        randBufLen += 2*maxSpikesEl ;
        
        /* call to 'myrand' */
        mxArray *randInputArray[2];                                                                         // input array for the MATLAB random number generator
        mxArray *randOutputArray[1];                                                                        // output array for the MATLAB random number generator
        
        randInputArray[0] = mxCreateDoubleMatrix(1, 2, mxREAL);                                             // size of the random number array
        randInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);                                             // fixed or variable random numbers
        (mxGetPr(randInputArray[0]))[0] = 1;
        (mxGetPr(randInputArray[0]))[1] = randBufLen;
        (mxGetPr(randInputArray[1]))[0] = noiseType;
        mexCallMATLAB(1, randOutputArray, 2, randInputArray, "Model.myrand");                                     // call MATLAB rand

        randNums = (double*) mxCalloc(randBufLen, sizeof(double));
        for (k=0; k<randBufLen; k++) {
            randNums[k] = (mxGetPr(randOutputArray[0]))[k];
        }
        mxDestroyArray(randInputArray[0]); 
        mxDestroyArray(randOutputArray[0]);    

    
        /* --- set initial values --- */
        double preRelease_initialGuessTimeBins[nSites];

        for (siteNo=0; siteNo<nSites; siteNo++) {

            double randNum = randNums[randBufIndex++] ;

            preRelease_initialGuessTimeBins[siteNo] = __max( - nbinsAc*nrep,                                // Initial preRelease_initialGuessTimeBins associated to nSites release sites
                ceil( (nSites/__max(Sout[0],0.1) + t_rd_init)*log(randNum)/dtAc ) );
            
        }

        /* Now Sort the initial preRelease times and associate the farthest to zero as the site which has also generated a spike */
        mxArray *sortInputArray[1];                                                                         // input array for the MATLAB sort function
        mxArray *sortOutputArray[1];                                                                        // output array for the MATLAB sort function

        sortInputArray[0] = mxCreateDoubleMatrix(1, nSites, mxREAL);                            
        for (siteNo=0; siteNo<nSites; siteNo++) {
            (mxGetPr(sortInputArray[0]))[siteNo] = preRelease_initialGuessTimeBins[siteNo];
        }
        mexCallMATLAB(1, sortOutputArray, 1, sortInputArray, "sort");                                       // call MATLAB sort

        preReleaseTimeBinsSorted = (double*) mxCalloc(nSites, sizeof(double));
        for (siteNo=0; siteNo<nSites; siteNo++) {
            preReleaseTimeBinsSorted[siteNo] = (mxGetPr(sortOutputArray[0]))[siteNo];
        }

        mxDestroyArray(sortInputArray[0]); 
        mxDestroyArray(sortOutputArray[0]);

        for (siteNo=0; siteNo<nSites; siteNo++) {
            oneSiteRedock_ptr[siteNo]    = - t_rd_init*log(randNums[randBufIndex++]);                       // Initial redocking time associated to nSites release sites
            prevReleaseTimes_ptr[siteNo] = ( (double) preReleaseTimeBinsSorted[siteNo] )*dtAc;              // Consider the inital prevReleaseTimes to be the preReleaseTimeBinsSorted*dtAc
        }
        kInit = (int) preReleaseTimeBinsSorted[0];                                                          // The position of first spike, also where the process is started - continued from the past
        /* ============ */
        // To ensure that the stimulus is preceded by at least 10 ms of silence
        // This is important for the electric model in order to be in a randomized initial state at t=0
        kInit = __min( kInit, -ceil(10e-3/dtAc)) ;
        /* ============ */
        currRefracPeriod = (double) kInit*dtAc;
        prevRedockPeriod = t_rd_init;                                                                       // dynamic mean redocking time
        currRedockPeriod = prevRedockPeriod;                                                                // dynamic mean redocking time
    }
    
    /* ---------------------------------------------------------------- */
    /* ----- variables and preparation for JDE2017 electric model ----- */
    /* ---------------------------------------------------------------- */

    /* ----- parameters ----- */
    const double alpha     = 0.80;                                                                          // noise shaping parameter alpha
    const double sigma_per = sigmaPer ;
    const double sigma_cen = sigmaCen ;

    /* ----- model variables ----- */
    double  tLastSpike     = -10000 ;                                                                       // time of last spike occurrence in [s]
    int     isARP          =      0 ;                                                                       // logical 'false' as to whether the neurons are in an absolute refractory period
    int     neurNo;                                                                                         // selected neuron (0=peripheral neuron; 1=central neuron)
    long    i, iInit ;                                                                                      // counter

    double *tLastSpike_ptr = &tLastSpike;
    int    *isARP_ptr      = &isARP;

    /* for coupling=2 */
    double  tLastRelease      = -1000 ;
    double *tLastRelease_ptr  = &tLastRelease ;
    double  tTransmitter      = 40e-6 ;                     // duration of transmitter-released current
    double  iTransmitter      = 3e-3 ;                    
    double  iVesTmp = 0 ;

    /* allocate space for temporary copies (JDE2017 model) */
    double *vMemTmp    = (double*) mxCalloc(2, sizeof(double));
    double *iSubTmp    = (double*) mxCalloc(2, sizeof(double));
    double *iSupraTmp  = (double*) mxCalloc(2, sizeof(double));
    
    /* ----------------------------------------------------------- */
    /* ----- compute spike times for both models in parallel ----- */
    /* ----------------------------------------------------------- */

    int acousticSpikeThisStep    = 0;
    int electricSpikeThisStep    = 0;
    int electricSpikeThisStepTmp = 0;

    int    timestep_BEZ2018(int,long,double,double,double,double*,double,double,double,double*,double,int*,int*,double*,double*,double*,double*,double*,int*,double*,double*,double*,double*,long*);
    int    timestep_JDE2017(long,double,double,double,double,double,double,double,double,double,int*,double*,double*,double*,double*);
    double spikeAdaptation_BEZ2018(double, double, double, double, double);
    double spikeAdaptation_JDE2017(int*, double*, double*, double);

    k     = kInit;
    iInit = round(kInit*dtAc/dtEl) ;
    i     = iInit ;
    
    /* calculate noise currents for JDE2017 model: */
    void    calcShapedNoise(double*,int,double,double,double,double);
    double *iNoiseP_ptr = (double*) mxCalloc(10*(nbinsAc*nrep+labs(kInit)+1), sizeof(double));
    double *iNoiseC_ptr = (double*) mxCalloc(10*(nbinsAc*nrep+labs(kInit)+1), sizeof(double));
    calcShapedNoise(iNoiseP_ptr, 10*(nbinsAc*nrep+labs(kInit)+1), dtEl, alpha, noiseType, sigma_per) ;
    calcShapedNoise(iNoiseC_ptr, 10*(nbinsAc*nrep+labs(kInit)+1), dtEl, alpha, noiseType, sigma_cen) ;
    
    while (k < nbinsAc*nrep){ // a loop to find the spike times for all the nbinsAc*nrep ; time = k * dtAc

        /* ------------------------------------------------------- */
        /* ----- propagate acoustic system to next time step ----- */
        /* ------------------------------------------------------- */

        /* call BEZ2018 model */
        acousticSpikeThisStep = timestep_BEZ2018(nSites, k, dtAc, t_rd_jump, t_rd_rest, preReleaseTimeBinsSorted, Sout[__max(0,k)], trelAc, tabsAc, randNums, tau,
                                                 t_rd_decay_ptr, rd_first_ptr, oneSiteRedock_ptr, elapsedTime_ptr, currRedockPeriod_ptr, prevRedockPeriod_ptr, Xsum_ptr, 
                                                 unitRateInterval_ptr, tLastRelease_ptr, currReleaseTimes_ptr, prevReleaseTimes_ptr, currRefracPeriod_ptr, randBufIndex_ptr);

        /* ------------------------------------------------------- */
        /* ----- propagate electric system to next time step ----- */
        /* ------------------------------------------------------- */   
        
        electricSpikeThisStep = 0;
        // int itmp = 0 ;
        while ( (i+1)*dtEl <= (k+1)*dtAc+dtEl/10 ) { // because of dtEl<dtAc, the electric system needs dtAc/dtEl steps to reach the same time as the acoustic part
            
            if ((2==coupling) &&  ( i*dtEl-tLastRelease <= tTransmitter )) {
                iVesTmp = iTransmitter ;
            } else {
                iVesTmp = 0.0 ;
            }
            
            /* call JDE2017 model */
            if ( iInit==i ) {   // initialize variables at first call (flag dtEl < 0):
                electricSpikeThisStepTmp = timestep_JDE2017(0,-1.0,0.0,0.0,0.0,0.0,tabsAc,trelAc,Cper,Ccen,isARP_ptr,tLastSpike_ptr,vMemTmp,iSubTmp,iSupraTmp);
            }
            if ( i<0 ) {        // no stimulus for t<0
                electricSpikeThisStepTmp = timestep_JDE2017(i,dtEl,0.0      ,iNoiseP_ptr[i-iInit],iNoiseC_ptr[i-iInit],iVesTmp,tabsAc,trelAc,Cper,Ccen,isARP_ptr,tLastSpike_ptr,vMemTmp,iSubTmp,iSupraTmp);
            } else {            // simulate stimulation for t>=0
                int ipst = (int) (fmod(i, nbinsEl));
                electricSpikeThisStepTmp = timestep_JDE2017(i,dtEl,pEl[ipst],iNoiseP_ptr[i-iInit],iNoiseC_ptr[i-iInit],iVesTmp,tabsAc,trelAc,Cper,Ccen,isARP_ptr,tLastSpike_ptr,vMemTmp,iSubTmp,iSupraTmp);
            }
                
            if (1==electricSpikeThisStepTmp) electricSpikeThisStep = 1;
            i += 1 ;
        } 
        
        /* ----------------------------- */
        /* ----- handle EAS spikes ----- */
        /* ----------------------------- */
        
        if (2==coupling) { /* COUPLING 2.0 */
            /* models are coupled - save the spike times of the electric model (acoustic model is ignored, as vesicle releases are fed into the electric model): */

            if ( 1==electricSpikeThisStep ) { // electric model is spiking
                spCountEl = registerSpike(spTimeEl, spCountEl, k*dtAc) ;
            } 
        } 
        else if (1==coupling) { /* COUPLING 1.0 */
            /* models are coupled - each spike triggers the adaptation process of the other model */

            if ( ( 1==electricSpikeThisStep ) && (1==acousticSpikeThisStep) ) { 
                /* both models are spiking: select the electric/acoustic spike randomly */
                if ( randNums[randBufIndex]<=.5 ) { spCountEl = registerSpike(spTimeEl, spCountEl, k*dtAc) ; } // for electric spike
                else                              { spCountAc = registerSpike(spTimeAc, spCountAc, k*dtAc) ; } // for acoustic spike
            } else if (1==electricSpikeThisStep) { 
                /* only electric model is spiking: trigger spike adaptation in acoustic model */
                spCountEl        = registerSpike(spTimeEl, spCountEl, k*dtAc) ;
                double randNum   = randNums[randBufIndex++];
                currRefracPeriod = spikeAdaptation_BEZ2018(trelAc, tabsAc, Sout[__max(0,k)], randNum, tLastSpike);
            } else if (1==acousticSpikeThisStep) { 
                /* only acoustic model is spiking: trigger spike adaptation in electric model */
                spCountAc        = registerSpike(spTimeAc, spCountAc, k*dtAc) ;
                tLastSpike       = spikeAdaptation_JDE2017(isARP_ptr, iSupraTmp, vMemTmp, k*dtAc);
            }
        } 
        else {
            /* models are uncoupled - simply save the electric/acoustic spike times: */
            if ( 1==electricSpikeThisStep ) { spCountEl = registerSpike(spTimeEl, spCountEl, k*dtAc) ; } // electric model is spiking
            if ( 1==acousticSpikeThisStep ) { spCountAc = registerSpike(spTimeAc, spCountAc, k*dtAc) ; } // acoustic model is spiking
        }

        /* ---------------------------------------------- */   
        /* ----- save the results of this time step ----- */
        /* ---------------------------------------------- */   
        if (k>=0) {
            trd_vector[k] = currRedockPeriod;

            /* store values of the JDE2017 model: */
            vMem[2*k  ]   = vMemTmp[0];
            vMem[2*k+1]   = vMemTmp[1];
            iSub[2*k  ]   = iSubTmp[0];
            iSub[2*k+1]   = iSubTmp[1];
            iSupra[2*k  ] = iSupraTmp[0];
            iSupra[2*k+1] = iSupraTmp[1];
            iNoise[2*k  ] = iNoiseP_ptr[k-kInit];
            iNoise[2*k+1] = iNoiseC_ptr[k-kInit];
            iVes[2*k  ]   = iVesTmp ;
            iVes[2*k+1]   = 0 ; 
        }

        /* -------------------------- */   
        /* ----- Error Catching ----- */
        /* -------------------------- */   
        if ( (randBufIndex+3*nSites)>randBufLen ) {
            mexPrintf("@ timestep k = %d / %d (kInit = %d): too few random numbers left, rerunning the function. \n", k, nbinsAc*nrep, (int) kInit );
            mexEvalString("pause(.001);"); 
            spCountAc = -1;
            spCountEl = -1;
            break;
        }
        if ( (spCountAc+1)>maxSpikesAc ) {
            mexPrintf("@ timestep k = %d / %d (kInit = %d): too many acoustically evoked spikes (%d > %d), rerunning the function. \n", k, nbinsAc*nrep, (int) kInit, (int) spCountAc, (int) maxSpikesAc );
            mexEvalString("pause(.001);"); 
            spCountAc = -1;
            break;
        }
        if ( (spCountEl+1)>maxSpikesEl ) {
            mexPrintf("@ timestep k = %d / %d (kInit = %d): too many electrically evoked spikes (%d > %d), rerunning the function. \n", k, nbinsAc*nrep, (int) kInit, (int) spCountEl, (int) maxSpikesEl );
            mexEvalString("pause(.001);"); 
            spCountEl = -1;
            break;
        }

        /* ------------------------------------------------------- */   

        k = k+1;

    };

    *nSpikesAcPtr = spCountAc;
    *nSpikesElPtr = spCountEl;

    mxFree(oneSiteRedock_ptr);
    mxFree(elapsedTime_ptr);
    mxFree(Xsum_ptr);
    mxFree(unitRateInterval_ptr);
    mxFree(currReleaseTimes_ptr);
    mxFree(prevReleaseTimes_ptr);

    mxFree(vMemTmp);
    mxFree(iSubTmp);
    mxFree(iSupraTmp);
    mxFree(iNoiseP_ptr);
    mxFree(iNoiseC_ptr);

    mxFree(randNums);
    mxFree(preReleaseTimeBinsSorted);

}

int registerSpike(double *spikeTimes, int spikeCount, double time) {
    if ( time>0 ) {
        spikeTimes[spikeCount] = time ;
        return spikeCount+1 ;  
    } else {
        return spikeCount ; // register only positive spike times
    }
}

int timestep_BEZ2018(int nSites, long k, double dtAc, double t_rd_jump, double t_rd_rest, double *preReleaseTimeBinsSorted,             // in
                     double Sout_k, double trelAc, double tabsAc, double *randNums, double tau,                                         // in
                     int *t_rd_decay, int *rd_first, double *oneSiteRedock, double *elapsedTime, double *currRedockPeriod,              // in / out
                     double *prevRedockPeriod,  double *Xsum, int *unitRateInterval, double *tLastRelease, double *currReleaseTimes,    // in / out
                     double *prevReleaseTimes, double *currRefracPeriod, long *randBufIndex)                                            // in / out
{                     
    int    siteNo;                     // release site index
    double trelAc_k;                   // mean relative refractory time
    int    oneSiteRedock_rounded;
    int    elapsedTime_rounded;
    int    acousticSpikeThisStep = 0;
    double Tref;

    double spikeAdaptation_BEZ2018(double, double, double, double, double);

    for (siteNo=0; siteNo<nSites; siteNo++) { // loop for all nSites

        if ( k > preReleaseTimeBinsSorted[siteNo] ) { // to be sure that for each site, the code start from its associated previus release time
            
            /* redocking times do not necessarily occur exactly at time step value
                * - calculate the number of integer steps for the elapsed time and redocking time */
            oneSiteRedock_rounded = (int) floor(oneSiteRedock[siteNo]/dtAc);
            elapsedTime_rounded  = (int) floor( elapsedTime[siteNo]/dtAc);

            if ( oneSiteRedock_rounded == elapsedTime_rounded ) { // a redocking event occurred!
                (* currRedockPeriod) = (* prevRedockPeriod) + t_rd_jump;     // Jump  trd by t_rd_jump if a redocking event has occurred
                (* prevRedockPeriod) = (* currRedockPeriod);
                (* t_rd_decay) = 0;                                          // Don't decay the value of currRedockPeriod if a jump has occurred
                (* rd_first)   = 1;                                          // Flag for when a jump has first occurred
            }
            
            elapsedTime[siteNo] = elapsedTime[siteNo] + dtAc;
        };
        /* --------------------------------------------------------------------------- */

        /* the elapsed time passes the one time redock (the redocking is finished),
            * In this case the synaptic vesicle starts sensing the input
            * for each site integration starts after the redocking is finished for the corresponding site)*/
        if ( elapsedTime[siteNo] >= oneSiteRedock[siteNo] ) {
            Xsum[siteNo] = Xsum[siteNo] + Sout_k/nSites;                    // there are nSites integrals each vesicle senses 1/nosites of the whole rate
        }
        /* --------------------------------------------------------------------------- */

        if ( (Xsum[siteNo] >= unitRateInterval[siteNo]) &&  (k >= preReleaseTimeBinsSorted[siteNo]) ) { // An event - a release  happened for the siteNo

            oneSiteRedock[siteNo]    = - (* currRedockPeriod) * log(randNums[(* randBufIndex)++]);
            currReleaseTimes[siteNo] =   prevReleaseTimes[siteNo] + elapsedTime[siteNo];
            tLastRelease[0]          =   currReleaseTimes[siteNo] ;
            elapsedTime[siteNo] = 0;               
            
            if (currReleaseTimes[siteNo] >= (* currRefracPeriod)) {         // A spike occured for the current event - release

                // if (currReleaseTimes[siteNo] >= 0) {                        // Register only non negative spike times
                    // acousticSpikeThisStep = 1;                              // set acoustic spike flag to 1
                // }

                /* Register also negative spike times for electric-acoustic coupling: */
                acousticSpikeThisStep = 1;                                     // set acoustic spike flag to 1
                
                double randNum = randNums[(* randBufIndex)++];
                (* currRefracPeriod) = spikeAdaptation_BEZ2018(trelAc, tabsAc, Sout_k, randNum, currReleaseTimes[siteNo]);
            }
            
            prevReleaseTimes[siteNo] = currReleaseTimes[siteNo];
            Xsum[siteNo] = 0;
            unitRateInterval[siteNo] = (int) (-log(randNums[(* randBufIndex)++])/dtAc);
        };
        /* --------------------------------------------------------------------------- */

    }; // end of loop for nSites
    
    /* Decay the adapative mean redocking time towards the resting value if no redocking events occurred in this time step */
    if ( (1 == (* t_rd_decay)) && (1 == (* rd_first)) ) {
        (* currRedockPeriod) = (* prevRedockPeriod) - (dtAc/tau) * ( (* prevRedockPeriod) - t_rd_rest );
        (* prevRedockPeriod) = (* currRedockPeriod);
    } else {
        (* t_rd_decay) = 1;
    }
    
    return acousticSpikeThisStep;
}

double spikeAdaptation_BEZ2018(double trelAc, double tabsAc, double Sout_k, double randNum, double time_of_spike) {

    double trelAc_k = __min(trelAc*100/Sout_k, trelAc);
    double Tref     = tabsAc - trelAc_k*log(randNum);            // Refractory periods
    
    double currRefracPeriod = time_of_spike + Tref;

    return currRefracPeriod;
}

int timestep_JDE2017(long i, double dtEl, double pEl, double iNoisePer, double iNoiseCen, double iVesTmp, double tabsAc, double trelAc, double Cper, double Ccen,   // in
                     int *isARP, double *tLastSpike, double *vMem, double *iSub, double *iSupra)                                                                    // in/out
{   
    
    /* model parameters (for peripheral/central neuron): */

    const double       gL[2] = {      1.1e-3 ,      2.7e-3 } ;  // membrane conductance g_L in [S]
    const double        C[2] = {        Cper ,        Ccen } ;  // membrane capacitance C in [F]
    const double     DelT[2] = {       10e-3 ,        3e-3 } ;  // slope factor deltaT in [V]
    const double          EL =        -80e-3 ;                  // resting potential E_L in [V]
    const double        Vthr =        -70e-3 ;                  // threshold potential v_threshold in [V]
    const double       Vpeak =         24e-3 ;                  // peak potential v_peak in [V]
    const double   tauSub[2] = {      250e-6 ,      250e-6 } ;  // subthreshold adaptation time constant tau_sub in [s]
    double             trel0 =      512.5e-6 ;
    double       tauSupraPer = 4500e-6 * trelAc/trel0 ;
    double       tauSupraCen = 2500e-6 * trelAc/trel0 ;
    const double tauSupra[2] = { tauSupraPer , tauSupraCen } ;  // suprathreshold adaptation time constant tau_supra in [s]
    const double        aSub =         2.0e-3 ;                  // subthreshold adaptation conductance a_sub in [S]
    const double      aSupra =         3.0e-3 ;                  // suprathreshold adaptation conductance a_supra in [S]

    const double       tDead =        tabsAc ;                  // dead time in [s]
    const double       beta  =          0.75 ;                  // inhibitory compression beta    

    /* technical variables: */
    int    n;
    double hv;
    double iInput[2] ;
    double vMem_dt;
    double iSub_dt;
    double iSupra_dt;

    int    electricSpikeThisStep = 0;

    /* function declaration: */
    double spikeAdaptation_JDE2017(int*, double*, double*, double);

    if ( dtEl<0 ) { // initialize the variables at the first call
        vMem[0] = vMem[1] = EL;
        iSub[0] = iSub[1] = iSupra[0] = iSupra[1] = 0;
        return 0;
    }

    /* check if ARP ends : */
    if (i*dtEl > (* tLastSpike) + tDead) {
        (* isARP) = 0 ;
    }

    /* compute stimulus input : */
    if (1 == (*isARP)) {
        iInput[0] = iNoisePer;
        iInput[1] = iNoiseCen;
    } 
    else if (pEl >= 0) {
        iInput[0] = iNoisePer - beta * pEl + iVesTmp;
        iInput[1] = iNoiseCen +        pEl ;
    } else {
        iInput[0] = iNoisePer -        pEl + iVesTmp ;
        iInput[1] = iNoiseCen + beta * pEl ;
    }

    /* compute values for next timestep : */
    for (n=0; n<=1; n++) { // n=0: peripheral neuron / n=1: central neuron

        // /* calculate values for next timestep using forward Euler : */
        hv  = - gL[n]*(vMem[n] - EL) + gL[n]*DelT[n]*exp( (vMem[n] - Vthr) / DelT[n] ) ;
        vMem_dt   = hv - iSub[n] - iSupra[n] + iInput[n];
        iSub_dt   = aSub   * (vMem[n] - EL) - iSub[n];
        iSupra_dt = aSupra * (vMem[n] - EL) - iSupra[n];

        vMem[n]   += dtEl * vMem_dt / C[n] ;
        iSub[n]   += dtEl * iSub_dt / tauSub[n] ;
        iSupra[n] += dtEl * iSupra_dt / tauSupra[n] ;

        /* detect spiking : */
        if ( vMem[n] >= Vpeak ) {
            if (*isARP) {
                vMem[n] = Vpeak ;                                       // ignore spiking during ARP (catch errors)
            } else {
                electricSpikeThisStep = 1;                              // set electric spike flag to 1    
            }
        }
    }

    /* if a spike occurred in this timestep, call spike handling routine */
    if ( 1 == electricSpikeThisStep ) {
        // update last spike time
        (* tLastSpike) = spikeAdaptation_JDE2017(isARP, iSupra, vMem, i*dtEl);
    }

    return electricSpikeThisStep;
}

double spikeAdaptation_JDE2017(int *isARP, double *iSupra, double *vMem, double time_of_spike) {

    const double Vres = -84e-3;                           // reset potential v_reset in [V]
    /* ---------------------------------------------- */
    const double b    =  90e-6;                           // spike-triggered offset for iSupra in [A]
    /* ---------------------------------------------- */
    double tLastSpike = time_of_spike;

    (* isARP) = 1 ;                                       // set neuron into absolute refractory period

    iSupra[0] += b ;                                      // shift suprathreshold adaptation currents by b
    iSupra[1] += b ;   

    vMem[0] = Vres ;                                      // reset both neurons to reset potential Vres
    vMem[1] = Vres ;                  

    return tLastSpike;

}

/* ============================= */
/* ======  error checking  ===== */
/* ============================= */

void check_input_arguments(double *pAcTmp, int pAcBins, double *pElTmp, int pElBins, double *nrepTmp, int nrep, double reptime, 
    double cf, double cohc, double cihc, double *speciesTmp, int species, double noiseType, double implnt, double spont, double tabsAc, double trelAc, 
    double Cper, double Ccen, double sigmaPer, double sigmaCen, double dtAc, double dtEl, double *couplingTmp, int coupling) {
    
    /* general */
    if (pAcBins==1) mexErrMsgTxt("pAc must be a row vector\n");
    if (pElBins==1) mexErrMsgTxt("pEl must be a row vector\n");
    if (nrepTmp[0]!=nrep) { mexPrintf("nrep = %f \n",nrepTmp[0]); mexErrMsgTxt("nrep must an integer.\n"); }
    if (nrep<1) { mexPrintf("nrep = %f \n",nrepTmp[0]); mexErrMsgTxt("nrep must be greater than 0.\n"); }
    if ( (reptime<(pAcBins-0.1)*dtAc) || (reptime<(pElBins-0.1)*dtEl) ) {
        mexPrintf("reptime = %f s \n",reptime) ;
        mexPrintf("stimulus duration (ac) = %f s \n",pAcBins*dtAc) ;
        mexPrintf("stimulus duration (el) = %f s \n",pElBins*dtEl) ;
        mexErrMsgTxt("reptime should be equal to or longer than the stimulus duration.\n");
    }
    if ( fabs(pAcBins*dtAc - pElBins*dtEl) > __min(dtAc,dtEl) ) {
        mexPrintf("pAc and pEl must have the same duration.\n");
        mexPrintf("pAc must be sampled at %.0f kHz and pEl must be sampled at %.0f kHz.\n",1e-3/dtAc,1e-3/dtEl); 
        mexErrMsgTxt("\n");
    }
    if (couplingTmp[0]!=coupling) { mexPrintf("coupling = %f \n",couplingTmp[0]); mexErrMsgTxt("coupling must be an integer.\n"); }
    if ( (coupling<0) || (coupling>2) ) { mexPrintf("coupling = %d \n",coupling); mexErrMsgTxt("coupling must be 0, 1 or 2\n"); }

    /* for IHC arguments */
    if ((cohc<0) || (cohc>1)) { mexPrintf("cohc (= %1.1f) must be between 0 and 1\n",cohc); mexErrMsgTxt("\n"); }
    if ((cihc<0) || (cihc>1)) { mexPrintf("cihc (= %1.1f) must be between 0 and 1\n",cihc); mexErrMsgTxt("\n"); }
    if (speciesTmp[0]!=species) mexErrMsgTxt("species must an integer.\n");
    if (species<1 || species>3) mexErrMsgTxt("Species must be 1 for cat, or 2 or 3 for human.\n");
    if ((species==1) && ( (cf<124.9) || (cf>40.1e3) )) { mexPrintf("cf (= %1.1f Hz) must be between 125 Hz and 40 kHz for cat model\n",cf); mexErrMsgTxt("\n"); }
    if ((species>1)  && ( (cf<124.9) || (cf>20.1e3) )) { mexPrintf("cf (= %1.1f Hz) must be between 125 Hz and 20 kHz for human model\n",cf); mexErrMsgTxt("\n"); }
    
    /* for Synapse arguments */
    if ((spont<1e-4)||(spont>180))  { mexPrintf("spont = %f \n" ,spont ); mexErrMsgTxt("spont must in the range [1e-4,180]\n"); }
    if ((tabsAc<0)||(tabsAc>20e-3)) { mexPrintf("tabsAc = %f \n",tabsAc); mexErrMsgTxt("tabsAc must in the range [0,20e-3]\n"); }   
    if ((trelAc<0)||(trelAc>20e-3)) { mexPrintf("trelAc = %f \n",trelAc); mexErrMsgTxt("trelAc must in the range [0,20e-3]\n"); }   

    /* for membrane capacitances */
    if ((Cper<=0)||(Cper>1e-3)) { mexPrintf("Cper = %f \n",Cper); mexErrMsgTxt("Cper must in the range (0,1e-3]\n"); }   
    if ((Ccen<=0)||(Ccen>1e-3)) { mexPrintf("Ccen = %f \n",Ccen); mexErrMsgTxt("Ccen must in the range (0,1e-3]\n"); }   

    /* for noise current amplitudes */
    if ((sigmaPer<=0)||(sigmaPer>1e-3)) { mexPrintf("sigmaPer = %f \n",Cper); mexErrMsgTxt("sigmaPer must in the range (0,1e-3]\n"); }   
    if ((sigmaCen<=0)||(sigmaCen>1e-3)) { mexPrintf("sigmaCen = %f \n",Ccen); mexErrMsgTxt("sigmaCen must in the range (0,1e-3]\n"); }   
}

/* ============================ */
/* ======  IHC functions  ===== */
/* ============================ */

void IHCAN(double *pAc, double cf, int nrep, double dtAc, int nbinsAc, double cohc, double cihc, int species, double *ihcout) {	
    
    /*variables for middle-ear model */
	double megainmax;
    double *mey1, *mey2, *mey3, meout,c1filterouttmp,c2filterouttmp,c1vihctmp,c2vihctmp;
    double fp,C,m11,m12,m13,m14,m15,m16,m21,m22,m23,m24,m25,m26,m31,m32,m33,m34,m35,m36;

	/*variables for the signal-path, control-path and onward */
	double *ihcouttmp,*tmpgain;
	int    grd;

    double bmplace,centerfreq,gain,taubm,ratiowb,bmTaubm,fcohc,TauWBMax,TauWBMin,tauwb;
    double Taumin[1],Taumax[1],bmTaumin[1],bmTaumax[1],ratiobm[1],lasttmpgain,wbgain,ohcasym,ihcasym,delay;
	int    i,n,delaypoint,grdelay[1],bmorder,wborder;
	double wbout1,wbout,ohcnonlinout,ohcout,tmptauc1,tauc1,rsigma,wb_gain;
            
    /* Declarations of the functions used in the program */
	double C1ChirpFilt(double, double,double, int, double, double);
	double C2ChirpFilt(double, double,double, int, double, double);
    double WbGammaTone(double, double, double, int, double, double, int);

    double Get_tauwb(double, int, int, double *, double *);
	double Get_taubm(double, int, double, double *, double *, double *);
    double gain_groupdelay(double, double, double, double, int *);
    double delay_cat(double cf);
    double delay_human(double cf);

    double OhcLowPass(double, double, double, int, double, int);
    double IhcLowPass(double, double, double, int, double, int);
	double Boltzman(double, double, double, double, double);
    double NLafterohc(double, double, double, double);
	double ControlSignal(double, double, double, double, double);

    double NLogarithm(double, double, double, double);

    /* Allocate dynamic memory for the temporary variables */
	ihcouttmp = (double*) mxCalloc(nbinsAc*nrep,sizeof(double));
	mey1      = (double*) mxCalloc(nbinsAc,     sizeof(double));
	mey2      = (double*) mxCalloc(nbinsAc,     sizeof(double));
	mey3      = (double*) mxCalloc(nbinsAc,     sizeof(double));
	tmpgain   = (double*) mxCalloc(nbinsAc,     sizeof(double));
    
	/** Calculate the center frequency for the control-path wideband filter
	    from the location on basilar membrane, based on Greenwood (JASA 1990) */

	if (species==1) { // Cat frequency shift corresponding to 1.2 mm
        bmplace    = 11.9 * log10(0.80 + cf / 456.0);         // Calculate the location on basilar membrane from CF
        centerfreq = 456.0*(pow(10,(bmplace+1.2)/11.9)-0.80); // shift the center freq

        gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0);
    }

	if (species>1) { // Human frequency shift corresponding to 1.2 mm
        bmplace    = (35/2.1) * log10(1.0 + cf / 165.4);         // Calculate the location on basilar membrane from CF
        centerfreq = 165.4*(pow(10,(bmplace+1.2)/(35/2.1))-1.0); // shift the center freq
        
        gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0);
    }
    
	if(gain>60.0) gain = 60.0;  
    if(gain<15.0) gain = 15.0;
    
	/*====== Parameters for the control-path wideband filter =======*/
	bmorder = 3;
	Get_tauwb(cf,species,bmorder,Taumax,Taumin);
	taubm   = cohc*(Taumax[0]-Taumin[0])+Taumin[0];
	ratiowb = Taumin[0]/Taumax[0];
    /*====== Parameters for the signal-path C1 filter ======*/
	Get_taubm(cf,species,Taumax[0],bmTaumax,bmTaumin,ratiobm);
	bmTaubm  = cohc*(bmTaumax[0]-bmTaumin[0])+bmTaumin[0];
	fcohc    = bmTaumax[0]/bmTaubm;
    /*====== Parameters for the control-path wideband filter =======*/
	wborder  = 3;
    TauWBMax = Taumin[0]+0.2*(Taumax[0]-Taumin[0]);
	TauWBMin = TauWBMax/Taumax[0]*Taumin[0];
    tauwb    = TauWBMax+(bmTaubm-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);
	
	wbgain      = gain_groupdelay(dtAc,centerfreq,cf,tauwb,grdelay);
	tmpgain[0]  = wbgain; 
	lasttmpgain = wbgain;
  	/*===============================================================*/
    /* Nonlinear asymmetry of OHC function and IHC C1 transduction function*/
	ohcasym  = 7.0;    
	ihcasym  = 3.0;
  	/*===============================================================*/
    /*===============================================================*/
    /* PrewisARPing and related constants for the middle ear */
     fp = 1e3;  /* prewisARPing frequency 1 kHz */
     C  = TWOPI*fp/tan(TWOPI/2*fp*dtAc);
     
     if (species==1) { // Cat middle-ear filter - simplified version from Bruce et al. (JASA 2003)
         m11 = C/(C + 693.48);                    m12 = (693.48 - C)/C;            m13 = 0.0;
         m14 = 1.0;                               m15 = -1.0;                      m16 = 0.0;
         m21 = 1/(pow(C,2) + 11053*C + 1.163e8);  m22 = -2*pow(C,2) + 2.326e8;     m23 = pow(C,2) - 11053*C + 1.163e8; 
         m24 = pow(C,2) + 1356.3*C + 7.4417e8;    m25 = -2*pow(C,2) + 14.8834e8;   m26 = pow(C,2) - 1356.3*C + 7.4417e8;
         m31 = 1/(pow(C,2) + 4620*C + 909059944); m32 = -2*pow(C,2) + 2*909059944; m33 = pow(C,2) - 4620*C + 909059944;
         m34 = 5.7585e5*C + 7.1665e7;             m35 = 14.333e7;                  m36 = 7.1665e7 - 5.7585e5*C;
         megainmax=41.1405;
     };

     if (species>1) { // Human middle-ear filter - based on Pascal et al. (JASA 1998)
         m11 = 1/(pow(C,2)+5.9761e+003*C+2.5255e+007);  m12 = (-2*pow(C,2)+2*2.5255e+007);  m13 = (pow(C,2)-5.9761e+003*C+2.5255e+007);   
         m14 = (pow(C,2)+5.6665e+003*C);                m15 = -2*pow(C,2);					m16 = (pow(C,2)-5.6665e+003*C);
         m21 = 1/(pow(C,2)+6.4255e+003*C+1.3975e+008);  m22 = (-2*pow(C,2)+2*1.3975e+008);  m23 = (pow(C,2)-6.4255e+003*C+1.3975e+008);   
         m24 = (pow(C,2)+5.8934e+003*C+1.7926e+008);    m25 = (-2*pow(C,2)+2*1.7926e+008);  m26 = (pow(C,2)-5.8934e+003*C+1.7926e+008);
         m31 = 1/(pow(C,2)+2.4891e+004*C+1.2700e+009);  m32 = (-2*pow(C,2)+2*1.2700e+009);  m33 = (pow(C,2)-2.4891e+004*C+1.2700e+009);   
         m34 = (3.1137e+003*C+6.9768e+008);             m35 = 2*6.9768e+008;				m36 = (-3.1137e+003*C+6.9768e+008);
         megainmax=2;
     };

  	for (n=0;n<nbinsAc;n++) { // start of the loop    
        if (n==0)  { // Start of the middle-ear filtering section
	    	mey1[0] = m11*pAc[0];
            if (species>1) mey1[0] = m11*m14*pAc[0];
            mey2[0] = mey1[0]*m24*m21;
            mey3[0] = mey2[0]*m34*m31;
            meout   = mey3[0]/megainmax ;
        }
        else if (n==1) {
            mey1[1] = m11*(-m12*mey1[0] + pAc[1]       - pAc[0]);
            if (species>1) mey1[1] = m11*(-m12*mey1[0]+m14*pAc[1]+m15*pAc[0]);
			mey2[1] = m21*(-m22*mey2[0] + m24*mey1[1] + m25*mey1[0]);
            mey3[1] = m31*(-m32*mey3[0] + m34*mey2[1] + m35*mey2[0]);
            meout   = mey3[1]/megainmax;
		}
	    else {
            mey1[n] = m11*(-m12*mey1[n-1]  + pAc[n]         - pAc[n-1]);
            if (species>1) mey1[n]= m11*(-m12*mey1[n-1]-m13*mey1[n-2]+m14*pAc[n]+m15*pAc[n-1]+m16*pAc[n-2]);
            mey2[n] = m21*(-m22*mey2[n-1] - m23*mey2[n-2] + m24*mey1[n] + m25*mey1[n-1] + m26*mey1[n-2]);
            mey3[n] = m31*(-m32*mey3[n-1] - m33*mey3[n-2] + m34*mey2[n] + m35*mey2[n-1] + m36*mey2[n-2]);
            meout   = mey3[n]/megainmax;
		}; 	/* End of the middle-ear filtering section */   
     
		/* Control-path filter */
        wbout1 = WbGammaTone(meout,dtAc,centerfreq,n,tauwb,wbgain,wborder);
        wbout  = pow((tauwb/TauWBMax),wborder)*wbout1*10e3*__max(1,cf/5e3);
  
        ohcnonlinout = Boltzman(wbout,ohcasym,12.0,5.0,5.0); /* pass the control signal through OHC Nonlinear Function */
		ohcout = OhcLowPass(ohcnonlinout,dtAc,600,n,1.0,2);/* lowpass filtering after the OHC nonlinearity */
        
		tmptauc1 = NLafterohc(ohcout,bmTaumin[0],bmTaumax[0],ohcasym); /* nonlinear function after OHC low-pass filter */
		tauc1    = cohc*(tmptauc1-bmTaumin[0])+bmTaumin[0];  /* time -constant for the signal-path C1 filter */
		rsigma   = 1/tauc1-1/bmTaumax[0]; /* shift of the location of poles of the C1 filter from the initial positions */

		if (1/tauc1<0.0) mexErrMsgTxt("The poles are in the right-half plane; system is unstable.\n");

		tauwb = TauWBMax+(tauc1-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);

	    wb_gain = gain_groupdelay(dtAc,centerfreq,cf,tauwb,grdelay);
		
		grd = grdelay[0]; 

        if ((grd+n)<nbinsAc) tmpgain[grd+n] = wb_gain;

        if (tmpgain[n] == 0)   tmpgain[n] = lasttmpgain;	
		
		wbgain      = tmpgain[n];
		lasttmpgain = wbgain;
	 		        
        /*====== Signal-path C1 filter ======*/
		 c1filterouttmp = C1ChirpFilt(meout, dtAc, cf, n, bmTaumax[0], rsigma); /* C1 filter output */

        /*====== Parallel-path C2 filter ======*/
		 c2filterouttmp = C2ChirpFilt(meout, dtAc, cf, n, bmTaumax[0], 1/ratiobm[0]); /* parallel-filter output*/

	    /*=== Run the inner hair cell (IHC) section: NL function and then lowpass filtering ===*/
        c1vihctmp =  NLogarithm(cihc*c1filterouttmp,0.1,ihcasym,cf);
		c2vihctmp = -NLogarithm(c2filterouttmp*fabs(c2filterouttmp)*cf/10*cf/2e3,0.2,1.0,cf); /* C2 transduction output */
           
        ihcouttmp[n] = IhcLowPass(c1vihctmp+c2vihctmp,dtAc,3000,n,1.0,7);
   };  /* End of the loop */
    
    /* Stretched out the IHC output according to nrep (number of repetitions) */
   
    for(i=0;i<nbinsAc*nrep;i++) {
		ihcouttmp[i] = ihcouttmp[(int) (fmod(i,nbinsAc))];
  	};   
    
   	/* Adjust total path delay to IHC output signal */
    if (species==1)
        delay      = delay_cat(cf);
    if (species>1) {
        delay      = delay_cat(cf); /* signal delay changed back to cat function for version 5.2 */
    };
    delaypoint =__max(0,(int) ceil(delay/dtAc));    
         
    for(i=delaypoint;i<nbinsAc*nrep;i++) {        
		ihcout[i] = ihcouttmp[i - delaypoint];
  	};   

    /* Freeing dynamic memory allocated earlier */
    mxFree(ihcouttmp);
    mxFree(mey1); mxFree(mey2); mxFree(mey3);	
    mxFree(tmpgain);

} /* End of the SingleAN function */

/** Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level */
double Get_tauwb(double cf, int species, int order, double *taumax,double *taumin) {

  double Q10,bw,gain,ratio;
    
  if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
  if(species>1)  gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */

  if(gain>60.0) gain = 60.0;  
  if(gain<15.0) gain = 15.0;
   
  ratio = pow(10,(-gain/(20.0*order)));       /* ratio of TauMin/TauMax according to the gain, order */
  if (species==1) /* cat Q10 values */
  {
    Q10 = pow(10,0.4708*log10(cf/1e3)+0.4664);
  }
  if (species==2) /* human Q10 values from Shera et al. (PNAS 2002) */
  {
    Q10 = pow((cf/1000),0.3)*12.7*0.505+0.2085;
  }
  if (species==3) /* human Q10 values from Glasberg & Moore (Hear. Res. 1990) */
  {
    Q10 = cf/24.7/(4.37*(cf/1000)+1)*0.505+0.2085;
  }
  bw     = cf/Q10;
  taumax[0] = 2.0/(TWOPI*bw);
   
  taumin[0]   = taumax[0]*ratio;
  
  return 0;
}
double Get_taubm(double cf, int species, double taumax,double *bmTaumax,double *bmTaumin, double *ratio) {
    double gain,factor,bwfactor;
        
    if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
    if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/ /* older values */
    
    if(gain>60.0) gain = 60.0;  
    if(gain<15.0) gain = 15.0;

    bwfactor = 0.7;
    factor   = 2.5;

    ratio[0]  = pow(10,(-gain/(20.0*factor))); 

    bmTaumax[0] = taumax/bwfactor;
    bmTaumin[0] = bmTaumax[0]*ratio[0];     
    return 0;
}
/** Pass the signal through the signal-path C1 Tenth Order Nonlinear Chirp-Gammatone Filter */
double C1ChirpFilt(double x, double dtAc,double cf, int n, double taumax, double rsigma) { 
    static double C1gain_norm, C1initphase; 
    static double C1input[12][4], C1output[12][4];

    double ipw, ipb, rpa, pzero, rzero;
	double sigma0,fs_bilinear,CF,norm_gain,phase,c1filterout;
	int i,r,order_of_pole,half_order_pole,order_of_zero;
	double temp, dy, preal, pimg;

	COMPLEX p[11]; 
	
	/* Defining initial locations of the poles and zeros */
	/*======== setup the locations of poles and zeros =======*/
	  sigma0 = 1/taumax;
	  ipw    = 1.01*cf*TWOPI-50;
	  ipb    = 0.2343*TWOPI*cf-1104;
	  rpa    = pow(10, log10(cf)*0.9 + 0.55)+ 2000;
	  pzero  = pow(10,log10(cf)*0.7+1.6)+500;

	/*===============================================================*/     
         
     order_of_pole    = 10;             
     half_order_pole  = order_of_pole/2;
     order_of_zero    = half_order_pole;

	 fs_bilinear = TWOPI*cf/tan(TWOPI*cf*dtAc/2);
     rzero       = -pzero;
	 CF          = TWOPI*cf;
   
    if (n==0) {		  
        p[1].x = -sigma0;     
        p[1].y = ipw;
        p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;
        p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;
        p[2]   = compconj(p[1]);    p[4] = compconj(p[3]); p[6] = compconj(p[5]);
        p[7]   = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

        C1initphase = 0.0;
        for (i=1;i<=half_order_pole;i++) {
           preal     = p[i*2-1].x;
		   pimg      = p[i*2-1].y;
	       C1initphase = C1initphase + atan(CF/(-rzero))-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	   };

	    /*===================== Initialize C1input & C1output =====================*/

        for (i=1;i<=(half_order_pole+1);i++) {
		   C1input[i][3]  = 0; 
		   C1input[i][2]  = 0; 
		   C1input[i][1]  = 0;
		   C1output[i][3] = 0; 
		   C1output[i][2] = 0; 
		   C1output[i][1] = 0;
        }

	    /*===================== normalize the gain =====================*/
    
        C1gain_norm = 1.0;
        for (r=1; r<=order_of_pole; r++) {
		    C1gain_norm = C1gain_norm*(pow((CF - p[r].y),2) + p[r].x*p[r].x);
        }
   };
     
    norm_gain= sqrt(C1gain_norm)/pow(sqrt(CF*CF+rzero*rzero),order_of_zero);
	
	p[1].x = -sigma0 - rsigma;
	if (p[1].x>0.0) mexErrMsgTxt("The system becomes unstable.\n");
	p[1].y = ipw;
	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;
    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;
    p[2]   = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);
    p[7]   = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

    phase = 0.0;
    for (i=1;i<=half_order_pole;i++) {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
	       phase = phase-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	};

	rzero = -CF/tan((C1initphase-phase)/order_of_zero);
    if (rzero>0.0) mexErrMsgTxt("The zeros are in the right-half plane.\n");
	 
   /*%==================================================  */
	/*each loop below is for a pair of poles and one zero */
   /*%      time loop begins here                         */
   /*%==================================================  */
 
       C1input[1][3]=C1input[1][2]; 
	   C1input[1][2]=C1input[1][1]; 
	   C1input[1][1]= x;

        for (i=1;i<=half_order_pole;i++) {
            preal = p[i*2-1].x;
            pimg  = p[i*2-1].y;
            temp  = pow((fs_bilinear-preal),2)+ pow(pimg,2);

            dy = C1input[i][1]*(fs_bilinear-rzero) - 2*rzero*C1input[i][2] - (fs_bilinear+rzero)*C1input[i][3]
                    +2*C1output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg)
                    -C1output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);
            dy = dy/temp;

            C1input[i+1][3] = C1output[i][2]; 
            C1input[i+1][2] = C1output[i][1]; 
            C1input[i+1][1] = dy;

            C1output[i][2] = C1output[i][1]; 
            C1output[i][1] = dy;
        }

	   dy = C1output[half_order_pole][1]*norm_gain;  /* don't forget the gain term */
	   c1filterout= dy/4.0;   /* signal path output is divided by 4 to give correct C1 filter gain */
	                   
     return (c1filterout);
}  
/** Parallelpath C2 filter: same as the signal-path C1 filter with the OHC completely impaired */
double C2ChirpFilt(double xx, double dtAc,double cf, int n, double taumax, double fcohc) {
	static double C2gain_norm, C2initphase;
    static double C2input[12][4];  static double C2output[12][4];
   
	double ipw, ipb, rpa, pzero, rzero;

	double sigma0,fs_bilinear,CF,norm_gain,phase,c2filterout;
	int    i,r,order_of_pole,half_order_pole,order_of_zero;
	double temp, dy, preal, pimg;

	COMPLEX p[11]; 	
    
    /*================ setup the locations of poles and zeros =======*/
    sigma0 = 1/taumax;
    ipw    = 1.01*cf*TWOPI-50;
    ipb    = 0.2343*TWOPI*cf-1104;
    rpa    = pow(10, log10(cf)*0.9 + 0.55)+ 2000;
    pzero  = pow(10,log10(cf)*0.7+1.6)+500;
	/*===============================================================*/     
         
    order_of_pole    = 10;             
    half_order_pole  = order_of_pole/2;
    order_of_zero    = half_order_pole;

    fs_bilinear = TWOPI*cf/tan(TWOPI*cf*dtAc/2);
    rzero       = -pzero;
    CF          = TWOPI*cf;
   	    
    if (n==0) {		  
	    p[1].x = -sigma0;     
        p[1].y = ipw;
    	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;
        p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;
        p[2] = compconj(p[1]); 
        p[4] = compconj(p[3]); 
        p[6] = compconj(p[5]);
        p[7] = p[1]; 
        p[8] = p[2]; 
        p[9] = p[5]; 
        p[10]= p[6];

	    C2initphase = 0.0;
        for (i=1;i<=half_order_pole;i++) {
           preal       = p[i*2-1].x;
		   pimg        = p[i*2-1].y;
	       C2initphase = C2initphase + atan(CF/(-rzero))-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	    };

	    /*===================== Initialize C2input & C2output =====================*/

        for (i=1;i<=(half_order_pole+1);i++) {
            C2input[i][3]  = 0; 
            C2input[i][2]  = 0; 
            C2input[i][1]  = 0;
            C2output[i][3] = 0; 
            C2output[i][2] = 0; 
            C2output[i][1] = 0;
        }
    
        /*===================== normalize the gain =====================*/
    
        C2gain_norm = 1.0;
        for (r=1; r<=order_of_pole; r++) {
            C2gain_norm = C2gain_norm*(pow((CF - p[r].y),2) + p[r].x*p[r].x);
        }
    };
     
    norm_gain= sqrt(C2gain_norm)/pow(sqrt(CF*CF+rzero*rzero),order_of_zero);
    
	p[1].x = -sigma0*fcohc;
	if (p[1].x>0.0) mexErrMsgTxt("The system becomes unstable.\n");
	p[1].y = ipw;
	p[5].x = p[1].x - rpa; 
    p[5].y = p[1].y - ipb;
    p[3].x = (p[1].x + p[5].x) * 0.5; 
    p[3].y = (p[1].y + p[5].y) * 0.5;
    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);
    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

    phase = 0.0;
    for (i=1;i<=half_order_pole;i++) {
        preal = p[i*2-1].x;
        pimg  = p[i*2-1].y;
        phase = phase-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	};

	rzero = -CF/tan((C2initphase-phase)/order_of_zero);	
    if (rzero>0.0) mexErrMsgTxt("The zeros are in the right-hand plane.\n");
   
   /*%==================================================  */
   /*%      time loop begins here                         */
   /*%==================================================  */

    C2input[1][3]=C2input[1][2]; 
    C2input[1][2]=C2input[1][1]; 
    C2input[1][1]= xx;

    for (i=1;i<=half_order_pole;i++) {
        preal = p[i*2-1].x;
        pimg  = p[i*2-1].y;
        temp  = pow((fs_bilinear-preal),2)+ pow(pimg,2);
		   
        dy = C2input[i][1]*(fs_bilinear-rzero) - 2*rzero*C2input[i][2] - (fs_bilinear+rzero)*C2input[i][3]
                +2*C2output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg)
                -C2output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);
        dy = dy/temp;

        C2input[i+1][3] = C2output[i][2]; 
        C2input[i+1][2] = C2output[i][1]; 
        C2input[i+1][1] = dy;

        C2output[i][2] = C2output[i][1]; 
        C2output[i][1] = dy;
    };

    dy = C2output[half_order_pole][1]*norm_gain;
    c2filterout= dy/4.0;
    
    return (c2filterout); 
}   
/** Pass the signal through the Control path Third Order Nonlinear Gammatone Filter */
double WbGammaTone(double x,double dtAc,double centerfreq, int n, double tau,double gain,int order) {
    static double wbphase;
    static COMPLEX wbgtf[4], wbgtfl[4];

    double delta_phase,dtmp,c1LP,c2LP,out;
    int i,j;
    
    if (n==0) {
        wbphase = 0;
        for(i=0; i<=order;i++) {
            wbgtfl[i] = compmult(0,compexp(0));
            wbgtf[i]  = compmult(0,compexp(0));
        }
    }
  
    delta_phase = -TWOPI*centerfreq*dtAc;
    wbphase += delta_phase;
    
    dtmp = tau*2.0/dtAc;
    c1LP = (dtmp-1)/(dtmp+1);
    c2LP = 1.0/(dtmp+1);
    wbgtf[0] = compmult(x,compexp(wbphase));                 /* FREQUENCY SHIFT */
  
    for(j = 1; j <= order; j++)                              /* IIR Bilinear transformation LPF */
        wbgtf[j] = comp2sum(compmult(c2LP*gain,comp2sum(wbgtf[j-1],wbgtfl[j-1])), compmult(c1LP,wbgtfl[j]));
    out = REAL(compprod(compexp(-wbphase), wbgtf[order])); /* FREQ SHIFT BACK UP */
  
    for(i=0; i<=order;i++) wbgtfl[i] = wbgtf[i];
    return(out);
}
/** Calculate the gain and group delay for the Control path Filter */
double gain_groupdelay(double dtAc,double centerfreq, double cf, double tau,int *grdelay) { 
  double tmpcos,dtmp2,c1LP,c2LP,tmp1,tmp2,wb_gain;

  tmpcos = cos(TWOPI*(centerfreq-cf)*dtAc);
  dtmp2  = tau*2.0/dtAc;
  c1LP   = (dtmp2-1)/(dtmp2+1);
  c2LP   = 1.0/(dtmp2+1);
  tmp1   = 1+c1LP*c1LP-2*c1LP*tmpcos;
  tmp2   = 2*c2LP*c2LP*(1+tmpcos);
  
  wb_gain = pow(tmp1/tmp2, 1.0/2.0);
  
  grdelay[0] = (int)floor((0.5-(c1LP*c1LP-c1LP*tmpcos)/(1+c1LP*c1LP-2*c1LP*tmpcos)));

  return(wb_gain);
}
/** Calculate the delay (basilar membrane, synapse, etc. for cat) */
double delay_cat(double cf) {  
    double A0,A1,x,delay;

    A0    = 3.0;  
    A1    = 12.5;
    x     = 11.9 * log10(0.80 + cf / 456.0);      /* cat mapping */
    delay = A0 * exp( -x/A1 ) * 1e-3;
    
    return(delay);
}
/* Calculate the delay (basilar membrane, synapse, etc.) for human, based
        on Harte et al. (JASA 2009) */
double delay_human(double cf) {  
    double A,B,delay;

    A    = -0.37;  
    B    = 11.09/2;
    delay = B * pow(cf * 1e-3,A)*1e-3;
    
    return(delay);
}
/* Get the output of the OHC Nonlinear Function (Boltzman Function) */
double Boltzman(double x, double asym, double s0, double s1, double x1) {
    double shift,x0,out1,out;

    shift = 1.0/(1.0+asym);  /* asym is the ratio of positive Max to negative Max*/
    x0    = s0*log((1.0/shift-1)/(1+exp(x1/s1)));
	    
    out1 = 1.0/(1.0+exp(-(x-x0)/s0)*(1.0+exp(-(x-x1)/s1)))-shift;
	out  = out1/(1-shift);

    return(out);
}  /* output of the nonlinear function, the output is normalized with maximum value of 1 */
/* Get the output of the OHC Low Pass Filter in the Control path */
double OhcLowPass(double x,double dtAc,double Fc, int n,double gain,int order) {
    static double ohc[4],ohcl[4];

    double c,c1LP,c2LP;
    int i,j;

    if (n==0) {
        for(i=0; i<(order+1);i++) {
            ohc[i] = 0;
            ohcl[i] = 0;
        }
    }    
  
    c = 2.0/dtAc;
    c1LP = ( c - TWOPI*Fc ) / ( c + TWOPI*Fc );
    c2LP = TWOPI*Fc / (TWOPI*Fc + c);
  
    ohc[0] = x*gain;
    for(i=0; i<order;i++)  ohc[i+1] = c1LP*ohcl[i+1] + c2LP*(ohc[i]+ohcl[i]);
    for(j=0; j<=order;j++) ohcl[j]  = ohc[j];

    return(ohc[order]);
}
/* Get the output of the IHC Low Pass Filter  */
double IhcLowPass(double x,double dtAc,double Fc, int n,double gain,int order) {
    static double ihc[8],ihcl[8];
  
    double C,c1LP,c2LP;
    int i,j;

    if (n==0) {
        for(i=0; i<(order+1);i++) {
            ihc[i] = 0;
            ihcl[i] = 0;
        }
    }     
  
    C = 2.0/dtAc;
    c1LP = ( C - TWOPI*Fc ) / ( C + TWOPI*Fc );
    c2LP = TWOPI*Fc / (TWOPI*Fc + C);
  
    ihc[0] = x*gain;
    for(i=0; i<order;i++)  ihc[i+1] = c1LP*ihcl[i+1] + c2LP*(ihc[i]+ihcl[i]);
    for(j=0; j<=order;j++) ihcl[j]  = ihc[j];
    return(ihc[order]);
}
/* Get the output of the Control path using Nonlinear Function after OHC */
double NLafterohc(double x,double taumin, double taumax, double asym) {    
	double R,dc,R1,s0,x1,out,minR;

	minR = 0.05;
    R  = taumin/taumax;
    
	if(R<minR) minR = 0.5*R;
    else       minR = minR;
    
    dc = (asym-1)/(asym+1.0)/2.0-minR;
    R1 = R-minR;

    /* This is for new nonlinearity */
    s0 = -dc/log(R1/(1-minR));
	
    x1  = fabs(x);
    out = taumax*(minR+(1.0-minR)*exp(-x1/s0));
	if (out<taumin) out = taumin; 
    if (out>taumax) out = taumax;
    return(out);
}
/* Get the output of the IHC Nonlinear Function (Logarithmic Transduction Functions) */
double NLogarithm(double x, double slope, double asym, double cf) {
	double corner,strength,xx,splx,asym_t;
	    
    corner    = 80; 
    strength  = 20.0e6/pow(10,corner/20);
            
    xx = log(1.0+strength*fabs(x))*slope;
    
    if(x<0)	{
        splx   = 20*log10(-x/20e-6);
		asym_t = asym -(asym-1)/(1+exp(splx/5.0));
		xx = -1/asym_t*xx;
	};   
    return(xx);
}

/* ================================ */
/* ======  Synapse functions  ===== */
/* ================================ */

double Synapse(double *ihcout, double dtAc, double cf, int nbinsAc, int nrep, double spont, double noiseType, double implnt, double sampFreq, double *Sout) {

    /* --- Initalize Variables --- */
    int z;
    int b;
    int    resamp     = (int) ceil( 1/(dtAc*sampFreq) );
    int    delaypoint = (int) floor( 7500/(cf/1e3) );
    double incr       = 0.0; 
        
    /* Power-law function */
    const double alpha1   = 1.5e-6*100e3; // original
    const double beta1    = 5.0e-4;
    double I1 = 0.0;
    const double alpha2   = 1.0e-2*100e3; // original
    const double beta2    = 1.0e-1;
    double I2 = 0.0;
    const double binwidth = 1.0/sampFreq;
    int    k,j,indx,i;
    
    /* mapping function from IHCOUT to input to the PLA */
    const double cfslope   = pow(spont, 0.19) * pow(10, -0.87);
    const double cfconst   = 0.1 * pow(log10(spont), 2) + 0.56*log10(spont) - 0.84;
    const double cfsat     = pow(10, (cfslope*8965.5/1e3 + cfconst));
    const double cf_factor = __min(cfsat, pow(10, cfslope*cf/1e3 + cfconst))*2.0;
    const double multFac   = __max(2.95*__max(1.0, 1.5 - spont/100), 4.3 - 0.2*cf/1e3);
    
    double *sout1      = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *sout2      = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *synSampOut = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *powerLawIn = (double*) mxCalloc( (long) ceil(  nbinsAc*nrep + 3*delaypoint                 ), sizeof(double) );
    double *mappingOut = (double*) mxCalloc( (long) ceil(  nbinsAc*nrep                                ), sizeof(double) );
    double *TmpSyn     = (double*) mxCalloc( (long) ceil(  nbinsAc*nrep + 2*delaypoint                 ), sizeof(double) );
    double *m1         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *m2         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *m3         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *m4         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *m5         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *n1         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *n2         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    double *n3         = (double*) mxCalloc( (long) ceil( (nbinsAc*nrep + 2*delaypoint) *dtAc*sampFreq ), sizeof(double) );
    
    mxArray	*randInputArray[5];
    mxArray *randOutputArray[1];
    double  *randNums;
    
    mxArray	*IhcInputArray[3];
    mxArray *IhcOutputArray[1];
    double  *sampIHC;
    double  *ihcDims;
  
    /*----------------------------------------------------------*/
    /*------- Generating a random sequence ---------------------*/
    /*----------------------------------------------------------*/
    randInputArray[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(randInputArray[0])= ceil((nbinsAc*nrep+2*delaypoint)*dtAc*sampFreq);
    randInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(randInputArray[1])= 1/sampFreq;
    randInputArray[2] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(randInputArray[2])= 0.9;          // Hurst index
    randInputArray[3] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(randInputArray[3])= noiseType;    // fixed or variable fGn
    randInputArray[4] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(randInputArray[4])= spont;        // high, medium, or low
    
    // mexCallMATLAB(1, randOutputArray, 5, randInputArray, "ffGn_ac");  randNums = mxGetPr(randOutputArray[0]);
    mexCallMATLAB(1, randOutputArray, 5, randInputArray, "Model.ffGn_ac");  randNums = mxGetPr(randOutputArray[0]);

    /*----------------------------------------------------------*/
    /*----------------------------------------------------------*/
    k = 0;
    for (indx=0; indx<nbinsAc*nrep; ++indx) {
        mappingOut[k] = pow(10, (0.9*log10(fabs(ihcout[indx])*cf_factor)) + multFac);
        if (ihcout[indx]<0) mappingOut[k] = - mappingOut[k];
        k=k+1;
    }
    for (k=0; k<delaypoint; k++)                                      powerLawIn[k] = mappingOut[0]            + 3.0*spont;
    for (k=delaypoint; k<nbinsAc*nrep+delaypoint; k++)                powerLawIn[k] = mappingOut[k-delaypoint] + 3.0*spont;
    for (k=nbinsAc*nrep+delaypoint; k<nbinsAc*nrep+3*delaypoint; k++) powerLawIn[k] = powerLawIn[k-1]          + 3.0*spont;
    
    /*----------------------------------------------------------*/
    /*------ Downsampling to sampFreq (Low) sampling rate ------*/
    /*----------------------------------------------------------*/
    IhcInputArray[0] = mxCreateDoubleMatrix(1, k, mxREAL);
    ihcDims          = mxGetPr(IhcInputArray[0]);
    for (i=0;i<k;++i) {
        ihcDims[i] = powerLawIn[i];
    }
    IhcInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(IhcInputArray[1]) = 1;
    IhcInputArray[2] = mxCreateDoubleMatrix(1, 1, mxREAL);  *mxGetPr(IhcInputArray[2]) = resamp;
    mexCallMATLAB(1, IhcOutputArray, 3, IhcInputArray, "resample");
    sampIHC = mxGetPr(IhcOutputArray[0]);
    
    mxFree(powerLawIn); 
    mxFree(mappingOut);

    /*----------------------------------------------------------*/
    /*----- Running Power-law Adaptation -----------------------*/
    /*----------------------------------------------------------*/
    k = 0;
    for (indx=0; indx<floor((nbinsAc*nrep+2*delaypoint)*dtAc*sampFreq); indx++) {
        sout1[k]  = __max( 0, sampIHC[indx] + randNums[indx] - alpha1*I1);
        sout2[k]  = __max( 0, sampIHC[indx] - alpha2*I2);
        
        if (implnt==1) { // ACTUAL Implementation
            I1 = 0; 
            I2 = 0;
            for (j=0; j<k+1; ++j) {
                I1 += (sout1[j])*binwidth/((k-j)*binwidth + beta1);
                I2 += (sout2[j])*binwidth/((k-j)*binwidth + beta2);
            }
        } // end of actual
        
        if (implnt==0) { // APPROXIMATE Implementation
            if (k==0) {
                n1[k] = 1.0e-3*sout2[k];
                n2[k] = n1[k]; 
                n3[0] = n2[k];
            } else if (k==1) {
                n1[k] =  1.992127932802320*n1[k-1] + 1.0e-3*(sout2[k] - 0.994466986569624*sout2[k-1]);
                n2[k] =  1.999195329360981*n2[k-1] + n1[k] - 1.997855276593802*n1[k-1];
                n3[k] = -0.798261718183851*n3[k-1] + n2[k] + 0.798261718184977*n2[k-1];
            } else {
                n1[k] = 1.992127932802320*n1[k-1] - 0.992140616993846*n1[k-2] + 1.0e-3*(sout2[k] - 0.994466986569624*sout2[k-1] + 0.000000000002347*sout2[k-2]);
                n2[k] = 1.999195329360981*n2[k-1] - 0.999195402928777*n2[k-2] + n1[k] - 1.997855276593802*n1[k-1] + 0.997855827934345*n1[k-2];
                n3[k] =-0.798261718183851*n3[k-1] - 0.199131619873480*n3[k-2] + n2[k] + 0.798261718184977*n2[k-1] + 0.199131619874064*n2[k-2];
            }
            I2 = n3[k];
            
            if (k==0) {
                m1[k] = 0.2*sout1[k];
                m2[k] = m1[k];	
                m3[k] = m2[k];
                m4[k] = m3[k];	
                m5[k] = m4[k];
            } else if (k==1) {
                m1[k] = 0.491115852967412*m1[k-1] + 0.2*(sout1[k] - 0.173492003319319*sout1[k-1]);
                m2[k] = 1.084520302502860*m2[k-1] + m1[k] - 0.803462163297112*m1[k-1];
                m3[k] = 1.588427084535629*m3[k-1] + m2[k] - 1.416084732997016*m2[k-1];
                m4[k] = 1.886287488516458*m4[k-1] + m3[k] - 1.830362725074550*m3[k-1];
                m5[k] = 1.989549282714008*m5[k-1] + m4[k] - 1.983165053215032*m4[k-1];
            } else {
                m1[k] = 0.491115852967412*m1[k-1] - 0.055050209956838*m1[k-2] + 0.2*(sout1[k]- 0.173492003319319*sout1[k-1]+ 0.000000172983796*sout1[k-2]);
                m2[k] = 1.084520302502860*m2[k-1] - 0.288760329320566*m2[k-2] + m1[k] - 0.803462163297112*m1[k-1] + 0.154962026341513*m1[k-2];
                m3[k] = 1.588427084535629*m3[k-1] - 0.628138993662508*m3[k-2] + m2[k] - 1.416084732997016*m2[k-1] + 0.496615555008723*m2[k-2];
                m4[k] = 1.886287488516458*m4[k-1] - 0.888972875389923*m4[k-2] + m3[k] - 1.830362725074550*m3[k-1] + 0.836399964176882*m3[k-2];
                m5[k] = 1.989549282714008*m5[k-1] - 0.989558985673023*m5[k-2] + m4[k] - 1.983165053215032*m4[k-1] + 0.983193027347456*m4[k-2];
            }
            I1 = m5[k];
        } // end of approximate implementation
        
        synSampOut[k] = sout1[k] + sout2[k];
        k = k+1;
    } // end of all samples

    mxFree(sout1); 
    mxFree(sout2);
    mxFree(m1); 
    mxFree(m2); 
    mxFree(m3); 
    mxFree(m4); 
    mxFree(m5); 
    mxFree(n1); 
    mxFree(n2); 
    mxFree(n3);
    
    /*-------------------------------------------------------------*/
    /*----- Upsampling to original (High 100 kHz) sampling rate ---*/
    /*-------------------------------------------------------------*/
    for(z=0; z<k-1; ++z) {
        incr = (synSampOut[z+1]-synSampOut[z])/resamp;
        for(b=0; b<resamp; ++b) {
            TmpSyn[z*resamp+b] = synSampOut[z]+ b*incr;
        }
    }
    for (i=0;i<nbinsAc*nrep;++i) {
        Sout[i] = TmpSyn[i+delaypoint];
    }
    
    mxFree(synSampOut); 
    mxFree(TmpSyn);
    mxDestroyArray(randInputArray[0]); 
    mxDestroyArray(randInputArray[1]); 
    mxDestroyArray(randInputArray[2]); 
    mxDestroyArray(randInputArray[3]);
    mxDestroyArray(randInputArray[4]);
    mxDestroyArray(randOutputArray[0]);
    mxDestroyArray(IhcInputArray[0]);  
    mxDestroyArray(IhcInputArray[1]); 
    mxDestroyArray(IhcInputArray[2]);
    mxDestroyArray(IhcOutputArray[0]);     
    
    return ( (long) ceil(nbinsAc*nrep) );
}

void calcShapedNoise(double *noisePtr, int length, double dt, double alpha, double noiseType, double sigma) {
    
    mxArray	*noiseIn[3], *noiseOut[1] ; // input and output arrays for MATLAB noise generator
    
    int      i;

    noiseIn[0] = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(noiseIn[0]) = length;           // length of the noise sequence
    noiseIn[1] = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(noiseIn[1]) = alpha;            // NoiseAlpha
    noiseIn[2] = mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(noiseIn[2]) = 1;

    mexCallMATLAB(1, noiseOut, 3, noiseIn, "Model.oneonfnoise"); // call to MATLAB noise generator

        
    for (i=0; i<length; i++) {
        noisePtr[i] = sigma * (mxGetPr(noiseOut[0]))[i];
    }

    mxDestroyArray(noiseIn[0]);
    mxDestroyArray(noiseIn[1]);
    mxDestroyArray(noiseIn[2]);
    mxDestroyArray(noiseOut[0]);

}