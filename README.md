# Single-fiber-EAS-model

This is the code for the EAS single fiber model described in:
  
   Kipping, D., and Nogueira, W. (under review). "A computational model of a single auditory nerve fiber for electric-acoustic stimulation", 
   to appear in JARO - Journal of the Association for Research in Otolaryngology
 
Please cite this paper if you publish any research results obtained with this code or any modified versions of this code.
 
## Acknowledgements

This implementation of the EAS model is based on code of I. Bruce et al. and S. Joshi et al. described in:
 
   Bruce, I.C., Erfani, Y., Zilany, M.S.A. (2018). "A phenomenological model of the synapse between the inner hair cell and auditory nerve:
   Implications of limited neurotransmitter release sites." Hear. Res. 360, 40–54. https://doi.org/10.1016/j.heares.2017.12.016

   Joshi, S.N., Dau, T., Epp, B. (2017). "A Model of Electrically Stimulated Auditory Nerve Fiber Responses with Peripheral and Central Sites of Spike 
   Generation." JARO - J. Assoc. Res. Otolaryngol. 18, 323–342. https://doi.org/10.1007/s10162-016-0608-2

## Instructions

The EAS model is implemented in c code which can be compiled as a Matlab .mex file. In order to compile the code, run 
  
  >>  Model.mexEASmodel 
  
from the Matlab command prompt. This will create a .mex file in the "Model" folder which can be called like a Matlab function:
  
  >>  [...] = Model.model_EAS2021(...) ;

This raw model has a large number of input and output arguments, which you can see by typing "help Model.model_EAS2021".
For a more convenient use, we provide a wrapper function that sets default values for optional input arguments:

  >>  [...] = Model.wrapper_EAS2021(...) ;
  
  >>  help Model.wrapper_EAS2021

The other functions in the "Model" folder are subroutines of the model and should not be changed.

The "Library" folder contains useful scripts for performing experiments:

  >>  Library.generateANpopulation_EAS2021:  generates a randomized parameter set for a fiber population as described in the manuscript
  
  >>  Library.computeFE:                     simulate statistics at a given electric stimulation level
  
  >>  Library.findElectricThreshold:         estimate electric threshold and simulate statistics at threshold
  
  >>  Library.constructNoise:                construct acoustic noise used for the Miller et al. (2009) study
  
  >>  Library.constructPulseTrain:           construct electric pulse train
  
  >>  Library.fitaudiogram2:                 fit cohc and cihc to a given audiogram (script from the BEZ2018 publication of the AS model)

For instance, to simulate threshold and statistics at threshold for a biphasic pulse using 100 stimulus repetitions at each stimulation level 
and a repetition time of 10 ms together with default settings for all other model parameters, use
  
  >>  [mu,RS,Lthr,Jthr] = Library.findElectricThreshold(stimE, 100, 10e-3) ;
 
where stimE is the biphasic pulse with an amplitude of 1A sampled at 1e6 Hz, mu is the estimated threshold in A, RS is the estimated relative spread, 
and Lthr & Jthr are the estimated mean spike latency and jitter (in s) at threshold.

************************
Daniel Kipping, Aug 2021
