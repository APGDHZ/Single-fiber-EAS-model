% model_EAS2021 - single fiber EAS model by D. Kipping & W. Nogueira (2021)
%
% [psthA,psthE,meanrate,varrate,Sout,trd_vec,trel_vec,vihc,vmem,isub,isupra,inoise,ives] 
%       = Model.model_EAS2021(stimA,stimE,nrep,reptime,cf,cohc,cihc,noiseType,implnt,spont,tabs,trel,Cper,Ccen,sigmaPer,sigmaCen,coupling) ;
%
% Input:
%   - stimA         acoustic stimulus in [pa], sampled at 1e5 Hz (!)
%   - stimE         electric stimulus in [A] , sampled at 1e6 Hz (!)
%   - nrep          number of repetitions for the PSTHs
%   - reptime       time between stimulus repetitions in s
%   - cf            AS model: characteristic frequency of the fiber in Hz
%   - cohc          AS model: OHC scaling factor (1 is normal OHC function; 0 is complete OHC dysfunction)
%   - cihc          AS model: IHC scaling factor (1 is normal IHC function; 0 is complete IHC dysfunction)
%   - noiseType     AS model: "1" for variable fGn and "0" for fixed (frozen) fGn
%   - implnt        AS model: "0" for approx. and "1" for actual implementation of the power-law functions
%   - spont         AS model: spontaneous firing rate in /s
%   - tabs          AS model: absolute refractory period in s (time constants of the ES model are scaled according to this variable)
%   - trel          AS model: baselines mean relative refractory period in s (time constants of the ES model are scaled according to this variable)
%   - Cper          ES model: membrane capacitance of the peripheral neuron in F
%   - Ccen          ES model: membrane capacitance of the central neuron in F
%   - sigmaPer      ES model: noise current amplitude of the peripheral neuron (in amperes)
%   - sigmaCen      ES model: noise current amplitude of the central neuron (in amperes)
%   - coupling      EAS model coupling variant: "0" for uncoupled EAS model; "1" for coupled EAS model; "2" for alternative EAS model
%
% Output:
%   - psthA         PSTH output of the AS model, 1 time bin = 1e-5 s
%   - psthE         PSTH output of the ES model, 1 time bin = 1e-5 s
%   - meanrate      AS model: analytical estimate of the mean firing rate in /s for each time bin
%   - varrate       AS model: analytical estimate of the variance in firing rate in /s for each time bin
%   - Sout          AS model: synapse output rate in /s  for each time bin (before the effects of redocking are considered)
%   - trd_vec       AS model: vector of the mean redocking time in s for each time bin
%   - trel_vec      AS model: vector of the mean relative refractor period in s for each time bin
%   - vihc          AS model: inner hair cell (IHC) relative transmembrane potential (in volts)
%   - vmem          ES model: transmembrane potentials of peripheral/central neuron (in volts)
%   - isub          ES model: subthreshold adaptation currents of peripheral/central neuron (in amperes)
%   - isupra        ES model: suprathreshold adaptation currents of peripheral/central neuron (in amperes)
%   - inoise        ES model: noise currents of peripheral/central neuron (in amperes)
%   - ives          Coupled EAS model: neurotransmitter release-triggered excitatory current to the peripheral neuron (in amperes)
%
% If more than 1 repetition is simulated (nrep>1), the outputs vmem,isub,...,ives of the ES model correspond to the last repetition.
%
% %%% D. Kipping, Aug. 2021 %%%
%
% - revised by DK in Jul 2022