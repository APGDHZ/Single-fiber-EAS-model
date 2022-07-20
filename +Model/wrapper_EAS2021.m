function [nspikes,psthA,psthE,meanrate,varrate,Sout,trd_vec,trel_vec,vihc,vmem,isub,isupra,inoise,ives] = wrapper_EAS2021(stimA,stimE,varargin)
% Wrapper function for the EAS2021 model by D. Kipping & W. Nogueira.
%
% Required input:
%   - stimA         acoustic stimulus in [pa], sampled at 1e5 Hz (!)
%   - stimE         electric stimulus in [A] , sampled at 1e6 Hz (!)
%
% Output:
%   - nspikes       total number of E+A spikes in all repetitions, nspikes = sum(psthA) + sum(psthE)
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
% (if more than 1 repetition is simulated (see below), the outputs vmem,isub,...,ives of the ES model correspond to the last of N runs.)
%
% Optional inputs (name-value-pairs):
%   - nrep          number of repetitions for the PSTHs [default=1] 
%   - reptime       time between stimulus repetitions in s [default=stimulus length] 
%   - cf            AS model: characteristic frequency of the fiber in Hz [default=1000] 
%   - cohc          AS model: OHC scaling factor (1 is normal OHC function; 0 is complete OHC dysfunction) [default=1] 
%   - cihc          AS model: IHC scaling factor (1 is normal IHC function; 0 is complete IHC dysfunction) [default=1] 
%   - noiseType     AS model: "1" for variable fGn and "0" for fixed (frozen) fGn [default=1] 
%   - implnt        AS model: "0" for approx. and "1" for actual implementation of the power-law functions [default=0] 
%   - spont         AS model: spontaneous firing rate in /s [default=1000] 
%   - tabs          AS model: absolute refractory period in s (time constants of the ES model are scaled according to this variable) [default=450.0e-6] 
%   - trel          AS model: baselines mean relative refractory period in s (time constants of the ES model are scaled according to this variable) [default=512.5e-6] 
%   - Cper          ES model: membrane capacitance of the peripheral neuron in F [default=856.96e-9] 
%   - Ccen          ES model: membrane capacitance of the central neuron in F [default=1772.4e-9] 
%   - sigmaPer      ES model: noise current amplitude of the peripheral neuron (in amperes) [default=8.697e-6] 
%   - sigmaCen      ES model: noise current amplitude of the central neuron (in amperes) [default=11.890e-6] 
%   - coupling      EAS model coupling variant: "0" for uncoupled EAS model; "1" for coupled EAS model; "2" for alternative EAS model [default=1] 
%   - model         function handle of the EAS model .mex file [default=@Model.model_EAS2021] 
%
% %%% D. Kipping, Aug. 2021 %%%
%
% - revised by DK in Jul 2022

fsAc = 1e5 ;
fsEl = 1e6 ;

default_nrep      = 1 ;
default_reptime   = max( length(stimA)/fsAc , length(stimE)/fsEl ) ;
default_cf        = 1e3 ;
default_cohc      = 1 ;
default_cihc      = 1 ;
default_noiseType = 1 ;
default_implnt    = 0 ;
default_spont     = 1e-3 ;
default_tabs      = 450.0e-6 ;
default_trel      = 512.5e-6 ;
default_cper      = 856.96e-9 ;
default_ccen      = 1772.4e-9 ;
default_coupling  = 1 ;
default_model     = @Model.model_EAS2021 ;
default_sigmaPer  =  8.697e-6 ;
default_sigmaCen  = 11.890e-6 ;

p = inputParser ;
addRequired(p, 'stimA')
addRequired(p, 'stimE')
addParameter(p, 'nrep', default_nrep)
addParameter(p, 'reptime', default_reptime)
addParameter(p, 'cf', default_cf)
addParameter(p, 'cohc', default_cohc)
addParameter(p, 'cihc', default_cihc)
addParameter(p, 'noiseType', default_noiseType)
addParameter(p, 'implnt', default_implnt)
addParameter(p, 'spont', default_spont)
addParameter(p, 'tabs', default_tabs)
addParameter(p, 'trel', default_trel)
addParameter(p, 'Cper', default_cper)
addParameter(p, 'Ccen', default_ccen)
addParameter(p, 'coupling', default_coupling)
addParameter(p, 'model', default_model)

addParameter(p, 'sigmaPer', default_sigmaPer)
addParameter(p, 'sigmaCen', default_sigmaCen)

parse(p,stimA,stimE,varargin{:});

[psthA,psthE,meanrate,varrate,Sout,trd_vec,trel_vec,vihc,vmem,isub,isupra,inoise,ives] = ...
    feval(p.Results.model, p.Results.stimA, p.Results.stimE, p.Results.nrep, p.Results.reptime, ...
    p.Results.cf, p.Results.cohc, p.Results.cihc, p.Results.noiseType, p.Results.implnt, p.Results.spont, ...
    p.Results.tabs, p.Results.trel, p.Results.Cper, p.Results.Ccen, p.Results.sigmaPer, p.Results.sigmaCen, p.Results.coupling) ;

nspikes = sum(psthA) + sum(psthE) ;

end

