function [FE,varFE,lat,jit,VS,FF,psth,nruns] = computeFE(stimA,stimE,nrep,reptime,baselineSR,varargin)
%

% for single pulses, use nPulses=1 and pulseRate=0 .
default_modelMode       = 'EAS';
default_nPulses         = 1 ;
default_pulseRate       = 1/reptime ;
default_ANFparameters   = struct() ;
default_maskerMode      = 'none' ;
default_maskerWindow    = [0 reptime] ;
default_analysisWindows = {[0 reptime]} ;    

p = inputParser ;
addRequired(p, 'stimA'     , @(x) ( all(isnumeric(x)) || isa(x,'function_handle') ))
addRequired(p, 'stimE'     , @(x) all(isnumeric(x)))
addRequired(p, 'nrep'      , @(x) (numel(x)==1 & isnumeric(x) & (x>=1) & (mod(x,1)==0)))
addRequired(p, 'reptime'   , @(x) (numel(x)==1 & isnumeric(x) & (x>0)))
addRequired(p, 'baselineSR', @(x) (numel(x)==1 & isnumeric(x)))
addParameter(p, 'modelMode'      , default_modelMode      , @(x) ismember(x,{'electric','EAS','acoustic'}))
addParameter(p, 'nPulses'        , default_nPulses        , @(x) (numel(x)==1 & isnumeric(x) & (x>=1) & (mod(x,1)==0)))
addParameter(p, 'pulseRate'      , default_pulseRate      , @(x) (numel(x)==1 & isnumeric(x) & (x>=0)))
addParameter(p, 'ANFparameters'  , default_ANFparameters  , @(x) isstruct(x))
addParameter(p, 'maskerMode'     , default_maskerMode     , @(x) ismember(x,{'none','subthreshold_masker','suprathreshold_masker'}))
addParameter(p, 'maskerWindow'   , default_maskerWindow   , @(x) (numel(x)==2 & all(isnumeric(x)) & (0<=x(1)) & (x(1)<x(2)) & (x(2)<=reptime) ))
addParameter(p, 'analysisWindows', default_analysisWindows, @(x) (iscell(x)) )
parse(p,stimA,stimE,nrep,reptime,baselineSR,varargin{:});

[FE,varFE,lat,jit,VS,FF,psth,nruns] = stats(p.Results.stimA, p.Results.stimE, ...
        p.Results.nrep, p.Results.reptime, p.Results.baselineSR, ...
        p.Results.modelMode, p.Results.nPulses, p.Results.pulseRate, ...
        p.Results.ANFparameters, p.Results.maskerMode, p.Results.maskerWindow, ...
        p.Results.analysisWindows) ;

end

function [FE,varFE,lat,jit,VS,FF,psth,nruns] = stats(stimA,stimE,nrep,reptime,baselineSR,modelMode,nPulses,pulseRate,ANFparameters,maskerMode,maskerWindow,analysisWindows)
    nW = numel(analysisWindows) ;
    if numel(pulseRate)==1
        pulseRate = repmat(pulseRate,nW,1) ;
    end
    analysisWindowsTotalLength = zeros(1,nW) ;
    for i=1:nW
        analysisWindowsTotalLength(i) = sum(diff(analysisWindows{i},1,2),1) ;
    end
    if numel(nPulses)==1 && nW > 1
        for i=1:nW
            nPulses(i) = max(0 , analysisWindowsTotalLength(i) * pulseRate(i)) ;
        end
    end
    indW  = cell(nW,1) ; 

    nProbe = nan(nW,nrep) ;
    nruns  = 0 ;
    
    for i=1:nrep

        [~,ApsthTmp,EpsthTmp] = Model.wrapper_EAS2021(stimA(),stimE,ANFparameters,'reptime',reptime) ;
        switch modelMode
            case 'electric'
                % only use output of the electric model
                psthTmp = EpsthTmp ;
            case 'acoustic'
                % only use output of the acoustic model
                psthTmp = ApsthTmp ;
            case 'EAS'
                % use outputs of both models
                psthTmp = ApsthTmp + EpsthTmp ;
            otherwise
                error('Mode must be either ''electric'', ''acoustic'' or ''EAS''')
        end
        
        if i==1
            % initialize time array
            t = (1:length(psthTmp))*1e-5 ;
            % initialize logical arrays for masker window
            iMasker = ( t >= maskerWindow(1) ) & ( t <= maskerWindow(2) ) ;
            % initialize psth arrays for every analysis window:
            indW = false(nW,length(t)) ;
            for j=1:nW
                analysisWindowTmp = analysisWindows{j} ;
                for k=1:size(analysisWindows{j},1) 
                    indW(j,:) = indW(j,:) |  ( (t>=analysisWindowTmp(k,1)) & (t<analysisWindowTmp(k,2)) ) ;
                end
            end
            psth = zeros(size(psthTmp)) ;
        end
        
        % count spikes in the masker / probe time window
        nMaskerTmp = sum(psthTmp(iMasker)) ;
        switch maskerMode
            case 'none'
                % count every run
                for j=1:nW
                    nProbe(j,i) = sum(psthTmp(indW(j,:))) ; % for calculation of FE and var(FE)
                end
                psth  = psth + psthTmp ; % for calculation of latency and jitter
                nruns = nruns + 1 ;
            case 'subthreshold_masker'
                % only count runs without spikes in the masker time window:
                if nMaskerTmp==0
                    for j=1:nW
                        nProbe(j,i) = sum(psthTmp(indW{j})) ; % for calculation of FE and var(FE)
                    end
                    psth  = psth + psthTmp ; % for calculation of latency and jitter
                    nruns = nruns + 1 ;                    
                end
            case 'suprathreshold_masker'
                % only count runs with spikes in the masker time window:
                if nMaskerTmp>0
                    for j=1:nW
                        nProbe(j,i) = sum(psthTmp(indW(j,:))) ; % for calculation of FE and var(FE)
                    end
                    psth  = psth + psthTmp ; % for calculation of latency and jitter
                    nruns = nruns + 1 ;                    
                end
            otherwise
                warning('Unknown maskerMode input.')
        end
    end
    
    FEtmps              = nan(nW,nrep) ;
    UncorrectedRateTmps = nan(nW,nrep) ;
    for j=1:nW
        UncorrectedRateTmps(j,:) = nProbe(j,:) / analysisWindowsTotalLength(j) ;
        FEtmps(j,:)              = ( UncorrectedRateTmps(j,:) - baselineSR ) * analysisWindowsTotalLength(j) / nPulses(j) ;
    end
    
    % compute firing efficiency and variance:
    FE    = nanmean(FEtmps,2) ;
    varFE = nanvar( FEtmps,[],2) ;
    
    % preallocate output arrays:
    lat = nan(nW,1) ;
    jit = nan(nW,1) ;
    VS  = nan(nW,1) ;
    FF  = nan(nW,1) ;
    
    if nargout > 2
        for j=1:nW
            tSpikeW = repelem(t(indW(j,:)), psth(indW(j,:))) ;
            tSpikeW = mod(tSpikeW, 1/pulseRate(j)) ;

            % compute latency and jitter
            lat(j,1) = nanmean(tSpikeW) ;
            jit(j,1) = nanstd( tSpikeW) ;

            % compute vector strength and fano factor
            if numel(tSpikeW) > 0
                spikePhasesW = 2*pi*pulseRate(j)*tSpikeW ; 
                VS(j,1) = sqrt( (sum(sin(spikePhasesW))).^2 + (sum(cos(spikePhasesW))).^2 ) / numel(tSpikeW) ;
                FF(j,1) = nanvar( UncorrectedRateTmps(j,:) ) / nanmean( UncorrectedRateTmps(j,:) ) ;
            else
                VS(j,1) = nan ;
                FF(j,1) = nan ;
            end
        end
    end
    
    
end