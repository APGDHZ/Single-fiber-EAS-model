function [mu,RS,Lthr,Jthr,levels,FEs,varFEs,Lats,Jits,VS,FF,baselineSR] = findElectricThreshold(stimE,nrep,reptime,varargin)

default_nPulses         = 1 ; % for single pulse stimuli
default_pulseRate       = 1/reptime ; % for single pulse stimuli
default_ANFparameters   = struct();
default_nrepStat        = nrep ;
default_startLevel      = 0.5e-3 ;
default_electricMasker  = zeros(size(stimE)) ;
default_acousticMasker  = zeros(1,max(1,floor(length(stimE)/10))) ;
default_maskerWindow    = [0 reptime] ;
default_analysisWindows = {[0 reptime]} ;
default_maskerMode      = 'none' ;
default_modelMode       = 'EAS' ;  

p = inputParser ;
addRequired(p, 'stimE')
addRequired(p, 'nrep')
addRequired(p, 'reptime')
addParameter(p, 'nPulses', default_nPulses)
addParameter(p, 'pulseRate', default_pulseRate)
addParameter(p, 'ANFparameters', default_ANFparameters)
addParameter(p, 'nrepStat', default_nrepStat)
addParameter(p, 'startLevel', default_startLevel)
addParameter(p, 'electricMasker', default_electricMasker)
addParameter(p, 'acousticMasker', default_acousticMasker)
addParameter(p, 'maskerWindow', default_maskerWindow)
addParameter(p, 'analysisWindows', default_analysisWindows)
addParameter(p, 'maskerMode', default_maskerMode)
addParameter(p, 'modelMode', default_modelMode)
parse(p,stimE,nrep,reptime,varargin{:});

parameters = {'ANFparameters'  , p.Results.ANFparameters, ...
              'modelMode'      , p.Results.modelMode, ...
              'nPulses'        , p.Results.nPulses, ...
              'pulseRate'      , p.Results.pulseRate, ...
              'maskerWindow'   , p.Results.maskerWindow, ...
              'analysisWindows', p.Results.analysisWindows, ...
              'maskerMode'     , p.Results.maskerMode} ;

[mu,RS,Lthr,Jthr,levels,FEs,varFEs,Lats,Jits,VS,FF,baselineSR] = findthreshold(...
    p.Results.acousticMasker, p.Results.electricMasker, p.Results.stimE, ...
    p.Results.nrep, p.Results.reptime, p.Results.nrepStat, ...
    p.Results.startLevel, parameters) ;
    
end

function [mu,RS,Lthr,Jthr,levels,FEs,varFEs,Lats,Jits,VS,FF,baselineSR] = findthreshold(maskerA,maskerE,probeE,nrep,reptime,nrepStat,startLevel,parameters)

% parameters for threshold algorithm:
step_db     = 2 ;       % level step for FE curve in dB 
nSat        = 5 ;
dFE         = 0.25 ;
nStop       = 25 ;      % maximum iterations for up & down sweeps
nLevels     = 10 ;      % 

% measure spont rate (5 s)
ind1 = find(contains(parameters(1:2:end), 'ANFparameters')) ;
ind2 = find(contains(parameters(1:2:end), 'modelMode')) ;
ind3 = find(contains(parameters(1:2:end), 'analysisWindows')) ;
parametersSR = parameters([(2*ind1-1):2*ind1 , (2*ind2-1):2*ind2]) ;
T = sum(diff(cell2mat(parameters{2*ind3}'),[],2)) ;
N = ceil(5/T) ;
baselineSR = Library.computeFE(zeros(1,10),zeros(1,100),N,T,0,parametersSR{:}) / T ;

% construct stimuli
nEp = length(probeE) ;
nEm = length(maskerE) ;
if isnumeric(maskerA)
    nA      = length(maskerA) ;
elseif isa(maskerA,'function_handle')
    nA      = length(maskerA()) ;
end
nAfinal = max([nA ceil(nEp/10) ceil(nEm/10)]) ;
fA      = @() [maskerA() zeros(1,nAfinal-nA)] ;

if nEp < 10*nAfinal
    probeE(10*nAfinal) = 0 ;
end
if nEm < 10*nAfinal
    maskerE(10*nAfinal) = 0 ;
end

% statistics at start level
levels = startLevel ;
[FEs,varFEs,Lats,Jits,VS,FF] = Library.computeFE(fA,maskerE+startLevel*probeE,nrep,reptime,baselineSR,parameters{:}) ;
nW = size(FEs,1) ;

% upward sweep
countSat  = 0 ;
countIter = 0 ;
while ( countSat < nSat )
    levels(end+1)= db2mag(+step_db)*levels(end) ;
    [FEs(:,end+1),varFEs(:,end+1),Lats(:,end+1),Jits(:,end+1),VS(:,end+1),FF(:,end+1)] = Library.computeFE(fA,maskerE+levels(end)*probeE,nrep,reptime,baselineSR,parameters{:}) ;
    if FEs(1,end) >= 1-dFE
        countSat = countSat+1 ;
    end
    countIter = countIter+1 ;
    if countIter == nStop
        % did not reach saturation level 1-dFE with allowed no. of steps
        mu = nan ;
        RS = nan ;
        Lthr = nan ;
        Jthr = nan ;
        return
    end
end

% downward sweep
countSat  = 0 ;
countIter = 0 ;
while ( countSat < nSat )
    if countIter==0
        levels(end+1)= db2mag(-step_db)*startLevel ;
    else
        levels(end+1)= db2mag(-step_db)*levels(end) ;
    end
    if levels(end)<0
        levels(end) = [] ;
        break
    end
    [FEs(:,end+1),varFEs(:,end+1),Lats(:,end+1),Jits(:,end+1),VS(:,end+1),FF(:,end+1)] = Library.computeFE(fA,maskerE+levels(end)*probeE,nrep,reptime,baselineSR,parameters{:}) ;
    if FEs(1,end) <= dFE
        countSat = countSat+1 ;
    end
    countIter = countIter+1 ;
    if countIter == nStop
        % did not reach saturation level dFE with allowed no. of steps
        mu = nan ;
        RS = nan ;
        Lthr = nan ;
        Jthr = nan ;
        return
    end
end

% sort levels and FEs
[levels,FEs,varFEs,Lats,Jits,VS,FF] = Library.sortAll(levels,'ascend',FEs,varFEs,Lats,Jits,VS,FF) ;

% fit FE curve
[~,pBest] = Library.fitFE(levels,FEs(1,:)) ;
mu        = pBest(1) ;
sigma     = pBest(2) ;

for iter = 1:4
    if ~isnan(mu)
        if iter==1
            % simulate FE curve between mu-5.0*sigma and mu+5.0*sigma
            range = mu + 5.0*sigma*[-1 1] ;
        elseif iter==2
            % simulate FE curve between mu-3.5*sigma and mu+3.5*sigma
            range = mu + 3.5*sigma*[-1 1] ;
        elseif iter>=3
            % simulate FE curve between mu-2.0*sigma and mu+2.0*sigma
            range = mu + 2.0*sigma*[-1 1] ;
        end
        levelsInRange = levels( (levels >= range(1)) & (levels <= range(2)) ) ;
        N = max(3,nLevels-length(levelsInRange)) ;
        levelTmp = nan(1,N) ;
        for i=1:N
            tmp = sort([range levelsInRange levelTmp]) ;
            [dlevel,ind] = nanmax(diff(tmp)) ;
            levelTmp(i)  = tmp(ind) + dlevel/2 ;
        end
        levels = [levels levelTmp] ;
        FEs    = [FEs    nan(nW,N)] ;
        varFEs = [varFEs nan(nW,N)] ;
        Lats   = [Lats   nan(nW,N)] ;
        Jits   = [Jits   nan(nW,N)] ;
        VS     = [VS     nan(nW,N)] ;
        FF     = [FF     nan(nW,N)] ;
        for i=1:N
            [FEs(:,end-N+i),varFEs(:,end-N+i),Lats(:,end-N+i),Jits(:,end-N+i),VS(:,end-N+i),FF(:,end-N+i)] = Library.computeFE(fA,maskerE+levelTmp(i)*probeE,nrep,reptime,baselineSR,parameters{:}) ;
        end

        % sort levels and FEs again
        [levels,FEs,varFEs,Lats,Jits,VS,FF] = Library.sortAll(levels,'ascend',FEs,varFEs,Lats,Jits,VS,FF) ;

        % fit FE curve again
        [~,pBest] = Library.fitFE(levels,FEs(1,:)) ;
        mu        = pBest(1) ;
        sigma     = pBest(2) ;
    end
end
RS = sigma / mu ;

% finally, calculate statistics at threshold
if nargout>2
    if ~isnan(mu)
        [FEs(:,end+1),varFEs(:,end+1),Lats(:,end+1),Jits(:,end+1),VS(:,end+1),FF(:,end+1)] = Library.computeFE(fA,maskerE+mu*probeE,nrepStat,reptime,baselineSR,parameters{:}) ;
        Lthr = Lats(:,end) ;
        Jthr = Jits(:,end) ;

        [levels,FEs,varFEs,Lats,Jits,VS,FF] = Library.sortAll([levels mu],'ascend',FEs,varFEs,Lats,Jits,VS,FF) ;
    else
        Lthr = nan(nW,1) ;
        Jthr = nan(nW,1) ;
    end
end

end