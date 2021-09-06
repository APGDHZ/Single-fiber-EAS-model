function [fun, pBest] = fitFE(lvls, FEs)
    ind  = ~isnan(lvls) & ~isnan(FEs) ;
    if ( sum(ind) < 2 ) || all(FEs < 0.5) || all(FEs > 0.5)
        fun   = @(x) nan ;
        pBest = nan(2,1) ;
        return
    else
        lvls = lvls(ind) ;
        FEs  = FEs(ind) ;
    end
    ind = FEs >= .05 & FEs <= .95 ;
    if sum(ind) < 2
        p0 = [median(lvls)      ; median(lvls)     ] ;
    else
        p0 = [median(lvls(ind)) ; median(lvls(ind))] ;
    end
    fun   = @(pars) sseval(pars,lvls,FEs) ;
    [pBest,~,exitflag] = fminsearch(fun,p0) ;
    
    if exitflag ~= 1
        pBest = nan(2,1) ;
        fun   = @(x) nan ;
    else
        mu    = pBest(1) ;
        sigma = pBest(2) ;
        fun   = @(x) normcdf(x,mu,sigma) ;
    end
end

function sse = sseval(pars,lvls,FEs)
    mu    = pars(1);
    sigma = pars(2);
    sse   = sum((FEs - normcdf(lvls,mu,sigma)).^2);
end
