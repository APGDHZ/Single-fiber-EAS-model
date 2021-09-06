function stimA = constructNoise(noiseDuration,rampDuration,level_db)

    fs = 1e5 ;
    
    N     = round(noiseDuration*fs) ;
    Nramp = round(rampDuration*fs) ;
    
    stimA = randn(1,N) ;
    stimA(          1:Nramp) = linspace(0,1,Nramp) .* stimA(          1:Nramp) ;
    stimA(end-Nramp+1:end  ) = linspace(1,0,Nramp) .* stimA(end-Nramp+1:end  ) ;
    
    if nargin>2
        level_Pa = 2e-5 * db2mag(level_db) ;
        stimA    = level_Pa * stimA ;
    end
    
end