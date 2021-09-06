function [stimE,nPulses] = constructPulseTrain(singlePulse,pulseRate,trainDuration)

    fs = 1e6 ;

    Npulse  = length(singlePulse) ;
    Ntrain  = round(trainDuration*fs) ;
    stimE   = zeros(1,Ntrain) ;
    nPulses = 0 ;
        
    if Npulse > Ntrain
        warning('trainDuration must be longer than singlePulse duration.')
        return
    end

    iPulses = max(1, round( 0 : (fs/pulseRate) : Ntrain  )) ;
    for i=iPulses
        if i+Npulse <= Ntrain
            stimE(i:(i+Npulse-1)) = singlePulse ;
            nPulses = nPulses + 1 ;
        end
    end
    
end