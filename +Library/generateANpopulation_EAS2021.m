function [sponts,tabss,trels,cfs,Cpers,Ccens]=generateANpopulation_EAS2021(nsponts,cfrange)

% parameters:
tabsmax = 1.5 * 461e-6 ;
tabsmin = 1.5 * 139e-6 ;
trelmax = 894e-6 ;
trelmin = 131e-6 ;

% for ffGn with tau=1e-6:
muPer  =   -6.1514 ;  % mu (peripheral neuron)
sigPer =    0.1947 ;  % sigma (peripheral neuron)
bmPer  = -164.0e-9 ;  % b/m (peripheral neuron)
muCen  =   -5.7547 ;  % mu (central neuron)
sigCen =    0.2010 ;  % sigma (central neuron)
bmCen  =  -32.7e-9 ;  % b/m (central neuron)

if nargin>1
    cfmin = cfrange(1) ;
    cfmax = cfrange(2) ;
end
if nargin<3
    cfmin =   125 ; % default value
    cfmax = 40000 ; % default value
end

% generate spont rates:
sponts_low  =  0.1 +  0.1 * randn(1,nsponts(1)) ;
sponts_med  =  4.0 +  4.0 * randn(1,nsponts(2)) ;
sponts_high = 70.0 + 30.0 * randn(1,nsponts(3)) ;
sponts_low  = truncate_rand(sponts_low , 1e-3,   0.2);
sponts_med  = truncate_rand(sponts_med ,  0.2,  18.0);
sponts_high = truncate_rand(sponts_high, 18.0, 180.0);
sponts = [sponts_low sponts_med sponts_high] ;

% generate refractory time constants:
refrand = rand(1,sum(nsponts)) ;
tabss   = (tabsmax - tabsmin)*refrand + tabsmin;
trels   = (trelmax - trelmin)*refrand + trelmin;

% generate CFs:
cfs     = 2.^( log2(cfmin) + rand(1,sum(nsponts))*log2(cfmax/cfmin) ) ;

% generate membrane capacitances:
Cper_rand = randn(1,sum(nsponts));
Ccen_rand = corr_randn(Cper_rand, 0.50); % correlation between Cper & Ccen = 0.50
Cpers     = 10.^( muPer + sigPer * Cper_rand ) - bmPer ;
Ccens     = 10.^( muCen + sigCen * Ccen_rand ) - bmCen ;
Cpers     = truncate_rand(Cpers, 10^(muPer-2*sigPer)-bmPer, 10^(muPer+2*sigPer)-bmPer) ;
Ccens     = truncate_rand(Ccens, 10^(muCen-2*sigCen)-bmCen, 10^(muCen+2*sigCen)-bmCen) ;

end

function Y = corr_randn(X,r)
    % For a normally distributed random array x,
    % generates an array of normally distributed random numbers y 
    % of the same size and correlation corr(x,y)=r  (0<=r<=1)
    % For a matrix x, each column is treated as a vector.
    Y = r*X + sqrt(1-r^2)*randn(size(X)) ;
end

function x = truncate_rand(x,xmin,xmax) 
    x = min(max(x,xmin),xmax) ;
end