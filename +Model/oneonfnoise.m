function gnoise=oneonfnoise(n,NoiseAlpha,varargin)

% gnoise=oneonfnoise(n,NoiseAlpha);
% creates n samples of 1/f^NoiseAlpha Gaussian noise .
% The noise is scaled accordingly to produce coect dynamic range of the
% neuron, for any given alpha (tested between 0 - 2).

if mod(n,2)==0
    f = [1:n/2]';
    p = 2*pi*rand(n/2,1);
    yc = f.^(-NoiseAlpha/2).*exp(i.*p.*f);
    s = [0+0i; yc; conj(yc(end-1:-1:1))];
else
    f = [1:(n-1)/2]';
    p = 2*pi*rand((n-1)/2,1);
    yc = f.^(-NoiseAlpha/2).*exp(i.*p.*f);
    s = [0+0i; yc; conj(flipud(yc))];
end
x = ifft(s);
gnoise = real(x);
gnoise = gnoise/std(gnoise);

switch nargin
    case 2
        % Scale the noise for this particular neuron model
        a =       5.073; %  (4.705, 5.442)
        b =       -2.61; %  (-2.91, -2.31)
        c =      0.3265; %  (-0.03942, 0.6925)
        d =    -0.00834; %  (-0.6009, 0.5842)
        sigma = a*exp(b*NoiseAlpha) + c*exp(d*NoiseAlpha);
        gnoise = gnoise*sigma;
    case 3
        gnoise = gnoise;
end

