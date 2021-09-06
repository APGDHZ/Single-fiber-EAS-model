function randnum = myrand(sz,type)
% Calls rand(sz(1),sz(2)) to produce a sz(1) * sz(2) matrix of uniformly 
% distributed random numbers, but additionally allows to specify a fixed or
% variable rng seed.
% For type=1: no rng seed/variable random numbers; 
% for type=0: rng seed used/fixed random numbers.
% 
% The type argument is optional, default=1.
%

%---- Check input arguments ---------- %

if ( (nargin < 1) || (nargin > 2) )
	error('Requires One to Two input arguments.')
end

if (numel(sz) ~= 2) || any(~isnumeric(sz)) || any(~isreal(sz)) || any(~isfinite(sz))
	error('sz must be an 1*2 array of finite real scalars.')
end

if ((sz(1) <= 0) || (sz(2) <= 0))
	error('Dimensions of the return vector must be positive.')
end

if (nargin > 1)
	if ( (numel(type) ~= 1) || ~isnumeric(type) || ~isreal(type) || ~isfinite(type) )
        error('type must be a finite real scalar.')
    end
else
    type = 1; % default is variable seed
end
	
%---- create random numbers ---------- %
if (type==0) %fixed random numbers
    rng(123);
end
randnum = rand(sz);

end

