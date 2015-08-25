% OPTBANDW: Newey and Wests's Method of optimum Bandwidth Selection. 
%
% SYNTAX: bandwidth = optbandw(moments, kernel_type)
% 
% INPUTS
% moments     : A Txq matrix with the value of the moment conditions 
%               (q moment conditions with T observations each).
%
% kernel_type : The type of kernel, used by the HAC estimates. Set to:
%               'Bartlett', if the Bartlett kernel is used
%               'Parzen', if the Parzen kernel is used
%               #### If a "method" is not supplied by the user, the default is 'Bartlett'.
%
% OUTPUT
% bandwidth   : An integer scalar that denotes the optimum bandwidth.    

function bandwidth = optbandw(moments, kernel_type)
[T,q] = size(moments);
% START WITH ERROR CHECK
if nargin<1, error('Insufficient number of inputs.'), end;
if T==1, error('You need more than on observation. Check the size of the moments.'), end;
if nargin<2, kernel_type = 'bartlett', disp('The default method (Bartlett) is used'),end;

% The main procedure starts here
b = moments*ones(q,1);

switch lower(kernel_type)
    case{'bartlett'}
        nu = 1;
        n = fix((T^(1/9))^2);
        cgamma = 1.1447;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';      % NW bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       
        j = (1:n)';                            % NW bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));  
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1));% NW bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
    case{'parzen'}
        nu  = 2;
        n   = fix((T^(1/25))^4);
        cgamma=2.6614;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';            % NW bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       % |
        j = (1:n)';                                  % |----> NW bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));        % |
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1)); % NW bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
        otherwise
        error('The kernel type must be "Bartlett", or "Parzen"');
end