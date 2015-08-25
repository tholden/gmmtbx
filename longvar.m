% LONGVAR: Estimates the long run covariance of the sample moment condition  
%
% SYNTAX: [S, bandw] = longvar(pm, center, method, bandw)
%
% INPUT
% pm     : A vector with the values of the moment condition.
% center : A dummy variable that controls the form of the data that is 
%          used in the covariance matrix calculation. Set this input 
%          equal to: 
%          0 --> if each individual element of the moments' vector is already centered                
%          1 --> if you want to center the moments            
% method : There are three methods available for calculating the
%          covariance matrix. Set this input equal to:
%          'HACC_B' --> for HACC with Bartlett kernel
%          'HACC_P' --> for HACC with Parzen kernel1 
%          'SerUnc' --> for Serially uncorrelated 
% bandw  : The bandwidth used in the HACC estimates. It must be a
%          non-negative integer (or zero, in the case of serially 
%          uncorrelated moments). 
%          IF THIS INPUT IS empty ([]), the program will autoamtically 
%          calculate the optimal bandwidth, using Newey and Wests's Method 
%          of Bandwidth Selection
%
% OUTPUT
% S      : The estimated long run covariance matrix.
% bandw  : Returns the value of the bandwidth used (this can help us check 
%          the value of the "optimal" bandwidth that the program has generated)  

function [S, bandw] = longvar(pm, center, method, bandw)
[T,~] = size(pm);

% Center data if necessary 
if center == 1;
     m = mean(pm);
    pm = pm - meshgrid(m,ones(T,1));
end;

% Calculate the optimum bandwidth if it is not given by the user
if strcmpi( method, 'hacc_b' ) && isempty(bandw)
    bandw = optbandw(pm, 'Bartlett');
elseif strcmpi( method, 'hacc_p' ) && isempty(bandw)
    bandw = optbandw(pm, 'Parzen');
end

% Calculate covariance matrix
switch lower(method)
    case('serunc')
        S = (1/T)*(pm'*pm);
    case('hacc_b')
        D = kernelest(T, bandw, 'Bartlett');
        tmp = chol( D ) * pm;
        % S = (1/T)*pm'*D*pm;
        S = (1/T)*(tmp' * tmp);
    case('hacc_p')
        D = kernelest(T, bandw, 'Parzen');
        tmp = chol( D ) * pm;
        % S = (1/T)*pm'*D*pm;
        S = (1/T)*(tmp' * tmp);
    otherwise
       error('Unknown method for moments'' variance.')
end