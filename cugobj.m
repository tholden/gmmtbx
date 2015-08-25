% CUGOBJ: Computes the CU GMM objective function
%
% SYNTAX: [obj, gradobj]=gobj(theta, popmom, data, W, varargin); 
%
% INPUT
% theta     : A vector with the estimated parameters.
% popmom    : An m-file that calculates the moment conditions and its gradient.
%             The file must be of the type: [mom, gradmom] = popmom(data, theta, varargin)
% data      : A matrix containing the dataset used for estimation.
% W         : The weighting matrix used for calculating the objective function
% [varargin]: Additional parameters passed to the popmom function
%
% OUTPUT
% obj    : The value of the CU GMM objective function.
%
% SETABLE OPTIONS (use OPTSET):
%   center : A dummy variable, taking the values 0 or 1
%            Set to 0: The variance of the moments is calculated using the 
%                           uncentered moment conditions               
%            Set to 1: The variance of the moments is calculated using the 
%                          centered moment conditions      
%                          ##The default value is 0## 
%   method : Covariance matrix estimation method. Set this option equal to            
%           'HACC_B', for HAC with Bartlett kernel
%           'HACC_P', for HAC with Parzen kernel1 
%           'SerUnc', for Serially uncorrelated 
%           ##The default method is 'SerUnc'## 
%   bandw  : The bandwidth used in THE HAC estimation. The bandwidth must 
%           be a non-negative integer.    
%           ## If the bandwidth is not given by the user, and a HACC estimator 
%           has been selected by the user, the "optimal" bandwidth is automatically 
%           calculated using Newey and Wests's Method of Bandwidth Selection
% 
% For example, if you want to use a Bartlett type kernel with a bandwidth of 2 use the commands:
% optset('cugmmest', 'method', 'HACC_B'); 
% optset('cugmmest', 'bandw', 2); 


function obj = cugobj(theta, popmom, data, center, method, bandw, varargin);

center    = optget('cugmmest','center',0);
method    = optget('cugmmest','method','SerUnc');
bandw     = optget('cugmmest','bandw',0);

[pmc, dpmc] = feval(popmom, theta, data, varargin{:});   % Get the moments and their derivative
S = longvar(pmc, center, method, bandw);                      % Get the variance of the moments 
W = inv(S);
obs = size(pmc,1);
g = sum(pmc)';
obj = (1/obs)*g'*W*g;