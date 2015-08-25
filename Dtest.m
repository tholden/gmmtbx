function Dt = dtest(num_restr, unrestr_theta, restr_theta, unrestr_S, popmom, data, varargin);
% DTEST: Computes the D test for linear/nonlinear restrictions on the GMM estimates.
%
% SYNTAX
% Dt = dtest(num_restr, unrestr_theta, restr_theta, unrestr_S, 'popmom', data, varargin);
% 
% INPUTS                 
% num_restr       : The number of restrictins you are testing.
% unrestr_theta   : The unrestricted GMM estimates
% restr_theta     : The restricted GMM estimates  
% unrestr_S       : The moments' long-run covariance matrix, evaluated at the
%                   final unrestricted estimates
% popmom          : An M-file that calculates the moment conditions and their 
%                   gradient. The file must be of the form 
%                   [mom, gradmom] = popmom(theta,data,varargin)
%                   where mom is the moments and gradmom their gradient.                 
% data            : An Txm matrix with the dataset that the moments are based on 
%                  (m variables with T data-points each). 		
% [varargin]      : Additional parameters passed to the popmom function
%                   For example, in the case of a GIV estimation, 
%                   you can pass the vector of the instruments.
%
% OUTPUT
% Dt         : An output structure with the following fields
% Dt.value   : The value of the D test 
% Dt.prob    : The probability value of D  


% Get the weighting matrix
invS = inv(unrestr_S);

% Get the # of observations and the # of restrictions 
num_obs = size(data, 1);

% D test
unrestr_J = gobj(unrestr_theta, popmom, data, invS, varargin{:});
restr_J   = gobj(restr_theta, popmom, data, invS, varargin{:});
Dt.value   = (restr_J - unrestr_J);
Dt.prob  = 1-chi2cdf(Dt.value, num_restr);