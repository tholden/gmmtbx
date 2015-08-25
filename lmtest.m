function LM = lmtest(num_restr, restr_theta, unrestr_S, popmom, data, varargin);
% LMTEST: Computes the LM test for linear/nonlinear restrictions on the GMM estimates.
%
% SYNTAX
% LM = lmtest(num_restr, restr_theta, unrestr_S, 'popmom', data, varargin);
% 
%
% INPUTS
%                   
% num_restr       : The number of restrictins you are testing.
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
% LM         : An output structure with the following fields
% LM.value   : The value of the LM test. 
% LM.prob    : The probability value of LM 


% Get the weighting matrix
invS = inv(unrestr_S);
% Get the # of observations and the # of restrictions 
num_obs = size(data, 1);
                                              
% Evaluate the moments, and their gradient, at the restricted estimates.
[gtilde,Gtilde] = feval(popmom, restr_theta,data, varargin{:}); 
% LM test
gtilde = (sum(gtilde)/num_obs)';
restr_invGSG = inv(Gtilde'*invS*Gtilde);
gSG = gtilde'*invS*Gtilde;
LM.value = num_obs*gSG*restr_invGSG*gSG';
LM.prob = 1-chi2cdf(LM.value, num_restr);