function Wald = wtest(num_restr, unrestr_theta, unrestr_S, constrf, popmom, data, varargin)
% WTEST: Computes the Wald for linear/nonlinear restrictions on the GMM estimates.
%
% SYNTAX
% Wald = wtest(num_restr, unrestr_theta, unrestr_S, 'constrf', 'popmom', data, varargin);
% 
%
% INPUTS
%                   
% num_restr       : The number of restrictins you are testing.
% unrestr_theta   : The unrestricted GMM estimates
% unrestr_S       : The moments' long-run covariance matrix, evaluated at the
%                   final unrestricted estimates
% constrf         : A function that calculates the restrictions and their gradient. 
%                   The function must be of the form [c,ceq,GC,GCeq] = mycon(x),
%                   where 
%                      c = ...    Linear or Nonlinear inequalities at x
%                    ceq = ...    Linear or Nonlinear equalities at x
%                     GC = ...    Gradients of the inequalities
%                   GCeq = ...    Gradients of the equalities
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
% Wald          : An output structure with the following fields
% Wald.value    : The value of the LM test. 
% Wald.prob    : The probability value of LM 

% Get the weighting matrix
invS = inv(unrestr_S);
% Get the # of observations and the # of restrictions 
num_obs = size(data, 1);
                                              
% Evaluate the moments, and their gradient, at the unrestricted estimates.
[ghat,Ghat] = feval(popmom, unrestr_theta, data, varargin{:}); 
% Evaluate the constraint and their gradient at the unrestricted estimates
[omit1, r_unrestr, omit2, R_unrestr] = feval(constrf,unrestr_theta);
% Wald test
invGSG   = inv(Ghat'*invS*Ghat);
Wald.value  = num_obs*r_unrestr'*inv(R_unrestr*invGSG*R_unrestr')*r_unrestr; 
Wald.prob = 1-chi2cdf(Wald.value, num_restr);