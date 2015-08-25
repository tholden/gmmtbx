% RMSC: Computes the Relevance Moment Selection Criterion 
%
% SYNTAX: rmsc_stat = rmsc(p_vec, detVar_vec, q_vec, tau_vec)
%
% INPUT
% p_vec : A vector with the number of estimated parameters in each model you study
% detVar_vec  : A vector with the values of the determinant of the estimates' variance.
%               If you want to compare k different models(k>=1), detV must be a kx1 vector 
%               where its entries are determinant of the estimate's  variance for each model.
% qvec : A vector with the number of moments for each of the models you
%         want to compare.
% tau_vec  : A vector with the value of tau, that will be used for the penalty term of each model.
%            If the moments are martingale differences, set tau = T^0.5, otherwise set it equal to (T/b)^0.5
%            where T is the number of observations and b is the bandwidth used for the kernel estimation of 
%            the HAC covariance matrix
%
% OUTPUT
% rmsc_stat : A kx1 vector with the values of the RMSC statistic. The order
%             of the elements is the same as the order of Jvec, qvec, etc. 
%             The 1st entry is the RMSC of model_1, the 2nd is the RMSC of model_2, etc


function rmsc_stat = rmsc(p_vec, detVar_vec, q_vec, tau_vec)
% Error check
if nargin<4
    error('All the inputs are required for the estimation')
end

ratio = log(tau_vec)./tau_vec;
rmsc_stat = log(detVar_vec) + (q_vec - p_vec).*ratio;