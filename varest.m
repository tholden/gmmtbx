% VAREST: Computes the covariance matrix, std. errors, and confidence intervals of the GMM estimators.
%
% SYNTAX: [VAR, SD, CI] = varest(D, S, theta, T)
%
% INPUT
% D     : The gradient of the sample moment conditions.
% S     : The long-run covariance matrix estimator.
% theta : The GMM estimates.
% T     : The number os observations.
%
% OUTPUT
% VAR : The variance-covariance matrix of the estimators (size pxp, where p is the number of estimators).
% SD  : The 95% standard error [ = sqrt(Var/T) ] of the estimates (px1 vector).
% CI  : The 95% confidence intervals for the estimators(written in matrix form. An px2 matrix has the entries of these intervals).

function [VAR,SD,CI]=varest(D,S,theta,T);
    
      W = inv(S);
 invDWD = inv(D'*W*D);
    VAR = invDWD*D'*W*S*W*D*invDWD;
    SDA = diag(VAR);
     SD = (1/sqrt(T))*sqrt(SDA);

[nrth,ncth] = size(theta);
         CI = zeros(nrth,2);

for i=1:nrth
   CI(i,1) = theta(i,1)-(1.96*sqrt(VAR(i,i))/sqrt(T));
   CI(i,2) = theta(i,1)+(1.96*sqrt(VAR(i,i))/sqrt(T));
end