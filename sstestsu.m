% SSTESTSU: Sequence of structural stability tests with unknown break-point  
%
% SYNTAX: [ss_statvec, thetas] = sstestsu(options, data, popmom, thetafull, Sf, Dfull, pi_l, pi_u, varargin)
%
% INPUT
% options    :  Options structure defining the fminunc (and/or fmincon) optional parameters. 
%               Type optimset for more help 
% data       :  A Txm matrix with the dataset that the moments are based on 
%               (m variables with T observations each). 		
% popmom     :  An Matlab function that calculates the moment conditions and their 
%               gradient. The function must be of the form 
%               [mom, gradmom] = popmom(theta,data,varargin)
%               where mom is the moments and gradmom their gradient.                 
% theta_full :  A vector with the full sample GMM estimates.
% S_full     :  The long-run covariance matrix of the moments, 
%               evaluated at the full sample GMM estimates.
% D_full     :  The gradient of the moment conditions, evaluated at 
%               the full sample GMM estimates
% pi_l       :  Lower bound of the testing interval. This must be a positive number, 
%               between 0 and 1 (it must be bounded away from 0 and 1).
% pi_u       :  Upper bound of the testing interval. This must be a positive number, 
%               between 0 and 1 (it must be bounded away from 0 and 1).
% varargin   :  Additional arguments passed to the popmom function
%
% OUTPUT
% ss_statvec : A (pi_u*T - pi_l*T + 1)x5 matrix with the values of the statistics at each point in the 
%              [pi_l*T , pi_u*T] interval, and the value of the breakpoint. The columns are order as
%              the follows:
%               Column 1 |  Column 2 | Column 3 | Column 4 | Column 5
%              Wald-test |  LM-test  |  D-test  |  O-test  | break_point 
% thetas     : A 2P x (pi_u*T - pi_l*T + 1) matrix with the value of the
%              estimates before and after the break. Each column
%              corresponds to a different breakpoint. The first p rows are
%              the estimates before the break, and the last p are the
%              estimates after it.
%
% EXAMPLE on pi_l and pi_u
% Suppose that you have a dataset of 200 observations and you wish to test
% for a break between the observations 70 and 130. Then 
% pi_l = 70/200 = 0.35 and
% pi_u = 130/20 = 0.85

function [ss_statvec, thetas] = sstestsu(options, data, popmom, thetafull, Sf, Dfull, pi_l, pi_u, varargin)

[T,col] = size(data);
TL = round(pi_l*T);
TU = round(pi_u*T);
TS = TU-TL+1;

for bp = TL:TU;
   [theta1,theta2,Wss,LMss,Dss,O]=...
       sstests(options, data, popmom, thetafull, Sf, Dfull, bp, varargin{:});
i = bp-TL+1; 
ss_statvec(i,:) = [Wss.value LMss.value Dss.value O.value bp];
thet1(:,i) = theta1;
thet2(:,i) = theta2;
end
thetas = [thet1;thet2];