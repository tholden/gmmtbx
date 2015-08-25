% SSTESTS: Structural stability tests with known break-point  
%
% SYNTAX: [theta1, theta2, Wss, LMss, Dss, O, O1, O2] = ...
%    sstests(options, data, popmom, thetaf, Sfull, Df, bp, varargin)
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
% theta_f    :  A vector with the full sample GMM estimates.
% Sfull      :  The long-run covariance matrix of the moments, 
%               evaluated at the full sample GMM estimates.
% Df         :  The gradient of the moment conditions, evaluated at 
%               the full sample GMM estimates
% bp         :  The point of the dataset you wish to test for a break.
%               This musy be a positive integer between 0 and T (it must be bounded away from 0 and T).
% varargin   :  Additional arguments passed to the popmom function
%
% OUTPUT
% theta1     : The estimates before the break 
% theta2     : The estimates after the break 
% Wss        : A structure variable with two fields
%              Wss.value = The Wald statistic
%              Wss.prob  = The prob. of the Wald statistic
% LMss       : A structure variable with two fields
%              LMss.value = The LM statistic
%              LMss.prob  = The prob. of the LM statistic
% Dss        : A structure variable with two fields
%              Dss.value = The D statistic
%              Dss.prob  = The prob. of the D statistic
% O          : A structure variable with two fields
%              Oss.value = The O statistic
%              Oss.prob  = The prob. of the O statistic
% O1         : A structure variable with two fields
%              O1ss.value = The O1 statistic
%              O1ss.prob  = The prob. of the O1 statistic
% O2         : A structure variable with two fields
%              O2ss.value = The O2 statistic
%              O2ss.prob  = The prob. of the O2 statistic

function [theta1, theta2, Wss, LMss, Dss, O, O1, O2] = ...
    sstests(options, data, popmom, thetaf, Sfull, Df, bp, varargin)

center   = optget('gmmest','center');
method   = optget('gmmest','method');
bandw    = optget('gmmest','bandw');
itergmm  = optget('gmmest','itergmm');
tol      = optget('gmmest','tol');

% Check the number of inputs
if nargin<7
    error('The first 7 inputs must be provided')
end

% Get the observation where the break occurs
T   = size(data,1);
pit = bp/T;

% Get the number of observations
p = size(thetaf, 1);

% Create the subsets, before and after the break 
data1  = data(1:bp,:);
data2  = data(bp+1:T,:);
% If instruments are provided (i.e. the varargin is not empty), we have to
% create their subsets are well
if nargin == 8
    INSTR_USED  = varargin{1};
    INSTR_USED1 = INSTR_USED(1:bp,:);
    INSTR_USED2 = INSTR_USED(bp+1:T,:);
    % Get the first sub-sample estimates, and the variance of the moments   
    theta1  = gmmest(options, data1, popmom, thetaf, Sfull, INSTR_USED1);
    [m1u,D1u] = feval(popmom, theta1, data1, INSTR_USED1); 
    [m1r,D1r] = feval(popmom, thetaf, data1, INSTR_USED1);
    S1u = longvar(m1u, center, method, bandw);
    invS1u = inv(S1u);
    Q1u = gobj(theta1, popmom, data1, invS1u, INSTR_USED1);
    Q1r = gobj(thetaf, popmom, data1, invS1u, INSTR_USED1);
    % Get the second sub-sample estimates, and the variance of the moments  
    theta2  = gmmest(options, data2, popmom, thetaf, Sfull, INSTR_USED2);
    [m2u,D2u] = feval(popmom, theta2, data2, INSTR_USED2);
    [m2r,D2r] = feval(popmom, thetaf, data2, INSTR_USED2);
    S2u = longvar(m2u, center, method, bandw);
    invS2u = inv(S2u);
    Q2u = gobj(theta2, popmom, data2, invS2u, INSTR_USED2);
    Q2r = gobj(thetaf, popmom, data2, invS2u, INSTR_USED2);
else
    theta1  = gmmest(options, data1, popmom, thetaf, Sfull);
    [m1u,D1u] = feval(popmom, theta1, data1); 
    [m1r,D1r] = feval(popmom, thetaf, data1); 
    S1u = longvar(m1u, center, method, bandw);
    invS1u = inv(S1u);
    Q1u = gobj(theta1, popmom, data1, invS1u);
    Q1r = gobj(thetaf, popmom, data1, invS1u);
    % Get the second sub-sample estimates, and the variance of the moments  
    theta2  = gmmest(options, data2, popmom, thetaf, Sfull);
    [m2u,D2u] = feval(popmom, theta2, data2); 
    [m2r,D2r] = feval(popmom, thetaf, data2); 
    S2u = longvar(m2u, center, method, bandw);
    invS2u = invS2u;
    Q2u = gobj(theta2, popmom, data2, invS2u);
    Q2r = gobj(thetaf, popmom, data2, invS2u);
end

% Calculate the Wald test
VW1 = D1u'*invS1u*D1u;
VW2 = D2u'*invS2u*D2u;
VW  = (1/pit)*inv(VW1)+(1/(1-pit))*inv(VW2);
c   = theta1-theta2;
Wss.value = T*c'*inv(VW)*c;
Wss.prob = 1-chi2cdf(Wss.value,p);

% Calculate the LM test                
g1  = pit*(sum(m1r)/size(m1r,1))';
invS  = inv(Sfull);
LMv = Df'*invS*g1;
sc  = T/(pit*(1-pit));
LMss.value = sc*LMv'*inv(Df'*invS*Df)*LMv;
LMss.prob = 1-chi2cdf(LMss.value,p);

% Calculate the D test    
Dss.value   = (Q1r-Q1u)+(Q2r-Q2u);
Dss.prob = 1 - chi2cdf(Dss.value,p);

% Calculate O test      
q = size(m1r, 2);
O1.value = Q1u;
O2.value = Q2u;
O.value = Q1u + Q2u; 
dfo = 2*(q-p);
O.prob = 1-chi2cdf(O.value,dfo);
O1.prob = 1-chi2cdf(O1.value,dfo/2);
O2.prob = 1-chi2cdf(O2.value,dfo/2);