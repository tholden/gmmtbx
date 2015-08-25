% CUGMMEST: Main procedure for the Continuous Updated GMM estimation 
%
% SYNTAX 
%  [theta_final, S_final, J_test, probJ, bandw, var_theta, std_theta, conf_inter] = 
%		         cugmmest(options, data, popmom, startval, varargin);
% 
% INPUTS - OUTPUT
%
% The function uses a syntax similar to the one of gmmest. Please see the help file 
% of gmmest for more details concerning the inputs-output 
% -------------------------------------------------------------------------

function [theta_final, S_final, J_test, probJ, bandw, var_theta, std_theta, conf_inter] = ...
    cugmmest(options, data, popmom, startval, varargin);

% SIZE OF DATASET & STARTING VALUES
[dr,dc]     = size(data);
[stvr,stvc] = size(startval);

% ERROR CHECK
if nargin<5, error('The first four inputs (data, popmom, stval, W) must be provided by the user');end
if stvc~=1, error('The starting values must be a column vector');end
if stvc>dc, error('The system is un-identified. You must supply at least as many conditions as parameters.');end

% OPTIONS STRUCTURE FOR CU-GMM (DEFAULT VALUES) 
center    = optget('cugmmest','center',0);
method    = optget('cugmmest','method','SerUnc');
bandw     = optget('cugmmest','bandw',0);
itergmm   = optget('cugmmest','itergmm',50);
tol       = optget('cugmmest','tol',1e-006);

% CU estimator
theta_final = fminunc('cugobj', startval, options, popmom, data, center, method, bandw, varargin{:}); 
disp('The CU-GMM estimates have been computed.' );
[pmc,dpmc] = feval(popmom, theta_final,data, varargin{:}); 
[S_final, bandw] = longvar(pmc, center, method, bandw); 
J_test = gobj(theta_final, popmom, data, inv(S_final), varargin{:});
[pr,pc] = size(pmc); 
df = pc - stvr;
probJ = 1-chi2cdf(J_test, df);
[VAR,SD,CI] = varest(dpmc, S_final, theta_final, dr);
var_theta = VAR;
std_theta = SD;
conf_inter = CI;
final_moments = pmc;
final_moments_grad = dpmc;

% USER NOTIFICATIONS
if stvr == pc
    disp('The model is just identified; the probability value of the J test is meaningless.')
end

if isempty(optget('cugmmest', 'bandw')) & lower(optget('cugmmest', 'method'))~='serunc'
    message = sprintf('The optimum bandwidth, has been set to %4.0f', bandw);
    disp(message);
end