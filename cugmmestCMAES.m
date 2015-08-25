% CUGMMEST: Main procedure for the Continuous Updated GMM estimation 
%
% SYNTAX 
%  [theta_final, S_final, J_test, probJ, bandw, var_theta, std_theta, conf_inter] = 
%		         cugmmest(options, data, popmom, startval, varargin);
% 
% 
% INPUTS
%   options   : Options structure for fmincon (use OPTIMSET)
%               Set the options structure for fminunc, the unconstrained
%               minimization function used by MATLAB. To see the list of
%               possible options, enter optimset('fminunc').          
%   data      : A Txm matrix with the dataset that the moments are based on 
%               (m variables with T observations each). 		
%   popmom    : An Matlab function that calculates the moment conditions and their 
%               gradient. The function must be of the form 
%               [mom, gradmom] = popmom(theta,data,varargin)
%               where mom is the moments and gradmom their gradient.                 
%   stval     : A vector with an initial guess for the parameters. 	
%   Wstart    : The weighting matrix for the first step estimation. 	
%   [varargin]: Additional parameters passed to the popmom function
%               For example, in the case of a GIV estimation, 
%               you can pass the vector of the instruments.
%
% OUTPUT
%   theta_final    : A vector with the final GMM estimates.
%   J_test         : The J-test of overidentifying restrictions.
%   probJ          : The probability of the J-test.
%   S_final          The long-run covariance matrix of the moments, 
%                    evaluated at the final GMM estimates.
%   moments        : A matrix with the values of the moment conditions,
%                    evaluated at the final GMM estimates
%   moments_grad   : The gradient of the moment conditions, evaluated at 
%                    the final GMM estimates
%          bandw   : The bandwidth used to calculate S_final
%   var_theta      : The variance-covariance matrix of the estimates.
%   std_theta      : The standard error of each estimates, defined as 
%                    [Var(estimate)/T]^0.5.
%   conf_inter     : The 95% confidence interval of each estimate.
% -------------------------------------------------------------------------
% 
% SETABLE OPTIONS (use OPTSET):
%   center : A dummy variable, taking the values 0 or 1
%            Set to 0: The variance of the moments is calculated using the 
%                    uncentered moment conditions               
%            Set to 1: The variance of the moments is calculated using the 
%                    centered moment conditions      
%            ##The default value is 0## 
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
%   tol     : Tolerance criterion for the iteration procedure.
%           ##(The default is 1e-006).

function [theta_final, S_final, J_test, probJ, bandw, var_theta, std_theta, conf_inter] = ...
    cugmmestCMAES(options, data, popmom, startval, sigma, nonlcon, varargin)

% SIZE OF DATASET & STARTING VALUES
[dr,dc]     = size(data);
[stvr,stvc] = size(startval);

% ERROR CHECK
if nargin<4, error('The first four inputs (options, data, popmom, startval) must be provided by the user');end
if stvc~=1, error('The starting values must be a column vector');end
if stvc>dc, error('The system is un-identified. You must supply at least as many conditions as parameters.');end
if nargin<6
    nonlcon = [];
end
% OPTIONS STRUCTURE FOR CU-GMM (DEFAULT VALUES) 
center    = optget('gmmest','center',0);
method    = optget('gmmest','method','SerUnc');
bandw     = optget('gmmest','bandw',0);

% CU estimator
theta_final = cmaes(@(XV) parallel_wrapper( @(X) cugobj(X, popmom, data, center, method, bandw, varargin{:}) + 1e12 * sum( max( 0, nonlcon(X, popmom, data, center, method, bandw, varargin{:}) ) ), XV ), startval, sigma, options ); 
disp('The CU-GMM estimates have been computed.' );
if nargout > 1
	[pmc,dpmc] = popmom( theta_final,data, varargin{:}); 
	[S_final, bandw] = longvar(pmc, center, method, bandw); 
	J_test = gobj(theta_final, popmom, data, inv(S_final), varargin{:});
	[~,pc] = size(pmc); 
	df = pc - stvr;
	probJ = 1-chi2cdf(J_test, df);
	[VAR,SD,CI] = varest(dpmc, S_final, theta_final, dr);
	var_theta = VAR;
	std_theta = SD;
	conf_inter = CI;

	% USER NOTIFICATIONS
	if stvr == pc
		disp('The model is just identified; the probability value of the J test is meaningless.')
	end

	if isempty(optget('cugmmest', 'bandw')) && ~strcmpi( optget('cugmmest', 'method'),'serunc' )
		message = sprintf('The optimum bandwidth, has been set to %4.0f', bandw);
		disp(message);
	end
end