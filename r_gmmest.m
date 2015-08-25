% R_GMMEST: GMM estimation under linear or nonlinear resrictions 
% ---------------------------------------------------------------
% SYNTAX: [theta_final, J_test, probJ, S_final, final_moments,... 
%          final_moments_grad, var_theta, std_theta,conf_inter] = 
%          r_gmmest(options, data, 'popmom', stval, Wstart, varargin);
% 
%---------------------------------INPUTS------------------------------
% OPTIONS  : Options structure for fmincon (use OPTIMSET)
%            Set the options structure for fmincon, the constrained
%            minimization function used by MATLAB. To see the list of
%            possible options, enter optimset('fminunc'). Give an empty
%            matrix,[], if you want to use the default values.
%            IMPORTANT: The default value of MaxSQPIter is Inf;                                                                       
%            this may result to an infinite loop during the                               
%            estimation. This "bug" is discussed in Matlab's 
%            "Technical Solution" webpage, Solution #: 35698.
% DATA      : An mxn matrix with the data that the estimation will use 
%            (n variables with m data-points each). 		
% POPMOM    : An M-file that calculates the moment conditions and their 
%            gradient. The file must be of the form 
%            [m, dm]=popmom(theta,data,varargin)
%            where m is the moments and dm their gradient.                 
% STVAL     : A vector with an initial guess for the parameters. 	
% WSTART    : The weighting matrix for the first step estimation. 	
% [varargin]: Additional parameters passed to the popmom function
%
% ---------------------------------OUTPUT---------------------------------
% THETA_FINAL       : A vector with the final GMM estimates.
% J_TEST            : The J-test of overidentifying restrictions.
% PROBJ             : The prob of the J-test.
% S_FINAL           : The weighting matrix computed at the final estimates. 
% FINAL_MOMENTS     : A matrix with the values of the moment conditions,
%                     evaluated at the final GMM estimates
% FINAL_MOMENTS_GRAD: The gradient of the moment conditions, evaluated at 
%                     the final GMM estimates
% BANDW             : The value of bandwidth used to calculate S_final
% VAR_THETA         : The variance-covariance matrix of the estimators.
% STD_THETA         : The standard deviations of the estimators.
% CONF_INTER        : The 95% confidence interval of each estimator.
%
% -------------------------------------------------------------------------
% SETABLE OPTIONS FOR R_GMMEST (use OPTSET)
% center : A dummy variable, taking the values 0 or 1
%          Set to 0: The variance of the moments is calculated 
%                    using the uncentered moment conditions               
%          Set to 1: Center the moments and then calculate 
%                    their variance      
%          ##The default value is 0## 
% method : Covariance matrix estimation method. Set this option equal to            
%          'HACC_B', for HACC with Bartlett kernel
%          'HACC_P', for HACC with Parzen kernel1 
%          'SerUnc', for Serially uncorrelated 
%          ##The default method is 'SerUnc'## 
% bandw  : The desired bandwidth, used in HACC estimation. 
%           The bandwidth must be a non-negative integer.    
%           ## If the bandwidth is not given by the user, the 
%           optimal bandwidth is automatically calculated, using 
%           Newey and Wests's Method of Bandwidth Selection
% itergmm : Maximum number for the iterated GMM estimator.
%           ##(The default is 100).
% tol     : Tolernace criterion for the termination of the GMM iterations
% LB      : Lower Bound for theta. This vector must have the same
%           dimensions with theta. Set it to [] if there are no bounds. If
%           there is bound only in one of the parameters (e.g. theta2>10,
%           theta1 unbounded, set it to [-Inf, 10].
% UB      : Upper Bound for theta. It uses the same rules as the Lowr Bound
%           above.
% NONLCON : A function that calculates a linear/nonlinear, equality
%           constraint and its derivative. In order for this function to be
%           compatible with fmincon, it must return 4 arguments. However,
%           in our case we need only two (value of the constraint and its
%           derivative). To ensure compatibility, write the function in the
%           following way:
%           [c1 value c2 derivative] = myfunction(theta, varargin)
                                                                                                                

function [theta_final, J_test, probJ, S_final, final_moments,... 
          final_moments_grad,bandw, var_theta, std_theta,conf_inter] = r_gmmest(options, data, popmom, startval, Wstart, varargin);
  
% SIZE OF DATASET & STARTING VALUES
[dr,dc]     = size(data);
[stvr,stvc] = size(startval);

% ERROR CHECK
if nargin<5, error('The first four inputs (data, popmom, stval, W) must be provided by the user');end
if stvc~=1, error('The starting values must be a column vector');end
if stvc>dc, error('The system is un-identified. You must supply at least as many conditions as parameters.');end

% OPTIONS STRUCTURE FOR GMM (DEFAULT VALUES) 
center    = optget('r_gmmest','center',0);
method    = optget('r_gmmest','method','SerUnc');
bandw     = optget('r_gmmest','bandw',[]);
itergmm   = optget('r_gmmest','itergmm',50);
tol       = optget('r_gmmest','tol',1e-006);
LB        = optget('r_gmmest','LB',[]);
UB        = optget('r_gmmest','UB',[]);
NONLCON   = optget('r_gmmest','NONLCON',[]); 
                                                                                                                                                                                                         
% First step estimator
theta = fmincon('gobj', startval, [], [], [], [], LB, UB, NONLCON, options, popmom, data, Wstart, varargin{:});

% The iterative, restricted, GMM estimation starts here
for i=2:itergmm
    pmc = feval(popmom, theta,data, varargin{:});
    S = longvar(pmc, center, method, bandw);
    [thetanew,fval] = fmincon('gobj', theta,[],[],[],[], LB, UB, NONLCON, options, popmom, data, inv(S), varargin{:});
    if norm(abs(theta-thetanew))<tol
        result = sprintf('The algorithm converged to a solution. The optimal estimator was achieved in iteration %2.0f .', i);
        disp(result);
        break
    end
    theta = thetanew;
end
    

% OUTPUT
if itergmm == 1
    theta_final = theta;
    pmc = feval(popmom, theta,data, varargin{:});
    S   = longvar(pmc, center, method, bandw);
else
    theta_final = thetanew;
end

[pmc,dpmc]  = feval(popmom, theta_final,data, varargin{:}); 
W_final     = inv(S);
[S_final, bandw] = longvar(pmc, center, method, bandw); 
S_final     = S_final;
J_test      = gobj(theta_final, popmom, data, inv(S_final), varargin{:});
[pr,pc]     = size(pmc); 
df          = pc - stvr;
probJ       = 1-chi2cdf(J_test, df);
[VAR,SD,CI] = varest(dpmc, S_final, theta_final, dr);
var_theta   = VAR;
std_theta   = SD;
conf_inter  = CI;
final_moments       = pmc;
final_moments_grad  = dpmc;

% USER NOTIFICATIONS
if isempty(optget('r_gmmest', 'bandw')) & lower(optget('r_gmmest', 'method'))~='serunc'
    message = sprintf('The optimum bandwidth, has been set to %4.0f', bandw);
    disp(message);
end