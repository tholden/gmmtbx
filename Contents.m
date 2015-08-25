% GMM Toolbox for Matlab
% Last Revision: June 08, 2004
%
% MAIN FUNCTIONS
% --------------
%    GMMEST: Main procedure for the unrestricted GMM estimation 
%  R_GMMEST: GMM estimation under linear or nonlinear resrictions 
%      GOBJ: Computes the GMM objective function and its gradient
%    CUGOBJ: Computes the CU GMM objective function
%   LONGVAR: Estimates the long run covariance of the sample moment 
%            condition  
% KERNELEST: Computes a matrix where its entries are the Bartlett or 
%            Parzen weights
%  OPTBANDW: Newey and Wests's Method of optimum Bandwidth Selection. 
%    VAREST: Computes the covariance matrix, std. errors, and confidence 
%            intervals of the estimates
% 
% HYPOTHESIS TESTING 
% ------------------
% WTEST    : Computes the Wald test for linear/nonlinear restrictions on 
%            the GMM estimates
% LMTEST   : Computes the LM test for linear/nonlinear restrictions on 
%            the GMM estimates
% DTEST    : Computes the D test for linear/nonlinear restrictions on the 
%            GMM estimates.
% SSTESTS  : Structural stability tests with known break-point  
% SSTESTSU : Sequence of structural stability tests with unknown break-point
% UBSS_STAT: Sup-, Av-, and Exp- functionals of the unknown break-point 
%            structural stability tests 
% MSC      : Computes the Moment Selection Criterion 
% RMSC     : Computes the Relevance Moment Selection Criterion 
% 
% GRAPHICAL USER INTERFACE
% ------------------------
% GMMGUI     : GUI object for interactive estimation of restricted, unrestricted, 
%              and Continuous Updated GMM models.
% GMMTESTING : GUI object to perfrom Wald, LM, amd D tests
% MOMSEL     : GUI object to calculate the Moment Selection Criterion
% RELMOMSEL  : GUI object to calculate the Relevance Moment Selection Criterion
% STRSTAB    : GUI object to calculate structural stability tests
%
% HELPER FUNCTIONS
% ----------------
% OPTGET   : Utility to get previously set function default values (by P. Fackler and M. Miranda) 
% OPTSET   : Utility to set function options (by P. Fackler and M. Miranda)
%
%  Copyright (c) 2004   Kostas Kyriakoulis.