% GOBJ: Computes the GMM objective function and its gradient       
%
% SYNTAX: [obj, gradobj]=gobj(theta, popmom, data, W, varargin); 
%
% INPUT
% theta     : A vector with the estimated parameters.
% popmom    : An m-file that calculates the moment conditions and its gradient.
%             The file must be of the type: [mom, gradmom] = popmom(data, theta, varargin)
% data      : A matrix containing the dataset used for estimation.
% W         : The weighting matrix used for calculating the objective function
% [varargin]: Additional parameters passed to the popmom function
%
% OUTPUT
% obj    : The value of the GMM objective function.
% gradobj: The gradient of the objective function.

function obj = gobj(theta, popmom, data, W, varargin)
pmc = popmom( theta, data, varargin{:});
obs = size(pmc,1);
g = sum(pmc)';
obj = (1/obs)*g'*W*g;
% if nargout>1
%     gradobj = 2*g'*W*dpmc;
% end
        