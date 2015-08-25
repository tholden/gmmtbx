% MSC: Computes the Moment Selection Criterion 
%
% SYNTAX: msc_stat = msc(method, b, p, Jvec, qvec, Tvec)
%
% INPUT
% method: The type of the bonus term. This input can be 
%         'BIC'    : The Bayesian Inf. Criterion is used for the penalty term
%         'HanQuin': The Hannan-Quinn criterion is used for the penalty term
% b     : A finite constant, greater than 2, used by the HanQuin method.
%         Give an empty matrix, [], if you use the BIC method.
% p     : The number of parameters
% Jvec  : A vector with the values of the J-test for each of the models you
%         want to compare. If you want to compare k different models
%         (k>=1), the Jvec must be a kx1 vector where its entries are the
%         values of the J-test for each model.
% qvec : A vector with the number of moments for each of the models you
%         want to compare.
% Tvec  : A vector with the number of observations for each of the models you
%         want to compare.
%
% OUTPUT
% msc_stat : A kx1 vector with the values of the MSC statistic. The order
%            of the elements is the same as the order of Jvec, qvec, etc. 
%            The 1st entry is the MSC of model_1, the 2nd is the MSC of model_2, etc

function msc_stat = msc(method, b, p, Jvec, qvec, Tvec) 

% Error check
if nargin < 6
    error('All the inputs are required for the estimation')
end

switch lower(method)
    case 'bic'
        msc_stat = Jvec - (qvec - p).*log(Tvec);
    case 'hanquin'
        if isempty(b) | b < 2+1e-8
            error('The Hannan and Quinn method requires a numeric value for b, greater than 2');
        end
        msc_stat = Jvec - (qvec - p).*b.*log( log(Tvec) );
    otherwise
        error('Unknown penelty method')
end