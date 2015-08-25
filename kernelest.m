% KERNELEST: Computes a matrix where its entries are the Bartlett or Parzen weights 
%           
% SYNTAX: K = kernelest(T, b)
%
% INPUT
% T          : The desired dimension for the weighting matrix to be computed
% b          : The bandwidth of the kernel. It must be a non-negative integer
% kernel_type: The function currently supports 2 different types of kernels.
%              Set to:
%              'Bartlett' for Bartlett kernel
%              'Parzen'   for Parzen kernel
%              
% OUTPUT
% K: A symmetric matrix where its entries are the Bartlett or Parzen weights 
%   (based on the method that the user selected)

function K = kernelest(T,b,kernel_type);
% ERROR CHECK
if b<0, error('The bandwidth must be non-negative.');end
integ = int16(b);
if integ~=b
    warning('The bandwidth you entered is not an integer; rounding bandwidth to nearest integer.');
    b=round(b);
end
if b > T
    error('The bandwidth (i.e. b) must be less than T, the dimansion of M.'); 
end

% The kernel estimation starts here
W2 = zeros(T-b-1,1);

a = (1:b)/(b+1);

if b == 0
    K = eye(T);
else
    switch lower(kernel_type)
        case{'bartlett'}
            W1 = [1;(1-a)'];
        case{'parzen'}
            W1 = [1 (1-6*a.^2 + 6*a.^3).*(a<=0.5) + 2*(1-a).^3.*(a>0.5)]';
        otherwise
            error('Incorrect kernel_type choice.')
    end
    W  = [W1;W2];
    K = (toeplitz(W))';
end