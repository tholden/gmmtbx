% UBSS_STAT: Sup-, Av-, and Exp- functionals of the unknown break-point str. stability tests  
%
% SYNTAX: [Sup, bpt, Av, Exp] = ubss_stat(stat)
%
% INPUT
% stat    :  This matrix must have a similar format with ss_statvec, that appears in the function sstestsu.
%            Type 'help sstestsu' for more details
%            IMPORTANT: Each row of the matrix represents the data-point on which the struct. stability
%            test were performed. The last column of the stat matrix must be this datapoint
% OUTPUT
% Sup     : The supremum of the sequence of the tests 
% bpt     : The datapoint were the supremum occured
% Av      : The average of the sequence of the tests
% Exp     : The exponential average of the sequence of the tests
% 
% EXAMPLE 
% For illustration purposes, suppose that you have a dataset of 20
% observations and you performed a sequence of str. stability tests on all
% observations from the 14th till the 18th. You have calculated the
% following sequence of Wald and D statistics:
%
% Wald = [5 9 2 8 6]';        
%   LM = [3 7 4 12 5]';
% 
% To calculate the Sup-, AV-, and Exp- functionals of these statistice, you
% first need to create a vector with the data-points on which the tests
% were performed. Since we tested all observations from the 14th till the
% 18th, this vector will be
% 
% DataPoints = (14:18)';
%
% In order to use this function, create a matrix that holds the tests and
% the data-points, i.e.
%
% demo_tests = [Wald LM DataPoints];
%
% and use this function by simply entering the following command
%
% [demoSup demoBp demoAv demoExp] = ubss_stat(demo_tests)
%
% What you should get is:
%
% demoSup =
%     9    12
%
% demoBp =
%    15    17
%
% demoAv =
%    6.0000    6.2000
%
% demoExp =
%    3.5813    4.5231
%
% In the above results, the first column corresponds to the Wald test,
% whereas the second to the LM (since this is the order we input them
% originally). For example: 
% The supremum of the LM test is 12 and occurs on observation 17
% The average of the Wald is 6.0000.
% etc.


function [Sup, bpt, Av, Exp] = ubss_stat(stat);
[Sup ind] = max(stat(:,1:end-1));
bpt = stat(ind, end)';
Av = mean(stat(:,1:end-1));
Exp_a = mean(exp(stat(:,1:end-1).*0.5));
Exp = log(Exp_a);