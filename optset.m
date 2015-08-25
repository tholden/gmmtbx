% OPTSET Utility to set function options 
% USAGE
%   optset(funcname,optname,optvalue)
% INPUTS
%   funcname : name of function
%   optname  : name of option
%   optval   : option value
%
% If optname='defaults' the current setting of the options will be
%    cleared. The next time the function is called, the default
%    options will be restored.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function optvalue = optset(funcname,optname,optvalue)

global options_

funcname = lower( funcname );
optname = lower( optname );
if ~isstruct( options_ )
    options_ = struct;
end
if ~isfield( options_, funcname )
    options_.( funcname ) = struct;
end
FunctionOptions = options_.( funcname );
FunctionOptions.( optname ) = optvalue;
if strcmp( optname, 'defaults' )
    FunctionOptions = rmfield( FunctionOptions, optname );
end
options_.( funcname ) = FunctionOptions;
end