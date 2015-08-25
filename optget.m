% OPTGET Utility to get previously set function default values 
% USAGE
%   optvalue=optget(funcname,optname,optvalue);
% INPUTS
%   funcname : name of function
%   optname  : name of option
%   optval   : option value
% OUTPUT
%   optval   : the current value of the option
%
% If the named field is not already defined, it will be set to
% optvalue, but optvalue has no effect if the field has already 
% been set. Use OPTSET to change a previously set field.
%
% optget(funcname) returns the current values of the options structure.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function optvalue = optget(funcname,optname,optvalue)

global options_
funcname = lower( funcname );
optname = lower( optname );
if isstruct( options_ )
    if isfield( options_, funcname )
        FunctionOptions = options_.( funcname );
        if isfield( FunctionOptions, optname )
            optvalue = FunctionOptions.( optname );
        end
    end
end
optset(funcname,optname,optvalue);

end