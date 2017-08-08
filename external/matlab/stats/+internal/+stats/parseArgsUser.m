function [varargout]=parseArgsUser(pnames,dflts,delim,varargin)
%parseArgsUser Process parameter name/value pairs for statistics functions
%   [A,B,...,USER,SETFLAG] =
%       parseArgsUser(PNAMES,DFLTS,DELIM,'NAME1',VAL1,'NAME2',VAL2,...)
%   accepts a cell array PNAMES of valid parameter names, a cell array
%   DFLTS of default values for the parameters named in PNAMES, the
%   name of a "user parameter delimiter", and additional parameter name/value
%   pairs.  Returns:
%
%        parameter values A,B,... in the same order as the names in PNAMES
%        a single cell array USER of all input args that follow the name
%           specified by DELIM.
%        a structure indicating which parameters were set explicitly
%
%   Outputs corresponding to entries in PNAMES that are not set explicitly
%   by the name/value pairs are set to the corresponding value from DFLTS.
%   USER is set to [] if the name in  DELIM does not appear in the name/value
%   pairs.  Unrecognized name/value pairs are an error.
%
%   [A,B,...,USER,SETFLAG,EXTRA] =
%       parseArgsUser(PNAMES,DFLTS,DELIM,'NAME1',VAL1,'NAME2',VAL2,...)
%   has nargout equal to length(PNAMES)+3. All unrecognized name/value
%   pairs are returned in a single cell array EXTRA.
%
%   Example:
%       pnames = {'color' 'linestyle', 'linewidth'}
%       dflts  = {    'r'         '_'          '1'}
%       delim = 'userargs'
%       varargin = {'linew' 2 'nonesuch' [1 2 3] 'linestyle' ':' ...
%          'userargs' 'pretty much' 'anything here' [1 2 3] {1 2 3}}
%       [c,ls,lw,ua] = internal.stats.parseArgsUser(pnames,dflts,varargin{:})  % error
%       [c,ls,lw,ua,us,ur] = internal.stats.parseArgsUser(pnames,dflts,varargin{:}) % ok
%       % c is 'r', ls is ':', and lw is 2
%       % ua is {'pretty much' 'anything here' [1 2 3] {1 2 3}}
%       % ur is {'nonesuch' [1 2 3]}

%   Copyright 2011-2014 The MathWorks, Inc. 


% We always create (nparams+2) outputs:
%    nparams varargs for values corresponding to names in pnames
%    one more for all user values
%    one more for a vector of flags indicating which parameters were set explicitly
% If they ask for one more (nargout == nparams+3), it's for unrecognized
% names/values

% Initialize some variables
nparams = length(pnames);
nargs = length(varargin);
allpnames = [pnames {delim}]; % (nparams+1)th element is the user args delimiter
userargs = {};

% Search for the DELIM parameter only
for j=1:2:nargs-1
    pname = varargin{j};
    if ischar(pname)
        tf = strncmpi(pname,allpnames,length(pname));
        if sum(tf)>1
            tf = strcmpi(pname,allpnames);
        end
        if tf(end)
            userargs = varargin(j+1:end); % user args
            varargin(j:end) = [];         % remove before further parsing
            break
        end
    end
end

% Process all remaining arguments, with or without "extra" arguments
if nargout>=nparams+3
    [varargout{1:nparams},setflag,extra] = internal.stats.parseArgs(pnames,dflts,varargin{:});
    varargout{nparams+3} = extra;
else
    [varargout{1:nparams},setflag] = internal.stats.parseArgs(pnames,dflts,varargin{:});
end
varargout{nparams+1} = userargs;
varargout{nparams+2} = setflag;
