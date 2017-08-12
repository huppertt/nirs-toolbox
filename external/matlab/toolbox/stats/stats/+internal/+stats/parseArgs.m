function [varargout]=parseArgs(pnames,dflts,varargin)
%parseArgs Process parameter name/value pairs for statistics functions
%   [A,B,...] = parseArgs(PNAMES,DFLTS,'NAME1',VAL1,'NAME2',VAL2,...)
%   In typical use there are N output values, where PNAMES is a cell array
%   of N valid parameter names, and DFLTS is a cell array of N default
%   values for these parameters. The remaining arguments are parameter
%   name/value pairs that were passed into the caller. The N outputs
%   [A,B,...] are assigned in the same order as the names in PNAMES.
%   Outputs corresponding to entries in PNAMES that are not specified
%   in the name/value pairs are set to the corresponding value from DFLTS. 
%   Unrecognized name/value pairs are an error.
%
%   Each element of PNAMES is either a string with a parameter name or a
%   cell array of strings with possible names for the same parameter.
%
%   [A,B,...,SETFLAG] = parseArgs(...), where SETFLAG is the N+1 output
%   argument, also returns a structure with a field for each parameter
%   name. The value of the field indicates whether that parameter was
%   specified in the name/value pairs (true) or taken from the defaults
%   (false).
%
%   [A,B,...,SETFLAG,EXTRA] = parseArgs(...), where EXTRA is the N+2 output
%   argument, accepts parameter names that are not listed in PNAMES. These
%   are returned in the output EXTRA as a cell array.
%
%   Example:
%       pnames = {'color' 'linestyle', 'linewidth'}
%       dflts  = {    'r'         '_'          '1'}
%       varargin = {'linew' 2 'linestyle' ':'}
%       [c,ls,lw] = internal.stats.parseArgs(pnames,dflts,varargin{:})
%       % On return, c='r', ls=':', lw=2
%
%       [c,ls,lw,sf] = internal.stats.parseArgs(pnames,dflts,varargin{:})
%       % On return, sf = [false true true]
%
%       varargin = {'linew' 2 'linestyle' ':' 'special' 99}
%       [c,ls,lw,sf,ex] = internal.stats.parseArgs(pnames,dflts,varargin{:})
%       % On return, ex = {'special' 99}

%   Copyright 2010-2014 The MathWorks, Inc.

[varargout{1:nargout}] = statslib.internal.parseArgs(pnames,dflts,varargin{:});


