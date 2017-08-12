function [selected,locs] = getParamVal(inputs,kwds,argname,multok)
%getParamVal Select from a finite set of choices
%   [PVAL,PLOC] = getParamVal(PVAL,KWDS,PNAME) process input parameters
%   that can take values from a finite set. On input, PVAL is a candidate
%   parameter value (string), KWDS is a cell array of strings defining all
%   legal parameter values, and PNAME is the name of the parameter (used
%   only if needed in an error message). getParamVal tries to match PVAL
%   to the KWDS values, ignoring case and allowing unambigous abbreviations.
%   On output, PVAL is the matched value and PLOC is its index in KWDS.
%
%   [PVAL,PLOC] = getParamVal(PVAL,KWDS,PNAME,true) allows for multiple
%   selections. Here PVAL can be a cell array on input and is a cell
%   array on output.
%
%   Example:
%       onoff = getParamVal(pval,{'on' 'off'},'Display')


%   Copyright 2010-2013 The MathWorks, Inc.

if nargin<4
    multok = false;
end

[selected,locs] = statslib.internal.getParamVal(inputs,kwds,argname,multok);

