function [rnd, ofmode] = conv2hdlroundoverflow(this)
%CONV2HDLROUNDOVERFLOW Converts the Round & Overflow Mode to HDL Round and
%HDL Overflow modes

%   OUT = CONV2HDLROUNDOVERFLOW(ARGS) <long description>

%   Copyright 2007 The MathWorks, Inc.

if strcmpi(this.arithmetic, 'fixed')
    rnd = get(this, 'RoundMode');
    ofmode = strcmpi(this.overflowmode, 'saturate');
 else % double
    rnd = 'floor';
    ofmode = false;
end

% [EOF]
