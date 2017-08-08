function [Ht, anum, aden] = iirxform(Ho,fun,varargin)
%IIRXFORM IIR Transformations

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% This should be private
error(message('signal:dfilt:singleton:iirxform:FreqXfmNotSupported', get( Ho, 'FilterStructure' )));


% [EOF]
