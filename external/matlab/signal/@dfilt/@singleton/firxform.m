function Ht = firxform(Ho,fun,varargin)
%FIRXFORM FIR Transformations

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

error(message('signal:dfilt:singleton:firxform:FreqXfmNotSupported', get( Ho, 'FilterStructure' )));


% [EOF]
