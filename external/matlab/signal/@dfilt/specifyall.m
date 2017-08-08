function specifyall(this,flag)
%SPECIFYALL  Fully specify fixed-point filter settings.
%   SPECIFYALL(Hd) sets the properties of the object so that you can specify
%   all the settings in the filter Hd directly. All autoscale options that
%   are true are set to false and all "modes" (such as 'OutputMode',
%   'ProductMode', 'AccumMode') are set to 'SpecifyPrecision'.
%
%   SPECIFYALL is intended as a shortcut to changing all these properties
%   one-by-one.
%
%   SPECIFYALL(Hd,false) is essentially the opposite of SPECIFYALL(Hd). All
%   autoscale options are set to true and "modes" are set to their default
%   values.
%
%   SPECIFYALL(Hd,true) is equivalent to SPECIFYALL(Hd).
%
%   % Example: Provide full control over the fixed-point settings of an FIR
%   % filter implemented with the direct-form structure.
%   f = fdesign.lowpass('N,Fp,Fst,Ast',12,0.4,0.5,20); 
%   Hd = design(f,'equiripple'); Hd.Arithmetic = 'fixed';
%   specifyall(Hd);
%
%   See also DFILT/REFFILTER, DFILT/DOUBLE, DFILT/COPY.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.



% [EOF]
