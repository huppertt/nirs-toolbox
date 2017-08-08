%QREPORT Quantization report.
%   R = QREPORT(Hd) returns the logging report R stored in the filter
%   object Hd. The logging capability is integrated to the 'filter' method
%   with the 'Arithmetic' property set to 'fixed'. It is triggered when the
%   'Logging' FI preference is 'on'. The stored report corresponds to the
%   last simulation. It is overwritten each time the filter command is
%   executed.
%
%   QREPORT provides users with a way to instrument fixed-point filters and
%   gives insight on how the filter is responding to a given stimulus. The
%   report object R contains a structure-specific list of internal signals
%   of the filter. Each signal contains:
%   - the minimum and maximum values that were recorded during the last
%     simulation (values are logged before quantization),
%   - the range of the signal,
%   - the number of overflows.
%
%   Notice that this method requires Fixed-Point Designer.
%
%   % EXAMPLE #1: Quantization report of a Direct-Form FIR filter.
%   fipref('LoggingMode', 'on');
%   Hd = design(fdesign.lowpass, 'equiripple');
%   Hd. arithmetic = 'fixed';
%   y = filter(Hd, rand(100,1));
%   R = qreport(Hd)
%
%   % EXAMPLE #2: Quantization report of a Direct-Form II, Second-Order Sections IIR filter.
%   fipref('LoggingMode', 'on');
%   Hd = design(fdesign.lowpass, 'ellip');
%   Hd. arithmetic = 'fixed';
%   y = filter(Hd, rand(100,1));
%   R = qreport(Hd)
%
%   See also DFILT/FUNCTIONS

%   Author(s): V. Pellissier
%   Copyright 2005-2011 The MathWorks, Inc.



% [EOF]
