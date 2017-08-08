%DESIGN   Design the filter from the specifications.
%   H = DESIGN(D)  Design the filter H from the specifications in D.  The
%   algorithm will be chosen from those available for the given
%   specification. 
%
%   DESIGN(...) launches FVTool to visualize the designed filter.
%
%   H = DESIGN(D, 'FIR') forces an FIR design.  This method will error if
%   there are no FIR designs available.
%
%   H = DESIGN(D, 'IIR') forces an IIR design.  This method will error if
%   there are no IIR designs available.
%
%   H = DESIGN(D, 'ALL') designs multiple filters using all of the
%   available design methods for the specifications.
%
%   H = DESIGN(D, 'ALLFIR') designs multiple filters using all of the
%   available FIR design methods for the specifications.
%
%   H = DESIGN(D, 'ALLIIR') designs multiple filters using all of the
%   available IIR design methods for the specifications.
%
%   H = DESIGN(D, METHOD) forces the design method specified by the
%   string METHOD to be used.  METHOD must be one of the strings
%   returned by <a href="matlab:help fdesign/designmethods">designmethods</a>. Use designmethods(D,'default') to determine 
%   which algorithm is used by default.
%
%   H = DESIGN(D, METHOD, PARAM1, VALUE1, PARAM2, VALUE2, etc.) specifies
%   design method specific options.  Use <a href="matlab:help fdesign/help">help(D, METHOD)</a> for more
%   information on optional inputs.
%
%   H = DESIGN(D, METHOD, OPTS) specifies design method specific options
%   using the structure OPTS. OPTS is usually obtained from <a href="matlab:help fdesign/designopts">designopts</a>,
%   modified, and then specified as an input to DESIGN. Use <a href="matlab:help fdesign/help">help(D, METHOD)</a>
%   for more information on optional inputs.
%
%   When the structure returned by the <a href="matlab:help fdesign/designopts">designopts</a> method contains a 'SystemObject'
%   field, then H = DESIGN(D, ..., 'SystemObject', true) will implement the
%   filter design using a filter System object H.
%
%   % Example #1 - Design the filter using the default (equiripple) method
%   %              and visualize the filter's response in FVTool.
%   d = fdesign.lowpass('Fp,Fst,Ap,Ast',.2, .22, 1, 60);
%   Hd = design(d)
%   fvtool(Hd)
%   designmethods(d,'default')
%
%   % Example #2 - Design all possible IIRs for the specifications given.
%   Fs = 48e3; % Sampling frequency is 48 kHz
%   d = fdesign.lowpass('Fp,Fst,Ap,Ast',10000,11000,0.5,80,Fs);
%   designmethods(d, 'iir')
%   Hd = design(d, 'alliir'); % Hd is a vector of IIR filters
%   fvtool(Hd)
%
%   % Example #3 - Design a Kaiser window highpass FIR filter.
%   d  = fdesign.highpass('Fst,Fp,Ast,Ap',0.35, 0.4,74,1);
%   Hd = design(d, 'kaiserwin');
%   fvtool(Hd)
%   
%   % Example #4 - Design a 50th order equiripple FIR filter with a
%   %              sloped stopband. Type help(d, 'equiripple') for more
%   %              info (*)
%   d = fdesign.lowpass('N,Fc,Ap,Ast',50, 0.4, 0.8, 80);
%   Hd = design(d, 'equiripple', 'StopbandShape', 'linear', 'StopbandDecay', 40);
%   fvtool(Hd)
%
%   % Example #5 - Design a lowpass filter using the interpolated FIR
%   %              method (*)
%   d = fdesign.decimator(4,'lowpass','Fp,Fst,Ap,Ast',.2, .22, 1 ,60);
%   opts = designopts(d,'ifir');
%   opts.JointOptimization = true;
%   Hd = design(d,'ifir',opts); % Hd has two FIR filters cascaded
%   fvtool(Hd)
%
%   % Example #6 - Design a lowpass Butterworth filter and implement it
%   %               using a filter System object (*)
%   d = fdesign.lowpass('Fp,Fst,Ap,Ast',.2, .22, 1, 60);
%   Hd = design(d,'butter','SystemObject',true)
%
%   %(*) DSP System Toolbox required
%
%   See also FDESIGN, FDESIGN/DESIGNMETHODS, FDESIGN/DESIGNOPTS, FDESIGN/HELP.

%   Copyright 2005-2011 The MathWorks, Inc.

% [EOF]
