%AUTOSCALE  Dynamic range scaling.
%   AUTOSCALE(Hd,X) provides dynamic range scaling for each node of the
%   filter Hd. This method runs the signal X through the filter in
%   floating-point and uses the maximum and minimum data obtained from that
%   simulation to set fraction lengths such that the simulation range is
%   covered and the precision is maximized. Word lengths are not changed.
%
%   HNEW = AUTOSCALE(Hd,X) If an output is requested a new filter is
%   generated with the scaled fraction lengths and the original filter is
%   not changed.
%
%   % EXAMPLE: Autoscale a bandpass IIR elliptic filter
%   Hd = design(fdesign.bandpass, 'ellip');
%   Hd = convert(Hd, 'latticearma');
%   Hd.Arithmetic = 'fixed';
%   x = rand(100,10);
%   Hd(2) = autoscale(Hd,x);
%   hfvt = fvtool(Hd,'Analysis','magestimate','Showreference','off');
%   legend(hfvt,'Before Autoscaling', 'After AutoScaling')
%
%   See also the <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\iirautoscaledemo.html'])">Fixed-Point Scaling of an Elliptic IIR Filter</a>
%   and the <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\iirfloat2fixeddemo.html'])">Floating-point to Fixed-Point Conversion of IIR Filters</a>
%   demos for more examples.
%
%   See also DFILT

%   Author(s): V. Pellissier
%   Copyright 2006-2010 The MathWorks, Inc.



% [EOF]
