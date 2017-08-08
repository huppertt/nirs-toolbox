function measure(this)
%MEASURE  Measure the frequency response characteristics of a filter.
%   MEASURE(H) displays measuremens of various quantities from the
%   frequency response of the filter. Values that are measured can include
%   the actual passband ripple, the minimum stopband attenuation, the
%   frequency point at which the filter's gain is 3 dB below the nominal
%   passband gain, etc. In order for MEASURE to work, the filter H needs to
%   be designed using FDESIGN. See examples below.
%
%   M = MEASURE(H) returns the measurements in such a way that they can be
%   queried programmatically. For example, to query the 3 dB point, one
%   would type M.F3dB. Type GET(M) to see the full list of properties that
%   can be queried. Note that different filter responses will generate
%   different measurements.
%
%   For designs that do not specify some of the frequency constraints, it
%   may not be possible to determine corresponding magnitude measurements.
%   In these cases, a frequency value can be passed in to MEASURE in order
%   to determine the magnitude measurements that corresponds to such value.
%   See example #3 below for an illustration of this.
%
%   % Example #1: Measure an equiripple FIR design of a bandpass filter
%   f = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',0.3,0.4,0.6,0.7,80,1,75);
%   H = design(f,'equiripple');
%   measure(H)
%
%   % Example #2: Determine the minimum stopband attenuation of a
%   % multistage decimation filter
%   f = fdesign.decimator(4,'lowpass','Fp,Fst,Ap,Ast',0.2,0.22,1,75);
%   H = design(f,'multistage');
%   M = measure(H); M.Astop
%
%   % Example #3: 
%   f = fdesign.lowpass('N,F3dB,Ast',8,0.5,80);
%   h = design(f,'cheby2');
%   measure(h,'Fpass',0.4)
%
%   See also DFILT/INFO, DFILT/COST.

%   Author(s): R. Losada
%   Copyright 2006-2010 The MathWorks, Inc.



% [EOF]
