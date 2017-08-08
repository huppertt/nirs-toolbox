function eml_freqz_plot(s,b,a,h,w,isTF,varargin)
% A supporting function for an Embedded MATLAB Library Function

% Copyright 2009-2012 The MathWorks, Inc.

if isTF
  phi = phasez(b,a,varargin{:});
else
  % Input is SOS matrix
  phi = phasez(b,varargin{:});
end
  
data(:,:,1) = h;
data(:,:,2) = phi;
% Turn off "obsolete" warning before calling freqzplot.  When freqzplot gets
% obsoleted, it will move into signal/private, so eml_freqz_plot function
% will continue to work.
ws = warning('off','signal:freqzplot:obsoleteFunction');
freqzplot(data,w,s,'magphase');
warning(ws);
