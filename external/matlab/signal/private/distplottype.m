function plotType = distplottype(n, esttype)
%DISTPLOTTYPE Helper function for plotting power and psd estimates of
%   distortion functions

%   Copyright 2013 The MathWorks, Inc.
if n>0
  plotType = 'none';
elseif strcmp(esttype,'psd')
  plotType = 'psd';
else
  plotType = 'power';
end
