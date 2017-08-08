function [Vp,Yp] = nlminputnoutput(this,L,M)
%NLMINPUTNOUTPUT   

%   Copyright 2006-2011 The MathWorks, Inc.

% Where are going to filter so work with copy
Hc = copy(this); 
Hc.PersistentMemory = false; % Make sure to avoid initial conditions

[vp,Vp] = nlmgeninput(this,M,L);

% Rescale input to be in input quantizer range
[vp,Vp] = nlmrescaleinput(get_filterquantizer(this),vp,Vp);

v = [vp; vp]; % Avoid transients?
 
v = quantizeinput(get_filterquantizer(this),v);

if isprop(this,'FromSysObjFlag') && this.FromSysObjFlag
  if isfi(v) 
    v.fimath = [];
  end
  y = step(this.ContainedSysObj, v);
  release(this.ContainedSysObj);
else
  y = filter(Hc,v);      % Applying v to the system under test
end

y = double(y);         % Cast

y = y(M+1:2*M,:);        % Selecting and transforming the last period
  
Yp = fft(y);

% [EOF]
