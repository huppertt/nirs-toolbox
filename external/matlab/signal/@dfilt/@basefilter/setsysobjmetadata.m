function setsysobjmetadata(this,Hs)
%SETSYSOBJMETADATA Set metadata of generated filter System object

%   Copyright 2011 The MathWorks, Inc.

fdtmp = getfdesign(this);
if ~isempty(fdtmp)
  fd = copy(fdtmp);
else
  fd = fdtmp;
end

fmtmp = getfmethod(this);
if ~isempty(fmtmp)   
  fm = copy(fmtmp);
  addsysobjdesignopt(fm); 
  fm.SystemObject = true;
else
  fm = fmtmp;
end

mtmp = get(this, 'privMeasurements');
if ~isempty(mtmp)
  m = copy(mtmp);
else
  m = mtmp;
end

fdesign      = fd;
fmethod      = fm;
measurements = m;
designmethod = get(this, 'privdesignmethod');

setMetaData(Hs,fdesign, fmethod, measurements, designmethod)
