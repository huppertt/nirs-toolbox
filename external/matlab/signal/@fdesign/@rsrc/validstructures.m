function vs = validstructures(this, varargin)
%VALIDSTRUCTURES Return the valid structure for the design method.

%   Copyright 2005-2011 The MathWorks, Inc.

% Parse the SystemObject input
[varargin, sysObjFlag] = parsesysobj(this, 'validstructures', varargin{:}); 

if this.InterpolationFactor > this.DecimationFactor
    vs_str = {'firsrc'};
else
    vs_str = {'firsrc'};
end

if isempty(varargin)
  d = designmethods(this);
  dIntersect = d;
  if sysObjFlag
      dSysObj = designmethods(this,'SystemObject',sysObjFlag);
      dIntersect = intersect(d,dSysObj);
      hf = feval(getdesignobj(getcurrentspecs(this),dIntersect{1}));
      sysObjSupportedStructs = getsysobjsupportedstructs(hf);
      vs_str = intersect(vs_str,sysObjSupportedStructs); 
  end  
else
  d = varargin{1};
  dIntersect = d;
  if sysObjFlag
      dSysObj = designmethods(this,'SystemObject',sysObjFlag);
      dIntersect = intersect(d,dSysObj);
      if isempty(dIntersect)
        vs = cell(1,0);
        return;
      end
      hf = feval(getdesignobj(getcurrentspecs(this),dIntersect{1}));
      sysObjSupportedStructs = getsysobjsupportedstructs(hf);
      vs_str = intersect(vs_str,sysObjSupportedStructs); 
  end  
end


% If we are not given a design method, return a structure.
if isempty(varargin)    
    vs = struct;
    for indx = 1:length(dIntersect)
        vs.(dIntersect{indx}) = vs_str;
    end
else
    if isempty(dIntersect)
      vs = cell(1,0);
    else
      if ~any(strcmpi(d,varargin{1}))
        error(message('signal:fdesign:rsrc:validstructures:InvalidMethod', varargin{1}));
      end
      vs = vs_str;
    end
end

% [EOF]
