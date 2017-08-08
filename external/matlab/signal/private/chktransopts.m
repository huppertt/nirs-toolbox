function varargout = chktransopts(x, opts, varargin)
%CHKTRANSOPTS check arguments for bilevel waveform transitions
%
%   out position   option                     default
%       1          'Tolerance'                2 
%       2          'MidPct' or 'PctRefLevels' 50 or [10 90]
%       3          'StateLevels'              statelevels(x)
%       4          'Polarity'                 'positive' 
%    4 and 5       'Aberration'               'postshoot' and 3
%
%   Duplicated/unrecognized options generate an error
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011-2013 The MathWorks, Inc.
  try
    p = varargin;
    [tol, p] = extractTol(p);
    varargout{1} = tol;
    
    if any(strcmpi(opts,'MidPct'))
      [varargout{2}, p] = extractMprl(tol,p);
    else
      [varargout{2}, p] = extractPrl(tol,p);
    end

    if any(strcmpi(p,'StateLevels'))
      [varargout{3}, p] = extractStateLevs(p);
    else
      varargout{3} = statelevels(x);
    end
    
    if any(strcmpi(opts,'Polarity'))
      [varargout{4}, p] = extractPol(p);
    end
    
    if any(strcmpi(opts,'Region'))
      [varargout{4}, p] = extractRgn(p);
      [varargout{5}, p] = extractFactor(p);
    end

    if ~isempty(p)
      if ischar(p{1})
        error(message('signal:chktransopts:UnexpectedOption',p{1}));
      else
        error(message('signal:chktransopts:UnexpectedInput'))
      end
    end
    
  catch ex
    throwAsCaller(ex);
  end
end

function [val, p] = findprop(p, prop, default)
  val = default;
  n = find(strcmpi(p, prop));
  if ~isempty(n)
    if numel(n) > 1
      error(message('signal:chktransopts:DuplicateProperty',prop))
    end
    if n>=numel(p)
      error(message('signal:chktransopts:MissingValue',prop))
    end
    val = p{n+1};
    p(n:n+1) = [];
  end
end

function [tol, p] = extractTol(p)
  [tol, p] = findprop(p,'tolerance',2);
  validateattributes(tol, ...
    {'double'},{'real','finite','positive','scalar','<',50},'','TOL');
end

function [mprl, p] = extractMprl(tol, p)
  [mprl, p] = findprop(p,'midpercentreferencelevel',50);
  if isequal(mprl,50)
    [mprl, p] = findprop(p,'midpct',50);
  end
  validateattributes(mprl, ...
    {'double'},{'real','finite','scalar'},'','L');
  if tol >= mprl
    error(message('signal:chktransopts:TolMustBeLessThanMidPctRefLev'));
  elseif 100 <= mprl+tol
    error(message('signal:chktransopts:TolAndMidPctRefLevTooBig'));
  end
end

function [prl, p] = extractPrl(tol, p)
  [prl, p] = findprop(p,'percentreferencelevels',[10 90]);
  if isequal(prl,[10 90])
    [prl, p] = findprop(p,'pctreflevels',[10 90]);
  end
  validateattributes(prl, ...
    {'double'},{'real','finite','size',[1 2]},'','PCTREFLEVELS');
  if prl(1)>=prl(2)
    error(message('signal:chktransopts:LowerMustBeLessThanUpperPctRefLev'));
  elseif tol >= prl(1)
    error(message('signal:chktransopts:TolMustBeLessThanLwrPctRefLev'));
  elseif 100 <= prl(2)+tol
    error(message('signal:chktransopts:TolAndUprPctRefLevTooBig'));
  end
end

function [stateLevs, p] = extractStateLevs(p)
  [stateLevs, p] = findprop(p,'statelevels',[0 2.3]);
  validateattributes(stateLevs, ...
    {'double'},{'real','finite','size',[1 2]},'','STATELEVELS');
  if stateLevs(1)>=stateLevs(2)
    error(message('signal:chktransopts:LowerMustBeLessThanUpperStateLev'));
  end
end

function [f, p] = extractFactor(p)
  [f, p] = findprop(p,'seekfactor',3);
  validateattributes(f, ...
    {'double'},{'real','finite','scalar','positive'},'','FACTOR');
end

function [polFlag, p] = extractPol(p)
  [polStr, p] = findprop(p,'polarity','positive');
  if strcmpi(polStr,'positive')
    polFlag = 1;
  elseif strcmpi(polStr,'negative')
    polFlag = -1;
  else
    error(message('signal:chktransopts:IllegalBooleanFlagValue','POLARITY','Positive','Negative'));
  end
end

function [rgn, p] = extractRgn(p)
  [rgn, p] = findprop(p,'region','postshoot');
  if ~strcmpi(rgn,'preshoot') && ~strcmpi(rgn,'postshoot')
    error(message('signal:chktransopts:IllegalBooleanFlagValue','REGION','Preshoot','Postshoot'));
  end
end