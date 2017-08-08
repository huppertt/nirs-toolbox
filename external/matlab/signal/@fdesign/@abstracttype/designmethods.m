function varargout = designmethods(this, varargin)
%DESIGNMETHODS Returns a cell of design methods.

%   Copyright 2004-2011 The MathWorks, Inc.

% Parse the SystemObject input
[varargin sysObjFlag] = parsesysobj(this, 'designmethods', varargin{:}); 

if any(strcmpi('default', varargin))
    d = defaultmethod(this);
    if any(strcmpi('full', varargin))
        switch lower(d) 
            case 'ellip'
                d = 'Elliptic';
            case 'butter'
                d = 'Butterworth';
            case 'cheby1'
                d = 'Chebyshev type I';
            case 'cheby2'
                d = 'Chebyshev type II';
            case 'kaiserwin'
                d = 'Kaiser window';
            case {'equiripple', 'window'}
                d = [upper(d(1)) d(2:end)];
            case 'firls'
                d = 'FIR least-squares';
            case 'fircls'
                d = 'FIR constrained least-squares';
            case 'ifir'
                d = 'Interpolated FIR';
            case 'iirlpnorm'
                d = 'IIR least p-norm';
            case 'freqsamp'
                d = 'Frequency sampling';
            case 'multistage'
                d = 'Multistage equiripple';
            case 'iirlinphase'
                d = 'IIR quasi-linear phase';
            case 'maxflat'
                d = 'Maximally flat';
            case 'ansis142'
                d = 'ANSI S1.42 weighting';
        end
    end
    d = {d};
    type = 'Default';
else
    [d, isfull, type] = thisdesignmethods(this, varargin{:}); 
    type = upper(type);
end

% If SystemObject has been passed as an input and its value is true, then
% remove structures that are not supported by System objects. 
if sysObjFlag
  idxVector = [];
  if isfull
    idx = strcmpi(varargin,'full');
    varargin(idx) = [];
    dtmp = thisdesignmethods(this, varargin{:});
  else
    dtmp = d;
  end
  for idx = 1:length(dtmp)
    hf = feval(getdesignobj(getcurrentspecs(this),dtmp{idx}));
    s = getvalidsysobjstructures(hf);
    if isempty(s)
      idxVector = [idxVector idx]; %#ok<AGROW>
    end
  end
  d(idxVector) = [];
end

if nargout,
    varargout = {d};
else
    if ~isempty(type)
        type = sprintf('%s ', type);
    end
    fprintf(1, '\n\n');
    if sysObjFlag
        fprintf('%sDesign Methods that support System objects for class %s (%s):\n', type, ...
          class(this), get(this, 'Specification'));
    else
      fprintf('%sDesign Methods for class %s (%s):\n', type, ...
          class(this), get(this, 'Specification'));
    end
    fprintf(1, '\n\n');
    for indx = 1:length(d),
        disp(d{indx});
    end
    fprintf(1, '\n');
end

% [EOF]
