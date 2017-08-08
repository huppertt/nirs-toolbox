function dopts = designoptions(this, dmethod, varargin)
%DESIGNOPTIONS

%   Copyright 2005-2014 The MathWorks, Inc.

if nargin < 2,
    error(message('signal:fdesign:abstractmultirate2:designoptions:notEnoughInputs', inputname( 1 )));
end

% Parse the SystemObject input
sysObjFlag = validatedesignoptionssysobjinput(this,varargin{:});

dopts = designoptions(this.CurrentFDesign, dmethod);

if isempty(fieldnames(dopts))
    return;
end

rcf = getratechangefactors(this);

if ~any(rcf == 1)
    % RSRC
    if ~any(strcmpi('firsrc', dopts.FilterStructure))
        dopts.FilterStructure = {'firsrc'};
    end
elseif rcf(1) < rcf(2)
    % Decimator
    if ~strcmpi(this.Response,'cic')
      if isempty(intersect({'firdecim', 'firtdecim'}, ...
          dopts.FilterStructure)) || ...
          isequal({'firdecim','firinterp','firtdecim'},...
          intersect(dopts.FilterStructure,{'firinterp','firdecim','firtdecim'}))
        dopts.FilterStructure = {'firdecim', 'firtdecim'};
      end
      if any(strcmpi(designmethods(this, 'iir'), dmethod))
        dopts.FilterStructure = {'iirdecim','iirwdfdecim'};
        % Remove SystemObject options - no multirate support of
        % allpass-based IIR structures
        if isfield(dopts, 'SystemObject')
            dopts = rmfield(dopts, {'SystemObject', 'DefaultSystemObject'});
        end
      end
    end
else
    % Interpolator
    if ~strcmpi(this.Response,'cic')      
      if isempty(intersect({'firinterp', 'fftfirinterp'}, ...
          dopts.FilterStructure))
        dopts.FilterStructure = {'firinterp', 'fftfirinterp'};
      elseif isequal({'firdecim','firinterp','firtdecim'},...
          intersect(dopts.FilterStructure,{'firinterp','firdecim','firtdecim'}))
        dopts.FilterStructure = {'firinterp'};
      end
      if any(strcmpi(designmethods(this, 'iir'), dmethod))
        dopts.FilterStructure = {'iirinterp','iirwdfinterp'};
        % Remove SystemObject options - no multirate support of
        % allpass-based IIR structures
        if isfield(dopts, 'SystemObject')
            dopts = rmfield(dopts, {'SystemObject', 'DefaultSystemObject'});
        end
      end
    end
end

if isfield(dopts,'FilterStructure')
  dopts.DefaultFilterStructure = dopts.FilterStructure{1};
  
  % If SystemObject parameter was passed as an input with a value of true,
  % then remove structures that are not supported by System objects.
  if sysObjFlag
    hf = feval(getdesignobj(getcurrentspecs(this),dmethod));
    supportedStructs = getsysobjsupportedstructs(hf);
    [s k ~] = intersect(dopts.FilterStructure,supportedStructs);
    if isempty(s)
      % Error out if none of the structures is supported by System objects
      error(message('signal:fdesign:basecatalog:NoSupportedSysObj',dmethod,'''SystemObject'''))
    end
    dopts.FilterStructure = dopts.FilterStructure(sort(k));
  end
end

% [EOF]
