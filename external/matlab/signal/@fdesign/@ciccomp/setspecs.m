function setspecs(this, M, N, varargin)
%SETSPECS Set the specs.

%   Copyright 2005-2011 The MathWorks, Inc.

if nargin > 1
  set(this, 'DifferentialDelay', M);
  if nargin > 2
    set(this, 'NumberOfSections', N);
  end
end

if ~isempty(varargin)
  inputCell = varargin;
  % Remove MAGUNITS string from inputCell
  for idx = 1:length(varargin)
    if ischar(inputCell{idx}) && any(strcmpi(inputCell{idx},{'linear','db','squared'}))
      inputCell(idx) = [];
    end
  end
  
  % Parse inputs
  if length(inputCell) >= 2 && isnumeric(inputCell{1}) && ischar(inputCell{2})
    % Input contains (..., R, 'Specstring', ...)
    set(this,'CICRateChangeFactor', inputCell{1});
    varargin(1) = [];
    
  elseif ~ischar(inputCell{1})
    
    switch lower(this.Specification)
      case 'fp,fst,ap,ast'
        % Default constructor with no specstring has been passed as an input
        if length(inputCell) == 6
          % Input contains (..., R,Fp,Fst,Ap,Ast,Fs,...)
          set(this,'CICRateChangeFactor', inputCell{1});
          varargin(1) = [];
        else
          if ~this.NormalizedFrequency
            warning(message('signal:fdesign:ciccomp:setspecs:NormalizedFreqAmbiguity'))
          elseif floor(inputCell{1})==inputCell{1}
            % Third input is an integer so it could be R
            if length(inputCell) == 1 || inputCell{2}<= 1
              % Fourth input is a normalized frequency value so third input is R
              set(this,'CICRateChangeFactor', inputCell{1});
              varargin(1) = [];
            end
          end
        end
      case {'n,fp,ap,ast','n,fc,ap,ast','n,fst,ap,ast','n,fp,fst'}
        lengthValue = 6;
        if strcmpi(this.Specification,'n,fp,fst')
          lengthValue = 5;
        end
        if length(inputCell) == lengthValue
          % Input contains (...,R,N,Fi,Ap,Ast,Fs, ...) or (...,R,N,Fp,Fst,Fs, ...)
          set(this,'CICRateChangeFactor', inputCell{1});
          varargin(1) = [];
        else
          if ~this.normalizedfrequency
            warning(message('signal:fdesign:ciccomp:setspecs:NormalizedFreqAmbiguity'))
          elseif length(inputCell) == 1
            % Input could be R or N, we cannot tell, so need to send a
            % warning. We consider the input to be N for backward
            % compatibility.
            warning(message('signal:fdesign:ciccomp:setspecs:NotFullySpecified'))
          elseif length(inputCell) == 2 && floor(inputCell{1}) == inputCell{1} ...
              && floor(inputCell{2}) == inputCell{2}
            % Inputs 3 and 4 are integers (...,R,N)
            set(this,'CICRateChangeFactor', inputCell{1});
            varargin(1) = [];
          elseif length(inputCell) >= 3 && floor(inputCell{1}) == inputCell{1} ...
              && floor(inputCell{2}) == inputCell{2} && inputCell{3} < 1
            if length(inputCell) == 5
              % Input could be (N,Fp,Ap,Ast,Fs) or (R,N,Fp,Ap,Ast). Need to
              % warn since we have an ambiguity. Assume that input is
              % (N,Fp,Ap,Ast,Fs) for backward compatibility.
              warning(message('signal:fdesign:ciccomp:setspecs:Ambiguous'))
            else
              % Input contains (R,N,Fp,...) no Fs has been specified since
              % length of inputCell is not 5.
              set(this,'CICRateChangeFactor', inputCell{1});
              varargin(1) = [];
            end
          end
        end
    end
  end
end

abstract_setspecs(this, varargin{:});

% [EOF]
