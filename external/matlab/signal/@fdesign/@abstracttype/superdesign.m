function varargout = superdesign(this, method, varargin)
%SUPERDESIGN   Design the filter.
%   This method is used to enable subclasses that overwrite design method to 
%   call design method of the super class.

%   Copyright 2009-2011 The MathWorks, Inc.

% This should be a private method.

if nargin > 1 && strcmpi(method,'systemobject')
  method = 'default';
  varargin = [{'SystemObject'} varargin];
end

sysObjArgs = {};
sysObjArgsDesign = {};

if nargin < 2
    method = 'default';
else
    [varargin sysObjFlag] = parsesysobj(this,'design',varargin{:});
    
    if sysObjFlag
      % Need sysObjArgs since methods might not take SystemObject as an
      % argument depending on the license.
      sysObjArgs = {'SystemObject',sysObjFlag}; 
      if ~isempty(varargin) && isstruct(varargin{1})
        % Cannot pass a structure and pv-pairs to the design method. So
        % just add the SystemObject field if options are in structure form
        % and leave sysObjArgsDesign empty.
        varargin{1}.SystemObject = true;
      else
         sysObjArgsDesign = {'SystemObject',sysObjFlag}; 
      end
    end
    validflags = {'all', 'allfir', 'alliir', 'default', 'fir', 'iir'};
    if ~isdesignmethod(this, method) && ~any(strcmpi(method, validflags))
        varargin = [{method}, varargin];
        method = 'default';
    end
end

method = lower(method);

switch method
    case 'all'
        d = designmethods(this,sysObjArgs{:});
        Hd = {};
        for indx = 1:length(d)
            try
                Hd{end+1} = feval(d{indx}, this, sysObjArgs{:}); %#ok<AGROW>
            catch me %#ok<NASGU>
                warning(message('signal:fdesign:abstracttype:superdesign:FailedDesign', d{ indx }));
            end
        end
        
        if ~sysObjFlag
          Hd = [Hd{:}];
        end

        varargout = {Hd};
    case {'alliir', 'allfir'}
        d = designmethods(this, method(4:end),sysObjArgs{:});
        if isempty(d)
            error(message('signal:fdesign:abstracttype:superdesign:invalidMethod', upper( method( 4:end ) ), this.SpecificationType));
        end
        Hd = {};
        for indx = 1:length(d)
            try
                Hd{end+1} = feval(d{indx}, this,sysObjArgs{:}); %#ok<AGROW>
            catch me %#ok<NASGU>
                warning(message('signal:fdesign:abstracttype:superdesign:FailedDesign', d{ indx }));
            end
        end
        if ~sysObjFlag
          Hd = [Hd{:}];
        end
        varargout = {Hd};
    case 'default'
        d = defaultmethod(this);
        if nargin > 1 && ~isempty(varargin) && ~isstruct(varargin{1}) && ~isnumeric(varargin{1})  && ...
                ~any(strcmpi(fieldnames(designopts(this,d)),varargin{1})),
            % An invalid method was specified
            error(message('signal:fdesign:abstracttype:superdesign:invalidDesignMethod', varargin{ 1 }));
        end
        if nargin > 1 && ~isempty(varargin) && rem(length(varargin),2) && ~isstruct(varargin{1}) ...
                && ~isnumeric(varargin{1}),
            error(message('signal:fdesign:abstracttype:superdesign:invalidPVpairs'));
        end
        Hd = feval(d, this, varargin{:},sysObjArgsDesign{:});
        varargout = {Hd};
    case {'fir', 'iir'}
        d = designmethods(this, method,sysObjArgs{:});
        if isempty(d)
            error(message('signal:fdesign:abstracttype:superdesign:invalidMethod', upper( method ), this.SpecificationType));
        end
        if strcmpi(method, 'fir') && any(strcmpi(d, 'equiripple'))
            d = 'equiripple';
        elseif strcmpi(method, 'iir') && any(strcmpi(d, 'ellip'))
            d = 'ellip';
        else
            d = d{1};
        end
        Hd = feval(d, this,sysObjArgs{:});
        varargout = {Hd};
    otherwise
        % hiddenmethods is there for backwards compatibility (for example
        % butter needs to work for 'N,Fc' (lowpass) although it is no
        % longer a valid designmethod (we changed it to 'N,F3dB'))
        if any(strcmpi(method, designmethods(this))) || ...
                any(strcmpi(method, hiddenmethods(this)))
            % Check for valid p-v pairs, but first remove any possible
            % options structure
            if nargin > 2,
                args = varargin;
                for k = 1:length(args),
                    if isstruct(args{k}),
                        args(k) = [];
                    end
                end
                if rem(length(args),2) && ~isnumeric(args{1}),
                    error(message('signal:fdesign:abstracttype:superdesign:invalidPVpairs'));
                end
            end
            Hd = feval(method, this, varargin{:},sysObjArgsDesign{:});
            varargout = {Hd};
            if any(strcmpi(method, hiddenmethods(this))),
                warning(message('signal:fdesign:abstracttype:superdesign:Obsolete', 'butter', 'butter', 'Fc', 'F3dB'));
            end
        else
            error(message('signal:fdesign:abstracttype:superdesign:invalidDesign', upper( method ), this.Specification));
        end
end

if ~nargout
    Hd = varargout{1};
    varargout = {};
    if this.NormalizedFrequency,
        inputs = {'NormalizedFrequency', 'On'};
    else
        inputs = {'Fs', this.Fs};
    end

    inputs = [inputs, {'DesignMask', 'on'}];

    h = fvtool(Hd, inputs{:});
    switch method
        case 'all'
            strs = designmethods(this, 'all', 'full');
        case {'allfir', 'alliir'}
            strs = designmethods(this, method(4:end), 'full');
        otherwise
            strs = {};
    end
    if ~isempty(strs), legend(h, strs{:}); end
end

% [EOF]
