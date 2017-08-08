classdef (CaseInsensitiveProperties=true) music < spectrum.abstractspecwwindow
  %MUSIC   Multiple signal classification (MUSIC) spectral estimator.
  %
  %   SPECTRUM.MUSIC is not recommended.  Use <a href="matlab:help pmusic">pmusic</a> and <a href="matlab:help rootmusic">rootmusic</a> instead.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS) returns a MUSIC pseudospectrum estimator
  %   in H with the number of complex sinusoids set to the numeric value
  %   specified in NSINUSOIDS.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS,SEGMENTLENGTH) returns a MUSIC
  %   pseudospectrum estimator with the number of samples in each segment set
  %   to the numeric value specified by SEGMENTLENGTH.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT) returns a
  %   MUSIC spectral estimator with the numeric value specified by
  %   OVERLAPPERCENT as the percentage of overlap between segments.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,WINNAME)
  %   specifies the window as a string. Use set(H,'WindowName') to get a list
  %   of valid <a href="matlab:set(spectrum.music,'WindowName')">windows</a>.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,...
  %   {WINNAME,WINPARAMETER}) specifies the window in WINNAME and the
  %   parameter value in WINPARAMETER both in a cell array.
  %
  %   NOTE: Depending on the window specified by WINNAME a window parameter
  %   property will be dynamically added to the MUSIC spectral estimator H.
  %   Type "help <WINNAME>" for more details.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,WINNAME,...
  %   THRESHOLD) specifies the cutoff in THRESHOLD for the signal and noise
  %   subspace separation.
  %
  %   H = SPECTRUM.MUSIC(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,WINNAME,...
  %   THRESHOLD,INPUTTYPE) specifies the type of input the MUSIC spectral
  %   estimator accepts. INPUTTYPE can be one of the following strings:
  %       'Vector'  (default)
  %       'DataMatrix'
  %       'CorrelationMatrix'
  %
  %   MUSIC pseudospectrum estimators can be passed to the following
  %   functions along with the data to perform that function:
  %       <a href="matlab:help spectrum/powerest">powerest</a>           - computes the powers and frequencies of sinusoids
  %       <a href="matlab:help spectrum/pseudospectrum">pseudospectrum</a>     - calculates the pseudospectrum
  %       <a href="matlab:help spectrum/pseudospectrumopts">pseudospectrumopts</a> - returns options to calculate the pseudospectrum
  %
  %   EXAMPLE: Spectral analysis of a signal containing complex sinusoids
  %            % and noise.
  %            n = 0:99;
  %            s = exp(1i*pi/2*n)+2*exp(1i*pi/4*n)+exp(1i*pi/3*n)+randn(1,100);
  %            h = spectrum.music(3,20);        % Create a MUSIC spectral estimator.
  %            pseudospectrum(h,s);             % Calculate and plot the pseudospectrum.
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %SEGMENTLENGTH Property is of type 'spt_uint32 user-defined'
    SegmentLength
    %OVERLAPPERCENT Property is of type 'SignalSpectrumPercent user-defined'
    OverlapPercent
    %NSINUSOIDS Property is of type 'udouble user-defined'
    NSinusoids
    %SUBSPACETHRESHOLD Property is of type 'udouble user-defined'
    SubspaceThreshold
    %INPUTTYPE Property is of type 'SignalInputTypeList enumeration: {'Vector','DataMatrix','CorrelationMatrix'}'
    InputType = 'Vector';
  end  
  
  methods  % constructor block
    function this = music(varargin)
      
      narginchk(0,7);
      
      % Set the properties of the object.
      % this = spectrum.music;
      this.EstimationMethod = 'Multiple Signal Classification (MUSIC)';
      initialize(this,varargin{:});
      
    end  % music
    
    
    function set.SegmentLength(obj,value)
      
      validateattributes(value,{'numeric'},{'real','finite','>=',0})
      value = double(int32(value));
      
      obj.SegmentLength = value;
      
    end
    
    function set.OverlapPercent(obj,value)
      
      validateattributes(value,{'numeric'},{'real','finite','>=',0,'<=',100})
      value = double(value);

      obj.OverlapPercent = value;
      
    end
    
    function set.NSinusoids(obj,value)
      % User-defined DataType = 'udouble user-defined'
      
      validateattributes(value,{'numeric'},{'real','finite','integer','>=',0})
      value = double(value);
      
      obj.NSinusoids = setnsinusoids(obj,value);
    end
    
    function set.SubspaceThreshold(obj,value)
      % User-defined DataType = 'udouble user-defined'
      
      validateattributes(value,{'numeric'},{'real','finite','>=',0})
      value = double(value);
      
      obj.SubspaceThreshold = value;
    end
    
    function set.InputType(obj,value)
      % Enumerated DataType = 'SignalInputTypeList enumeration: {'Vector','DataMatrix','CorrelationMatrix'}'
      value = validatestring(value,{'Vector','DataMatrix','CorrelationMatrix'},'','InputType');
      obj.InputType = value;
    end
    
    function [x,nchans,msg] = checkinputsig(this,x)
      %CHECKINPUTSIG   Return the input vector column'ized & number of channels.
      %
      % This is a private method.
      
      msg = '';
      [lenX,nchans] = size(x);
      xIsMatrix = ~any([lenX,nchans]==1);
      
      if strcmpi(this.InputType, 'Vector') & xIsMatrix,
        msg = getString(message('signal:spectrum:MultichannelDatamatricesIsNotSupportedWhenTheInputTy'));
        return;
        
      elseif ~xIsMatrix,
        x = x(:);   % Column'ize it.
        [lenX,nchans] = size(x);
      end
      
      
    end
    
    function varargout = powerest(this,x,Fs)
      %POWEREST   Computes the powers and frequencies of sinusoids.
      %   POW = POWEREST(H,X) returns the vector POW containing the estimates
      %   of the powers of the complex sinusoids contained in the data
      %   represented by X.  H must be a MUSIC object.
      %
      %   X can be a vector or a matrix. If it's a vector it is a signal, if
      %   it's a matrix it may be either a data matrix such that X'*X=R, or a
      %   correlation matrix R.  How X is interpreted depends on the value of the
      %   spectrum object's (H) InputType property.
      %
      %   [POW,W] = POWEREST(...) returns in addition a vector of frequencies W
      %   of the sinusoids contained in X.  W is in units of rad/sample.
      %
      %   [POW,F] = POWEREST(...,Fs) uses the sampling frequency Fs in the
      %   computation and returns the vector of frequencies, F, in Hz.
      %
      %   EXAMPLES:
      %      n = 0:99;
      %      s = exp(1i*pi/2*n)+2*exp(1i*pi/4*n)+exp(1i*pi/3*n)+randn(1,100);
      %      H = spectrum.music(3);
      %      [P,W] = powerest(H,s);
      %
      %   See also ROOTEIG, PMUSIC, PEIG, PMTM, PBURG, PWELCH and CORRMTX.
      
      narginchk(2,3);
      
      if nargin < 3,
        Fs = 1;
      end
      
      P = [this.NSinusoids, this.SubspaceThreshold];
      
      if strcmpi(this.InputType,'CorrelationMatrix'),
        [w,pow] = rootmusic(x,P,'corr',Fs);
      else
        [w,pow] = rootmusic(x,P,Fs);
      end
      
      % If input is single, W will also be single so cast it to double
      w = double(w);
      varargout = {pow,w};
      
    end
    
    function psd(this,varargin)
      %PSD  Overloaded PSD method to produce a meaningful error message.
      
      ClassName = regexp(class(this),'\.','split');
      ClassName = ClassName{end};
      error(message('signal:spectrum:music:psd:NotSupported', ClassName));
      
    end
    
    
    function psdopts(this,varargin)
      %PSDOPTS  Overloaded PSDOPTS method to produce a meaningful error message.
      
      ClassName = regexp(class(this),'\.','split');
      ClassName = ClassName{end};
      error(message('signal:spectrum:music:psdopts:NotSupported', ClassName));
      
    end
    
    function varargout = pseudospectrum(this,x,varargin)
      %PSEUDOSPECTRUM  Pseudospectrum estimate via MUSIC.
      %   Type help spectrum/pseudospectrum for help.
      
      
      error(nargchk(2,10,nargin,'struct'));
      [x,lenX,nchans] = checkinputsigdim(x); % Column'izes x if it's a row vector.
      
      validatesegmentlength(this,nchans);
      
      hopts = uddpvparse('dspopts.pseudospectrum',...
        {'pseudospectrumopts',this,x},varargin{:});
      
      % Call pseudospectrum of concrete class.
      [hopts,opts,P] = setupinputs(this,lenX,hopts);
      [Sxx W] = thispseudospectrum(this,x,opts,P);
      
      % This is necessary because I can't call pmusic without specifying Fs -
      % which would result in the freq with the desired normalized units.
      if hopts.NormalizedFrequency,  % Use normalized frequency.
        W = psdfreqvec('npts',opts{1},'Fs',[],'Range',hopts.SpectrumRange);
      end
      
      %
      % Create a dspdata object and center-DC if necessary.
      %
      % p* fcns return spectrums w/ pos freq only, so call centerDC if requested.
      centerDCcache = hopts.CenterDC;
      hopts.CenterDC=false;
      
      % Define valid options for dspdata object.
      propName = getrangepropname(hopts);
      propValue= get(hopts,propName);
      dspdataOpts = {'CenterDC',hopts.CenterDC,'Fs',hopts.Fs,propName,propValue};
      
      hps = dspdata.pseudospectrum(Sxx,W,dspdataOpts{:});
      if centerDCcache,
        centerdc(hps);
      end
      hopts.CenterDC = centerDCcache;
      
      % Store a spectrum object in the data obj's metadata property.
      hps.Metadata.setsourcespectrum(this);
      
      if nargout == 0,
        plot(hps);
      else
        varargout{1} = hps;
      end
      
    end
    
    function hopts = pseudospectrumopts(this,x)
      %PSEUDOSPECTRUMOPTS   Create an options object for music and eigenvector
      %spectrum objects.
      
      % Construct default opts object
      hopts = dspopts.pseudospectrum;
      
      % Defaults.
      isrealX = true;
      
      % Parse input.
      if nargin == 2,
        isrealX = isreal(x);
      end
      
      if ~isrealX,
        hopts.SpectrumRange = 'whole';
      end
      
      
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      window_thisloadobj(this, s);
      
      this.SegmentLength      = s.SegmentLength;
      this.OverlapPercent     = s.OverlapPercent;
      this.NSinusoids         = s.NSinusoids;
      this.SubspaceThreshold  = s.SubspaceThreshold;
      this.InputType          = s.InputType;
      
      
    end
    
    
    function [Sxx,W] = thispseudospectrum(this,x,opts,P)
      %THISPSEUDOSPECTRUM   Calculate the pseudospectrum via MUSIC.
      %
      % This is a private method.
      
      narginchk(2,4);
      
      if strcmpi(this.InputType,'CorrelationMatrix'),
        [Sxx, W] = pmusic(x,P,'corr',opts{:});
      else
        [Sxx, W] = pmusic(x,P,opts{:});
      end
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
      
    end
    
  end  %% public methods
  
  
  methods (Hidden) %% possibly private or hidden
    function initialize(this,nsinusoids,segmentlength,overlappercent,...
        winName,threshold,varargin)
      %INITIALIZE   Initialize the object with defaults or user input.
      %
      % This is a private method.
      
      % Parse the input to the constructor and initialize the object.
      if nargin < 6,
        threshold = 0;
        if nargin < 5,
          winName = 'Rectangular';
          if nargin < 4,
            overlappercent = 50;
            if nargin < 3,
              if nargin < 2,
                nsinusoids = 2; % Valid for 1 real or 2 complex sinusoids
              end
              segmentlength = 2*nsinusoids;
            end
          end
        end
      end
      
      % Handle the new spec where FFTLength is no longer valid, but we must
      % support it for backwards compatibility.
      validInputTypeStrs = set(this,'InputType');
      
      inputType = 'Vector'; % default
      lenvars = length(varargin);
      if lenvars,
        % Check if new syntax (InputType is the 7th input arg) is being used.
        inputTypeIdxs = regexpi(validInputTypeStrs,varargin{1});
        
        if any([inputTypeIdxs{:}]),  % new syntax
          inputType = varargin{1};
        end
        
        if lenvars ==2,
          inputType = varargin{2};
        end
      end
      
      this.SegmentLength = segmentlength;
      this.OverlapPercent = overlappercent;
      this.NSinusoids = nsinusoids;
      this.SubspaceThreshold = threshold;
      this.InputType = inputType;
      
      setwindownamenparam(this,winName);  % Accepts string or cell array for winName.
      
    end
    
    
  end  
  
end  % classdef

function val = setnsinusoids(this,val)
% Allow integers only but allow 0.

flr_val = floor(val);
if ((val~=0) & flr_val==0) | ( (val~=0) & rem(val,flr_val) ),
  error(message('signal:spectrum:music:schema:MustBeInteger'));
end
end  % setnsinusoids

%--------------------------------------------------------------------------
function checkpercent(val)

if (val < 0) || (val > 100),
  error(message('signal:spectrum:abstractspectrum:schema:InvalidRange'));
end

end



%--------------------------------------------------------------------------
function [hopts,opts,P] = setupinputs(this,lenX,hopts)
% Set up the input arguments to the pmusic and peig functions.

% Generate window; by default it's a rectangular window.
this.Window.Length = this.SegmentLength;  % Window is a private property
win = generate(this.Window);

% Number of overlap samples.
NOverlap = overlapsamples(this);

% Setup options for pmusic and peig.
if hopts.NormalizedFrequency,    Fs = [];
else                             Fs = hopts.Fs;
end

% Determine numeric value of NFFT if it's set to a string.
nfft = calcnfft(hopts,lenX);

% Create cell array of options for the command-line function.
opts = {nfft,Fs,hopts.SpectrumRange,win,NOverlap};

% Complex sinusoids and threshold.
P = [this.NSinusoids, this.SubspaceThreshold];

end

%--------------------------------------------------------------------------
function validatesegmentlength(this,nchans)
% Verify that the segment (window) length equals the number of columns in
% the data matrix.

if strcmpi(this.InputType,'DataMatrix'),
  if this.SegmentLength ~= nchans,
    error(message('signal:spectrum:music:pseudospectrum:invalidSegmentLength', 'SegmentLength'));
  end
end


end


% [EOF]
