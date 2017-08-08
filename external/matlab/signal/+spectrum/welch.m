classdef (CaseInsensitiveProperties=true) welch < spectrum.periodogram
  %WELCH   Welch spectral estimator.
  %
  %   SPECTRUM.WELCH is not recommended.  Use <a href="matlab:help pwelch">pwelch</a> instead.
  %
  %   H = SPECTRUM.WELCH returns a Welch spectral estimator in H.
  %
  %   H = SPECTRUM.WELCH(WINNAME) returns a Welch spectral estimator in H
  %   with the string specified by WINNAME as the window. Use
  %   set(H,'WindowName') to get a list of valid <a href="matlab:set(spectrum.welch,'WindowName')">windows</a>.
  %
  %   H = SPECTRUM.WELCH({WINNAME,WINPARAMETER}) specifies the window in
  %   WINNAME and the window parameter value in WINPARAMETER in a cell array.
  %
  %   NOTE: Depending on the window specified by WINNAME a window parameter
  %   will be dynamically added to the Welch spectral estimator H. Type "help
  %   <WINNAME>" for more details.
  %
  %   Note also that the default window (Hamming) has a 42.5 dB sidelobe
  %   attenuation. This may mask spectral content below this value (relative
  %   to the peak spectral content). Choosing different windows will enable
  %   you to make tradeoffs between resolution (e.g., using a rectangular
  %   window) and sidelobe attenuation (e.g., using a Hann window). See
  %   WinTool for more details.
  %
  %   H = SPECTRUM.WELCH(WINNAME,SEGMENTLENGTH) specifies the length of each
  %   segment as SEGMENTLENGTH.  The length of the segment allows you to make
  %   tradeoffs between resolution and variance.  A long segment length will
  %   result in better resolution while a short segment length will result in
  %   more averages, and therefore decrease the variance.
  %
  %   H = SPECTRUM.WELCH(WINNAME,SEGMENTLENGTH,OVERLAPPERCENT) specifies the
  %   percentage of overlap between each segment.
  %
  %   Welch estimators can be passed to the following functions along with
  %   the data to perform that function:
  %       <a href="matlab:help spectrum/msspectrum">msspectrum</a>     - calculates the Mean-squared Spectrum (MSS)
  %       <a href="matlab:help spectrum/msspectrumopts">msspectrumopts</a> - returns options to calculate the MSS
  %       <a href="matlab:help spectrum/psd">psd</a>            - calculates the PSD
  %       <a href="matlab:help spectrum/psdopts">psdopts</a>        - returns options to calculate the PSD
  %
  %   EXAMPLE: Spectral analysis of a signal that contains a 200Hz cosine
  %            % plus noise.
  %            Fs = 1000;   t = 0:1/Fs:.296;
  %            x = cos(2*pi*t*200)+randn(size(t));
  %            h = spectrum.welch;                  % Create a Welch spectral estimator.
  %            psd(h,x,'Fs',Fs);                    % Calculate and plot the PSD.
  
  properties (AbortSet, SetObservable, GetObservable)
    %SEGMENTLENGTH Property is of type 'spt_uint32 user-defined'
    SegmentLength
    %OVERLAPPERCENT Property is of type 'SignalSpectrumPercent user-defined'
    OverlapPercent
  end
  
  
  methods  % constructor block
    function this = welch(varargin)
      
      narginchk(0,3);
      
      % Create default welch object.
      % this = spectrum.welch;
      
      winName = 'Hamming';
      if nargin >= 1,
        winName = varargin{1};
      end
      setwindownamenparam(this,winName);  % Accepts string or cell array for winName.
      
      % Parse the rest of the inputs.
      paramCell = {'SegmentLength','OverlapPercent'};
      valCell = {64,50};  % Default values for corresponding properties above.
      
      % Override default values with user input.  Exclude window and fftlength.
      if nargin>=2,
        valCell{1}=varargin{2};
        if nargin>=3,
          valCell{2}=varargin{3};
        end
      end
      
      % Set the properties of the object.
      this.EstimationMethod = 'Welch';
      
      for i = 1:length(paramCell)
        this.(paramCell{i}) = valCell{i};
      end
      
    end  % welch
    
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
    
  end   % set and get functions
  
  methods  %% public methods
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      window_thisloadobj(this, s);
      
      this.SegmentLength  = s.SegmentLength;
      this.OverlapPercent = s.OverlapPercent;
      
      
    end
    
    function [Sxx,W] = thismsspectrum(this,x,opts)
      %THISMSSPECTRUM   Mean-square Spectrum via Welch's method.
      %
      % This is a private method.
      
      narginchk(2,3);
      
      % Generate Window vector.
      this.Window.Length = this.SegmentLength;  % Window is a private property
      win = generate(this.Window);
      
      NOverlap = overlapsamples(this);
      
      [Sxx, W] = pwelch(x,...
        win,...
        NOverlap,...
        opts{:},...     % NFFT, Fs(?), and SpectrumType
        'ms');          % Window compensation produces correct peak heights.
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
    end
    
    function [Pxx,W] = thispsd(this,x,opts)
      %THISPSD Calculate the power spectral density via Welch's method.
      %
      % This is a private method.
      
      narginchk(2,3);
      
      % Generate window vector.
      this.Window.Length = this.SegmentLength;  % Window is a private property
      win = generate(this.Window);
      
      NOverlap = overlapsamples(this);
      
      % Calculate PSD.
      [Pxx, W] = pwelch(x,...
        win,...
        NOverlap,...
        opts{:}); % NFFT, Fs(?), and SpectrumType
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
      
    end
    
    function segLenFlag = useseglenfornfft(this)
      %USESEGLENFORNFFT   True for spectral techniques that use the segment
      %                   length as NFFT.
      
      segLenFlag = true;
      
    end
        
    function  validatenoverlap(this, NOverlap,N)
      %VALIDATENOVERLAP Validate the noverlap
      
      if NOverlap>N-1,
        error(message('signal:spectrum:welch:validatenoverlap:InvalidOverlapPercent', num2str( 100*(N - 1)/N ), N));
      end      
      
    end
    
  end  %% public methods
  
end  % classdef

