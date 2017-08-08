classdef periodogram < spectrum.abstractspecwwindow
  %PERIODOGRAM   Periodogram spectral estimator.
  %
  %   SPECTRUM.PERIODOGRAM is not recommended.  Use <a href="matlab:help periodogram">periodogram</a> instead.
  %
  %   H = SPECTRUM.PERIODOGRAM returns a periodogram spectral estimator in H.
  %
  %   H = SPECTRUM.PERIODOGRAM(WINNAME) returns a periodogram spectral
  %   estimator in H with the string specified by WINNAME as the window. Use
  %   set(H,'WindowName') to get a list of valid <a href="matlab:set(spectrum.periodogram,'WindowName')">windows</a>.
  %
  %   H = SPECTRUM.PERIODOGRAM({WINNAME,WINPARAMETER}) specifies the window
  %   in WINNAME and the window parameter value in WINPARAMETER both in a
  %   cell array.
  %
  %   NOTE: Depending on the window specified by WINNAME a window parameter
  %   will be dynamically added to the periodogram spectral estimator H. Type
  %   "help <WINNAME>" for more details.
  %
  %   Note that the default window (rectangular) has a 13.3 dB sidelobe
  %   attenuation. This may mask spectral content below this value (relative
  %   to the peak spectral content). Choosing different windows will enable
  %   you to make tradeoffs between resolution (e.g., using a rectangular
  %   window) and sidelobe attenuation (e.g., using a Hann window). See
  %   WinTool for more details.
  %
  %   Periodogram estimators can be passed to the following functions along
  %   with the data to perform that function:
  %       <a href="matlab:help spectrum/msspectrum">msspectrum</a>     - calculates the Mean-squared Spectrum (MSS)
  %       <a href="matlab:help spectrum/msspectrumopts">msspectrumopts</a> - returns options to calculate the MSS
  %       <a href="matlab:help spectrum/psd">psd</a>            - calculates the PSD
  %       <a href="matlab:help spectrum/psdopts">psdopts</a>        - returns options to calculate the PSD
  %
  %   EXAMPLE: Spectral analysis of a complex signal plus noise.
  %      Fs = 1000;   t = 0:1/Fs:.296;
  %      x = exp(1i*2*pi*200*t)+randn(size(t));
  %      h = spectrum.periodogram;      % Create a periodogram spectral estimator.
  %      psd(h,x,'Fs',Fs);              % Calculates and plots the two-sided PSD.
  
  %   Author(s): P. Pacheco
  
  
  
  methods  % constructor block
    function this = periodogram(winName)
      
      
      narginchk(0,1);
      
      % Create default periodogram object.
      % this = spectrum.periodogram;
      this  = this@spectrum.abstractspecwwindow;
      
      
      if nargin < 1,
        winName = 'rectangular';
      end
      
      % Set the properties of the object.
      this.EstimationMethod = 'Periodogram';
      
      setwindownamenparam(this,winName);  % Accepts string or cell array for winName.
      
      
    end  % periodogram
    
    
    function disp(this)
      %DISP Spectrum object display method.
      
      s  = get(this);
      fn = fieldnames(s);
      N  = length(fn);
      
      props = propstoaddtospectrum(this.Window);
      
      idx1 = find(strcmpi(fn,'WindowName'));
      idx2 = [];
      if ~isempty(props),
        for k=1:length(props),
          idx2(k) = find(strcmpi(fn,props{k}));
        end
      end
      % Reorder the fields so that windowname and windowparam are last.
      fn1 = fn([idx1 idx2]);
      fn([idx1 idx2]) = [];      
      fn = [fn; fn1];
      
      % Reorder so that EstimationMethod is first
      fn = reorderEstimationMethod(this,fn);

      for i=1:N,
        snew.(fn{i}) = getfield(s, fn{i});
      end
      
      % Make sure the two NFFT properties are together.
      props = reorderprops(this);
      if any(strcmpi(fn,'Nfft')),
        snew = reorderstructure(snew,props{:});
      end
      disp(snew)
      
    end
    
    function varargout = msspectrum(this,x,varargin)
      %MSSPECTRUM   Mean-square Spectrum estimate.
      %   Type help spectrum/msspectrum for help.
      
      narginchk(2,12);
      [x,lenX] = checkinputsigdim(x);% Column'izes x if a row vector.
      
      hopts = uddpvparse('dspopts.spectrum',{'psdopts',this,x},varargin{:});
      
      % Call msspectrum of concrete class.
      [hopts,opts] = saopts(this,lenX,hopts); % Opts for spectral analysis fcns.
      [Sxx, W] = thismsspectrum(this,x,opts);
      
      % %Initialize Confidence Interval if no computation is required
      CI = [];%
      % Create a dspdata object and center-DC if necessary.
      %
      % p* fcns return spectra w/ pos freq only, so centerDC if requested.
      centerDCcache = hopts.CenterDC;
      hopts.CenterDC=false;
      
      % Define valid options for dspdata object.
      propName = getrangepropname(hopts);
      propValue= get(hopts,propName);
      dspdataOpts = {'ConfLevel',hopts.ConfLevel,'ConfInterval',CI,...
        'CenterDC',hopts.CenterDC,'Fs',hopts.Fs,propName,propValue};
      
      hmss = dspdata.msspectrum(Sxx,W,dspdataOpts{:});
      if centerDCcache,
        centerdc(hmss);
      end
      hopts.centerDC = centerDCcache;
      
      % Calculation of Confidence Interval
      if(isnumeric(hopts.ConfLevel))
        CL = hopts.ConfLevel;
        Sxx = hmss.Data;
        W = hmss.Frequencies;
        
        CI = confinterval(this,x,Sxx,W,CL,hopts.Fs);
        dspdataOpts = {'ConfLevel',hopts.ConfLevel,'ConfInterval',CI,...
          'CenterDC',hopts.CenterDC,'Fs',hopts.Fs,propName,propValue};
        hmss = dspdata.msspectrum(Sxx,W,dspdataOpts{:});
      end
      
      % Store a spectrum object in the data obj's metadata property.
      hmss.Metadata.setsourcespectrum(this);
      
      if nargout == 0,
        plot(hmss);
      else
        varargout{1} = hmss;
      end
      
    end
    
    function hopts = msspectrumopts(this,varargin)
      %MSSPECTRUMOPTS   Create an options object for a spectrum object.
      %
      
      hopts = psdopts(this,varargin{:});
      
    end
    
    function [Sxx,W] = thismsspectrum(this,x,opts)
      %THISMSSPECTRUM   Mean-square Spectrum via periodogram.
      %
      % This is a private method.

      narginchk(2,3);
      
      % Generate window vector.
      this.Window.Length = length(x);
      win = generate(this.Window);
      
      [Sxx, W] = periodogram(x,...
        win,...
        opts{:},...     % NFFT, Fs(?), and SpectrumType
        'ms');          % Compensate for window so that it produces correct peak heights.
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
      
    end
    
    function [Pxx,W] = thispsd(this,x,opts)
      %THISPSD Calculate the power spectral density (PSD) via periodogram.
      %
      % This is a private method.

      narginchk(2,3);
      
      % Generate window.
      this.Window.Length = length(x);
      win = generate(this.Window);
      
      [Pxx, W] = periodogram(x,...
        win,...
        opts{:});  % NFFT, Fs(?), and SpectrumType
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
    end
    
  end  %% public methods
  
end  % classdef

