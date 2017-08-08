classdef (CaseInsensitiveProperties=true) abstractspectrum < matlab.mixin.SetGet & matlab.mixin.Copyable
  %spectrum.abstractspectrum class
  %    spectrum.abstractspectrum properties:
  %       EstimationMethod - Property is of type 'string' (read only)
  %
  %    spectrum.abstractspectrum methods:
  %       checkinputsig -   Return the input vector column'ized & number of channels.
  %       confinterval -  Confidence Interval
  %       copy -   Copy this object.
  %       disp -   Spectrum object display method.
  %       legendstring -   Return the legend string.
  %       loadobj -   Load this object.
  %       psd -   Power Spectral Density (PSD) estimate.
  %       psdopts -   Create an options object for a spectrum object.
  %       reorderprops -   List of properties to reorder.
  %       saopts -   Return options for the spectral analysis commdand-line functions.
  %       saveobj -   Save this object.
  %       thisloadobj -   Load this object.
  %       thispsd -   Power Spectral Density (PSD) estimate.
  %       thissaveobj -   Save this object.
  %       useseglenfornfft -   True for spectral techniques that use the segment
  
  
  properties (SetAccess=protected, Transient, AbortSet, SetObservable, GetObservable)
    %ESTIMATIONMETHOD Property is of type 'string' (read only)
    EstimationMethod = '';
  end
  
  
  methods
    
    function set.EstimationMethod(obj,value)
      % DataType = 'string'
      validateattributes(value,{'char'}, {'row'},'','EstimationMethod')
      obj.EstimationMethod = setestimationmethod(obj,value);
    end
    
    function [x,nchans,msg] = checkinputsig(this,x)
      %CHECKINPUTSIG   Return the input vector column'ized & number of channels.
      %
      % This is a private method.
      
      
      % If its a matrix error out.
      msg = '';
      [lenX,nchans] = size(x);
      xIsMatrix = ~any([lenX,nchans]==1);
      
      if xIsMatrix,
        msg = getString(message('signal:spectrum:MultichannelDatamatricesIsNotSupported'));
        return;
      else
        x = x(:);
        [lenX,nchans] = size(x);
      end
      
    end
    
    function CI = confinterval(this,x,Pxx,W,CL,fs) %#ok<INUSD,STOUT>
      %CONFINTERVAL  Confidence Interval
      
      
      error(message('signal:spectrum:abstractspectrum:confinterval:InternalError'));
      
      
    end
    
    
    function Hcopy = copytheobj(this)
      %COPYTHEOBJ   Copy this object.
      
      % LOADOBJ and COPY perform the same actions.
      Hcopy = loadobj(this);
      
      
    end
    
    
    function disp(this)
      %DISP   Spectrum object display method.
      
      s  = get(this);
      fn = fieldnames(s);
      N = length(fn);
      
      % Reorder so that EstimationMethod is first
      fn = reorderEstimationMethod(this,fn);
      
      for i=1:N,
        snew.(fn{i}) = s.(fn{i});
      end
      
      disp(snew);
      
    end
    
    
    function fn = reorderEstimationMethod(this,fn)
      
      i1 = find(strcmpi(fn,'EstimationMethod'));
      fnE = fn(i1);
      fn(i1) = [];
      fn = [fnE; fn];
      
    end
    
    
    
    function s = legendstring(this)
      %LEGENDSTRING   Return the legend string.
      %
      % This is a private method.
      
      s = [this.EstimationMethod,' Method'];
      
    end
    
    
    
    
    function varargout = psd(this,x,varargin)
      %PSD   Power Spectral Density (PSD) estimate.
      %   Type help spectrum/psd for help.
      
      narginchk(2,12);
      [x,lenX] = checkinputsigdim(x);% Column'izes x if a row vector.
      
      hopts = uddpvparse('dspopts.spectrum',{'psdopts',this,x},varargin{:});
      
      % Call psd of concrete class.
      [hopts,opts] = saopts(this,lenX,hopts); % Opts for spectral analysis fcns.
      [Pxx, W] = thispsd(this,x,opts);  % Produces spectrum with pos freq only!
      
      % %Initialize Confidence Interval if no computation is required
      CI = [];
      %
      % Create a dspdata object and center-DC if necessary.
      %
      % p* fcns return spectra w/ pos freq only, so call centerDC if requested.
      centerDCcache = hopts.CenterDC;
      hopts.CenterDC=false;            % Temporary
      
      % Define valid options for dspdata object.
      propName = getrangepropname(hopts);
      propValue= get(hopts,propName);
      dspdataOpts = {'ConfLevel',hopts.ConfLevel,'ConfInterval',CI,...
        'CenterDC',hopts.CenterDC,'Fs',hopts.Fs,propName,propValue};
      
      hpsd = dspdata.psd(Pxx,W,dspdataOpts{:});
      if centerDCcache,
        centerdc(hpsd);
      end
      hopts.CenterDC = centerDCcache;  % Restore original value.
      
      % Calculation of Confidence Interval
      if(isnumeric(hopts.ConfLevel))
        CL = hopts.ConfLevel;
        Pxx = hpsd.Data;
        W = hpsd.Frequencies;
        
        CI = confinterval(this,x,Pxx,W,CL,hopts.Fs);
        dspdataOpts = {'ConfLevel',hopts.ConfLevel,'ConfInterval',CI,...
          'CenterDC',hopts.CenterDC,'Fs',hopts.Fs,propName,propValue};
        hpsd = dspdata.psd(Pxx,W,dspdataOpts{:});
      end
      
      % Store a spectrum object in the data obj's metadata property.
      hpsd.Metadata.setsourcespectrum(this);
      
      if nargout == 0,
        plot(hpsd);
      else
        varargout{1} = hpsd;
      end
      
    end
    
    function hopts = psdopts(this,x)
      %PSDOPTS   Create an options object for a spectrum object.
      %
      
      % Construct default opts object
      hopts = dspopts.spectrum;
      
      % Defaults.
      isrealX = true;
      
      % Parse input.
      if nargin == 2,
        isrealX = isreal(x);
      end
      
      if ~isrealX,
        hopts.SpectrumType = 'twosided';
      end
      
      
    end
    
    function proplist = reorderprops(this)
      %REORDERPROPS   List of properties to reorder.
      
      proplist = {'EstimationMethod','Nfft'};
      
      
    end
    
    function [hopts,opts] = saopts(this,segLen,hopts)
      %SAOPTS   Return options for the spectral analysis commdand-line functions.
      %
      %   Sets up the options in the correct order to be passed in to the command
      %   line functions pwelch, pmusic, etc.
      
      opts = {hopts.SpectrumType};
      if ~hopts.NormalizedFrequency,
        opts = {hopts.Fs, opts{:}};  % Prepend Fs.
      end
      
      % If Welch use the segment length, instead of the input length.
      if useseglenfornfft(this),
        segLen = this.SegmentLength;
      end
      % Determine numeric value of NFFT if it's set to a string.
      nfft = calcnfft(hopts,segLen);
      
      opts = {nfft,opts{:}};     % Prepend NFFT.
      
      
    end
    
    function s = saveobj(this)
      %SAVEOBJ   Save this object.
      
      s = rmfield(get(this), 'EstimationMethod');
      
      s.class = class(this);
      
      s = setstructfields(s,thissaveobj(this,s));
      
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      
    end
    
    function varargout = thispsd(this,x,varargin)
      %THISPSD   Power Spectral Density (PSD) estimate.
      
      
      error(message('signal:spectrum:abstractspectrum:thispsd:InternalError'));
      
      
    end
    
    function s = thissaveobj(this,s)
      %THISSAVEOBJ   Save this object.
      
      
      % No op.
      
    end
    
    function segLenFlag = useseglenfornfft(this)
      %USESEGLENFORNFFT   True for spectral techniques that use the segment
      %                   length as NFFT.
      
      segLenFlag = false;
      
    end
    
    function varargout = set(obj,varargin)
      
      varargout = signal.internal.signalset(obj,varargin{:});
      varargout = {varargout};
      
    end
    
    function values = getAllowedStringValues(~,prop)
      % This function gives the the valid string values for object properties.
      
      switch prop
        case 'WindowName'
          values = {...
            'Bartlett'
            'Bartlett-Hanning'
            'Blackman'
            'Blackman-Harris'
            'Bohman'
            'Chebyshev'
            'Flat Top'
            'Gaussian'
            'Hamming'
            'Hann'
            'Kaiser'
            'Nuttall'
            'Parzen'
            'Rectangular'
            'Taylor'
            'Triangular'
            'Tukey'
            'User Defined'};
          
        case 'InputType'
          values = {...
            'Vector'
            'DataMatrix'
            'CorrelationMatrix'};
          
        case 'CombineMethod'
          values = {...
            'Adaptive'
            'Eigenvalue'
            'Unity'};
          
        case 'SpecifyDataWindowAs'
          values = {...
            'DPSS'
            'TimeBW'};
          
        case 'SamplingFlag'
          values = {...
            'symmetric'
            'periodic'};
          
        otherwise
          values = {};
      end
      
    end
    
    
  end  %% public methods
  
  
  methods (Static) %% static methods
    
    function this = loadobj(s)
      %LOADOBJ   Load this object.
      
      
      this = feval(s.class);
      
      thisloadobj(this, s);
      
    end
    
    
  end  %% static methods
  
end  % classdef

function checkpercent(val)

if (val < 0) || (val > 100),
  error(message('signal:spectrum:abstractspectrum:schema:InvalidRange'));
end

end  % checkpercent


%--------------------------------------------------------------------------
function str = setestimationmethod(~, str)

if ~license('checkout','signal_toolbox')
  error(message('signal:spectrum:abstractspectrum:schema:LicenseRequired'));
end
end  % setestimationmethod


% [EOF]
