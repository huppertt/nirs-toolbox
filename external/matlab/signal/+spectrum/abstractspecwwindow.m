classdef (CaseInsensitiveProperties=true) abstractspecwwindow < sigio.dyproputil & spectrum.abstractspectrum
  %spectrum.abstractspecwwindow class
  %   spectrum.abstractspecwwindow extends spectrum.abstractspectrum.
  %
  %    spectrum.abstractspecwwindow properties:
  %       EstimationMethod - Property is of type 'string' (read only)
  %       WindowName - Property is of type 'SignalSpectrumWindowList enumeration: {'Bartlett','Bartlett-Hanning','Blackman','Blackman-Harris','Bohman','Chebyshev','Flat Top','Gaussian','Hamming','Hann','Kaiser','Nuttall','Parzen','Rectangular','Taylor','Triangular','Tukey','User Defined'}'
  %
  %    spectrum.abstractspecwwindow methods:
  %       confinterval -  Confidence Interval for Periodogram and Welch methods.
  %       disp - Spectrum object display method.
  %       overlapsamples - Return the number of overlap samples.
  %       setwindowname -   Set function for the WindowName property.
  %       setwindownamenparam -   Sets the WindowName and parameter if specified.
  %       setwinobj - Sets the window property (object) of the response object.
  %       thisloadobj -   Load this object.
  %       validatenoverlap - Validate the noverlap
  %       window_thisloadobj -   Load this object.
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %WINDOWNAME Property is of type 'SignalSpectrumWindowList enumeration: {'Bartlett','Bartlett-Hanning','Blackman','Blackman-Harris','Bohman','Chebyshev','Flat Top','Gaussian','Hamming','Hann','Kaiser','Nuttall','Parzen','Rectangular','Taylor','Triangular','Tukey','User Defined'}'
    WindowName = 'rectangular';
  end
  
  properties (Access=protected, AbortSet, SetObservable, GetObservable, Hidden)
    %WINDOW Property is of type 'MATLAB array'
    Window = [];
    %WINDOWPARAMETERS Property is of type 'mxArray'
    WindowParameters = [];
  end
  
  
  methods
    function set.Window(obj,value)
      obj.Window = setwinobj(obj,value);
    end
    
    function set.WindowName(obj,value)
      % Enumerated DataType = 'SignalSpectrumWindowList enumeration: {'Bartlett','Bartlett-Hanning','Blackman','Blackman-Harris','Bohman','Chebyshev','Flat Top','Gaussian','Hamming','Hann','Kaiser','Nuttall','Parzen','Rectangular','Taylor','Triangular','Tukey','User Defined'}'
      [~,winNames] = findallwinclasses;
      value = validatestring(value,winNames,'','WindowName');
      obj.WindowName = setwindowname(obj,value);
    end
    
    function set.WindowParameters(obj,value)
      obj.WindowParameters = value;
    end
    
    function this = abstractspecwwindow
      this.setwindowname('Rectangular');
    end
    
    function CI = confinterval(this,x,Pxx,W,CL,fs)
      %CONFINTERVAL  Confidence Interval for Periodogram and Welch methods.
      %   CI = CONFINTERVAL(THIS,X,PXX,W,CL) calculates the confidence
      %   interval CI for spectrum estimate PXX based on confidence level CL. THIS is a
      %   spectrum object and W is the frequency vector. X is the data used for
      %   computing the spectrum estimate PXX.
      %
      %   Reference: D.G. Manolakis, V.K. Ingle and S.M. Kagon,
      %   Statistical and Adaptive Signal Processing,
      %   McGraw-Hill, 2000, Chapter 5
      
      name = this.EstimationMethod;
      L = length(x);
      
      switch lower(name)
        case{'periodogram'}
          k = 1;
        case{'welch'}
          SegLen = this.SegmentLength;
          Per = this.OverlapPercent;
          Noverlap = Per*SegLen/100;
          k = (L-Noverlap)/(SegLen-Noverlap);
          k = fix(k);
      end
      
      c = privatechi2conf(CL,k);
      CI = Pxx*c;
      
      % DC and Nyquist bins have only one degree of freedom for real signals
      if isreal(x)
        realConf = privatechi2conf(CL,k/2);
        CI(W == 0,:) = Pxx(W == 0) * realConf;
        if isempty(fs) || ~isnumeric(fs)
          CI(W==pi,:) = Pxx(W==pi) * realConf;
        else
          CI(W==fs/2,:) = Pxx(W==fs/2) * realConf;
        end
      end
      
    end
    
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
    
    function NOverlap = overlapsamples(this)
      %OVERLAPSAMPLES Return the number of overlap samples.
      
      
      N = this.SegmentLength;
      NOverlap = ceil((this.OverlapPercent/100) * N);
      
      validatenoverlap(this,NOverlap,N);
      
    end
    
    function winName = setwindowname(this,winName)
      %SETWINDOWNAME   Set function for the WindowName property.
      %
      % This function updates the Window private property that stores a sigwin
      % object whenever the property WindowName is set.
      %
      % This is a private method.
      
      
      this.Window = getwinobject(winName);
      
    end
    
    
    function setwindownamenparam(this,varargin)
      %SETWINDOWNAMENPARAM   Sets the WindowName and parameter if specified.
      %   This functions allows us to use either a string or cell array to
      %   set the WindowName property.
      
      winName = varargin{1};
      if iscell(winName), % Window parameter specified.
        this.WindowName = winName{1};
        paramName = propstoaddtospectrum(this.Window);
        for k=1:min(length(paramName),length(winName)-1),
          set(this,paramName{k},winName{k+1});
        end
      else
        this.WindowName = winName;
      end
      
    end
    
    function hwinObj = setwinobj(this,hwinObj)
      % SETWINOBJ Sets the window property (object) of the response object.
      
      
      if isempty(hwinObj),
        return;
      end
      
      p = this.WindowParameters;
      if ~isempty(this.Window)
        
        % Get a list of the properties added by the old window.
        props2remove = propstoaddtospectrum(this.Window);
        
        % Cache all of the old values in the WindowParameters structure.
        for indx = 1:length(props2remove)
          p.(props2remove{indx}) = this.(props2remove{indx});
        end
        this.WindowParameters = p;
      end
      
      % Remove the properties for the old window.
      rmprops(this, this.Window);
      
      % If there are no propstoadd, do nothing
      props2add = propstoaddtospectrum(hwinObj);
      if ~isempty(props2add),
        
        % Add the properties from the window object to the spectrum.
        hp = addprops(this,hwinObj,props2add{:});
        
        % Check if any of the properties that we have added for this window
        % have previously been added to the spectrum object.
        for indx = 1:length(props2add)
          if isfield(p, props2add{indx})
            set(this, props2add{indx}, p.(props2add{indx}));
          end
        end
      end
      
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      
      window_thisloadobj(this, s);
      
    end
    
    function  validatenoverlap(this, NOverlap,N)
      %VALIDATENOVERLAP Validate the noverlap
      
      % No op.
      
      
    end
    
    function window_thisloadobj(this, s)
      %WINDOW_THISLOADOBJ   Load this object.
      
      
      this.WindowName = s.WindowName;
      
      p = propstoaddtospectrum(this.Window);
      
      for indx = 1:length(p)
        set(this, p{indx}, s.(p{indx}));
      end
      
    end
    

    
  end  %% public methods
  
end  % classdef

