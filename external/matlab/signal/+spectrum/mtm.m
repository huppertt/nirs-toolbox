classdef (CaseInsensitiveProperties=true) mtm < spectrum.abstractspectrum & sigio.dyproputil
  %MTM   Thomson multitaper method (MTM) power spectral density (PSD) estimator.
  %
  %   SPECTRUM.MTM is not recommended.  Use <a href="matlab:help pmtm">pmtm</a> instead.
  %
  %   H = SPECTRUM.MTM(TIMEBW) returns an MTM PSD estimator with the
  %   time-bandwidth product set to TIMEBW. The time-bandwidth product is of
  %   the discrete prolate spheroidal sequences (or Slepian sequences) used
  %   as data windows.
  %
  %   H = SPECTRUM.MTM(DPSS,CONCENTRATIONS) returns an mtm spectral estimator
  %   with the discrete prolate spheroidal sequences and their concentrations
  %   set to DPSS and CONCENTRATIONS respectively.  Type "help dpss" for more
  %   information on these two input arguments.
  %
  %   NOTE: Specifying DPSS and CONCENTRATIONS when constructing the MTM
  %   estimator automatically changes the value of the SpecifyDataWindowAs
  %   property to 'DPSS' from its default value 'TimeBW'.
  %
  %   H = SPECTRUM.MTM(...,COMBINEMETHOD) specifies the algorithm for
  %   combining the individual spectral estimates. COMBINEMETHOD can be one
  %   of the following strings:
  %      'Adaptive'   - Thomson's adaptive non-linear combination
  %      'Eigenvalue' - linear combination with eigenvalue weights.
  %      'Unity'      - linear combination with unity weights.
  %
  %   MTM PSD estimators can be passed to the following functions along with
  %   the data to perform that function:
  %       <a href="matlab:help spectrum/psd">psd</a>     - calculates the PSD
  %       <a href="matlab:help spectrum/psdopts">psdopts</a> - returns options to calculate the PSD
  %
  %   EXAMPLES:
  %
  %   % Example 1: A cosine of 200Hz plus noise.
  %                Fs = 1000;   t = 0:1/Fs:.3;
  %                x = cos(2*pi*t*200)+randn(size(t));
  %                h = spectrum.mtm(3.5); % Specify the time-bandwidth product
  %                                       % when creating an MTM spectral estimator.
  %                psd(h,x,'Fs',Fs);      % Calculate and plot the PSD.
  %
  %   % Example 2: This is the same example as above, but we'll specify the
  %                % data tapers and their concentrations instead of the time BW product.
  %                Fs = 1000;   t = 0:1/Fs:.3;
  %                x = cos(2*pi*t*200)+randn(size(t));
  %                [E,V] = dpss(length(x),3.5);
  %                h = spectrum.mtm(E,V);    % Specify DPSS and concentrations
  %                                          % when creating the MTM spectral estimator.
  %                psd(h,x,'Fs',Fs);         % Calculate and plot the PSD.
  
  properties (SetObservable, GetObservable)
    %SPECIFYDATAWINDOWAS Property is of type 'SignalSpectrumDataWinSpecMode enumeration: {'DPSS','TimeBW'}'
    SpecifyDataWindowAs = 'TimeBW';
  end
  
  properties (AbortSet, SetObservable, GetObservable)
    %COMBINEMETHOD Property is of type 'SignalSpectrumCombineMethodList enumeration: {'Adaptive','Eigenvalue','Unity'}'
    CombineMethod = 'Adaptive';
  end
  
  methods
    function this = mtm(varargin)
      
      narginchk(0,4);
      
      % Create default MTM object.
      % this = spectrum.mtm;
      
      % Defaults.
      NW = 4;
      combinemethod = 'Adaptive';
      SpecifyDataWindowAs = this.SpecifyDataWindowAs;  % User default set in schema.
      E = [];
      V = [];
      
      if nargin >= 1,
        SpecifyDataWindowAs ='TimeBW';
        NW = varargin{1};
      end
      
      if ~any(size(NW) == 1);
        % DPSS and Concentrations were specified.
        SpecifyDataWindowAs = 'DPSS';
        E = NW;
        if nargin >= 2,
          V = varargin{2};
        else
          error(message('signal:spectrum:mtm:mtm:SignalErr'));
        end
        % Enable the code below to depend on the same number of inputs.
        varargin(2) = [];
      end
      
      % Override default values if user specified values.
      nargs = length(varargin);
      if nargs >= 2,
        combinemethod = varargin{2};
      end
      
      % Set the properties of the object.
      this.EstimationMethod = 'Thompson Multitaper';
      this.SpecifyDataWindowAs = SpecifyDataWindowAs;
      this.CombineMethod = combinemethod;
      
      if strcmpi(SpecifyDataWindowAs,'TimeBW'),
        this.TimeBW = NW;
      else
        this.DPSS = E;
        this.Concentrations = V;
      end
      
      
    end  % mtm
    
    function set.SpecifyDataWindowAs(obj,value)
      % Enumerated DataType = 'SignalSpectrumDataWinSpecMode enumeration: {'DPSS','TimeBW'}'
      value = validatestring(value,{'DPSS','TimeBW'},'','SpecifyDataWindowAs');
      obj.SpecifyDataWindowAs = setdatawindowspecmode(obj,value);
    end
    
    function set.CombineMethod(obj,value)
      % Enumerated DataType = 'SignalSpectrumCombineMethodList enumeration: {'Adaptive','Eigenvalue','Unity'}'
      value = validatestring(value,{'Adaptive','Eigenvalue','Unity'},'','CombineMethod');
      obj.CombineMethod = value;
    end
    
    
    function CI = confinterval(this,x,Pxx,~,CL,~)
      %CONFINTERVAL  Confidence Interval for MTM method.
      %   CI = CONFINTERVAL(THIS,X,PXX,W,CL,FS) calculates the confidence
      %   interval CI for spectrum estimate PXX based on confidence level CL. THIS is a
      %   spectrum object and W is the frequency vector. X is the data used for
      %   computing the spectrum estimate PXX.
      %
      %
      %   References:
      %     [1] Thomson, D.J."Spectrum estimation and harmonic analysis."
      %         In Proceedings of the IEEE. Vol. 10 (1982). Pgs 1055-1096.
      %     [2] Percival, D.B. and Walden, A.T., "Spectral Analysis For Physical
      %         Applications", Cambridge University Press, 1993, pp. 368-370.
      
      name = this.SpecifyDataWindowAs;
      N = length(x);
      switch lower(name)
        case{'timebw'}
          NW = this.TimeBW;
          k = min(round(2*NW),N);
          k = max(k-1,1);
        case{'dpss'}
          k = length(this.Concentrations);
      end
      
      c = privatechi2conf(CL,k);
      CI = Pxx*c;
      
      
    end
    
    function combineMethodStr = getcombinemethodstr(this)
      %GETCOMBINEMETHODSTR   Get string accepted by the pmtm function.
      %
      % This is a private method.
      
      % Convert CombineMethod enum type to strings accepted by the function.
      combinemethod = lower(this.CombineMethod);
      if strcmpi(combinemethod,'adaptive'),
        combineMethodStr = 'adapt';
        
      elseif strcmpi(combinemethod,'unity'),
        combineMethodStr = 'unity';
        
      elseif strcmpi(combinemethod,'eigenvector'),
        combineMethodStr = 'eigen';
      end
      
      
    end
    
    function s = legendstring(this)
      %LEGENDSTRING Return the legend string.
      %
      % This is a private method.
      
      s = 'Thompson Multitaper Method';
      
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      this.SpecifyDataWindowAs = s.SpecifyDataWindowAs;
      this.CombineMethod = s.CombineMethod;
      
      if strcmpi(this.SpecifyDataWindowAs,'TimeBW'),
        set(this, 'TimeBW', s.TimeBW);
      else
        set(this, ...
          'DPSS',           s.DPSS, ...
          'Concentrations', s.Concentrations);
      end
      
      
    end
    
    function [Pxx,W] = thispsd(this,x,opts)
      %THISPSD   Power spectral density (PSD) via MTM.
      %
      % OPTS = {NFFT, Fs(?), and SpectrumType}.
      %
      % This is a private method.
      
      narginchk(2,3);
      
      % Convert CombineMethod enum type to strings accepted by the function.
      combineMethod = getcombinemethodstr(this);
      opts = {opts{:},combineMethod};
      
      if strcmpi(this.SpecifyDataWindowAs,'TimeBW'),
        [Pxx, W] = pmtm(x,this.TimeBW,opts{:});
      else
        validatesizes(this,x); % Validate the size of E and V
        [Pxx, W] = pmtm(x,this.DPSS,this.Concentrations,opts{:});
      end
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
    end
    
  end  %% public methods
  
end  % classdef

function mode = setdatawindowspecmode(this,mode)
% Set function for the SpecifyDataWindowAs property.

propStr1 = 'TimeBW';
propStr2 = 'DPSS';
propStr3 = 'Concentrations';

if ~isprop(this,propStr1),
  % If it doesn't exist create it.
  dp = this.addprop(propStr1);
  this.(propStr1) = 4;
  dp.NonCopyable = false;
  % p1 = schema.prop(this,propStr1,'double');
  % set(this,propStr1,4);
else
  p1 = this.findprop(propStr1);
end

if ~isprop(this,propStr2),
  % If it doesn't exist create it.
  p2 = this.addprop(propStr2);
  p3 = this.addprop(propStr3);
  
  p2.NonCopyable = false;
  p3.NonCopyable = false;
  
  % Calculate default window tapers and their concentrations.
  nfft = 256;  % Arbitrary default value.
  [E,V] = dpss(nfft,this.TimeBW);
  set(this,'DPSS',E,'Concentrations',V);
  
else
  p2 = this.findprop(propStr2);
  p3 = this.findprop(propStr3);
end

% Enable/disable properties.
if strcmpi(mode,'TimeBW'),
  prop1State = 'on';
  prop2State = 'off';
  prop3State = 'off';
else
  prop1State = 'off';
  prop2State = 'on';
  prop3State = 'on';
  
end
enabdynprop(this,propStr1,prop1State);
enabdynprop(this,propStr2,prop2State);
enabdynprop(this,propStr3,prop3State);
end  % setdatawindowspecmode


%--------------------------------------------------------------------------
function validatesizes(this,x)
% Return error if size mismatch is found.

if size(this.DPSS,2) ~= length(this.Concentrations),
  error(message('signal:spectrum:mtm:thispsd:invalidNumCols', 'Concentration'));
end

if size(this.DPSS,1) ~= length(x),
  error(message('signal:spectrum:mtm:thispsd:invalidNumRows'));
end

end


% [EOF]
