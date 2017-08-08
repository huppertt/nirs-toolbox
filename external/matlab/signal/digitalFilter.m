classdef digitalFilter < dynamicprops & matlab.mixin.CustomDisplay
%digitalFilter Digital filter class
%   Design a digital filter D using the designfilt function and then apply
%   the filter to your data using the filter function.
%
%   D = <a href="matlab:help designfilt">designfilt</a>(RESP,P1,V1,P2,V2,...) designs a digital filter, D, with
%   response RESP and specifications P1,V1, ... PN,VN.
%
%   Y = <a href="matlab:help digitalFilter/filter">filter</a>(D,X) filters input data, X, with the digital filter, D,
%   and produces output data Y.
%
%   <a href="matlab:help fvtool">fvtool</a>(D) opens the filter visualization tool to analyze the response
%   of digital filter D.
%
%   designfilt(D) launches a Filter Design Assistant so that you can edit
%   the digital filter D.
%
%   There are several functions available to analyze your digital filter
%   and filter data. An exhaustive list is available in the documentation.
%   The most commonly used functions are listed below:
%
%   <a href="matlab:help digitalFilter/freqz">freqz</a>       - Frequency response
%   <a href="matlab:help digitalFilter/phasez">phasez</a>      - Phase response 
%   <a href="matlab:help digitalFilter/zerophase">zerophase</a>   - Zero-phase response
%   <a href="matlab:help digitalFilter/grpdelay">grpdelay</a>    - Group delay response
%   <a href="matlab:help digitalFilter/phasedelay">phasedelay</a>  - Phase delay response
%   <a href="matlab:help digitalFilter/impz">impz</a>        - Impulse response
%   <a href="matlab:help digitalFilter/stepz">stepz</a>       - Step response
%   <a href="matlab:help digitalFilter/zplane">zplane</a>      - Z-plane zero-pole plot
%   <a href="matlab:help digitalFilter/tf">tf</a>          - Convert digital filter to transfer function
%   <a href="matlab:help digitalFilter/zpk">zpk</a>         - Convert digital filter to zero-pole-gain representation
%   <a href="matlab:help digitalFilter/info">info</a>        - Information about digital filter
%   <a href="matlab:help digitalFilter/filt2block">filt2block</a>  - Generate Simulink filter block
%   <a href="matlab:help digitalFilter/single">single</a>      - Cast filter coefficients to single precision
%   <a href="matlab:help digitalFilter/double">double</a>      - Cast filter coefficients to double precision
%   <a href="matlab:help digitalFilter/issingle">issingle</a>    - True if digital filter contains single precision coefficients
%   <a href="matlab:help digitalFilter/isdouble">isdouble</a>    - True if digital filter contains double precision coefficients
%   <a href="matlab:help digitalFilter/filtfilt">filtfilt</a>    - Zero-phase forward and reverse filtering
%   <a href="matlab:help digitalFilter/fftfilt">fftfilt</a>     - Filter data with overlap-add method using FFT
%
%   % Example 1:
%   %   Design a lowpass FIR filter with passband frequency of 300 Hz,
%   %   stopband frequency of 350 Hz, passband ripple of 0.5 dB, and
%   %   stopband attenuation of 65 dB. The sample rate is 1 KHz.
%   %   Apply the filter to a vector of random data.
%
%   lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 300,...
%            'StopbandFrequency', 350, 'PassbandRipple', 0.5, ...
%            'StopbandAttenuation', 65, 'SampleRate', 1e3);
%
%   data = randn(1000,1);
%   y = filter(lpFilt,data);
%
%   % Example 2:
%   %   Design a highpass IIR filter with order 8, passband frequency of
%   %   75 KHz, and a passband ripple of 0.2 dB. Sample rate is 200 KHz.
%   %   Apply the filter to a vector of random data.
%
%   hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   data = randn(1000,1);
%   y = filter(hpFilt,data);
%
%   % Example 3:
%   %   Launch the Filter Design Assistant so that you can edit the filter 
%   %   designed in the example above.
%
%   designfilt(hpFilt)
%
%   See also <a href="matlab:help designfilt">designfilt</a>, <a href="matlab:help digitalFilter/filter">filter</a>.

%#ok<*MANU>

properties (Access = private)
  MetaData = [];
  DfiltObject = [];
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
methods (Hidden)
  function obj = digitalFilter(varargin)
    %figitalFilter Constructor method
    % Inputs:
    %         DFILT object
            
    obj.DfiltObject = varargin{1};
    obj.DfiltObject.FromDesignfilt = true;
    
    fm = getfmethod(obj.DfiltObject);
    fm.FromDesignfilt = true;
    setfmethod(obj.DfiltObject,fm);
        
    if nargin == 2
      obj.MetaData = varargin{2};
    end
    
    % Add dynamic properties based on design specifications
    fdes = getfdesign(obj.DfiltObject);
    
    [pnames,pvalues] = getSpecProperties(obj);
        
    % Design method, sample rate, and coefficients
    if fdes.NormalizedFrequency
      Fs = 2;
    else
      Fs = fdes.Fs;
    end
    
    [coeffPropNames, coeffPropValues] = getCoefficientPropNamesValues(obj);
    
    dm = signal.internal.DesignfiltProcessCheck.convertdesignmethodnames(...
      getdesignmethod(obj.DfiltObject),isfir(obj),'fdesignToDesignfilt');
    
    resp = getRespType(obj);
         
    pnames = [pnames; {'DesignMethod'; 'SampleRate'; ...
      'FrequencyResponse'; 'ImpulseResponse'}; coeffPropNames];
            
    pvalues = [pvalues; {dm; Fs; resp(1:end-3); resp(end-2:end);...
      coeffPropValues}];
    
    for idx = 1:numel(pnames)
      name = pnames{idx};
      obj.addprop(name); 
      % Set the properties before adding the set method since this will
      % make the properties read-only.
      obj.(name) = pvalues{idx};
      mobj = obj.findprop(name);
      mobj.SetMethod=@setProp;              
    end
  end
  
  %-----------------------------------------------------------------------
  function d = todfilt(obj)
    %TODFILT Convert digitalFilter object to a DFILT object
    d = copy(obj.DfiltObject);
  end  
end
%--------------------------------------------------------------------------
%  Public methods
%--------------------------------------------------------------------------
methods  
  %------------------------------------------------------------------------
  % Analysis methods
  %------------------------------------------------------------------------
  function [H, W] = freqz(obj,varargin)
  %FREQZ  Frequency response of digital filter
  %   [H,W] = freqz(D,N) returns the N-point complex frequency response
  %   vector H and the N-point frequency vector W in radians/sample of the
  %   digital filter, D. The frequency response is evaluated at N points
  %   equally spaced around the upper half of the unit circle. If N is not
  %   specified, it defaults to 8192.
  %
  %   freqz(D) with no output argument launches FVTool and shows the
  %   magnitude and phase responses of filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/freqz">signal/freqz</a>.
  %
  %   See also FVTOOL.
    try        
      if nargout
        c = getCoeffsForAnalysis(obj);
        [H,W] = freqz(c{:},varargin{:});
      else
        freqz(obj.DfiltObject,varargin{:});
      end
    catch ME
        throwAsCaller(ME);
    end
  end  
  %------------------------------------------------------------------------
  function [Phi, W] = phasez(obj,varargin)
  %PHASEZ Phase response of digital filter 
  %   [Phi,W] = phasez(D,N) returns vectors Phi and W containing the phase
  %   response of the digital filter, D, and the frequencies (in radians)
  %   at which it is evaluated. The phase response is evaluated at N points
  %   equally spaced around the upper half of the unit circle. If you do
  %   not specify N, it defaults to 8192.
  %
  %   phasez(D) with no output argument launches FVTool and shows the
  %   phase response of filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/phasez">signal/phasez</a>.
  %
  %   See also FVTOOL.

    try      
      if nargout
        c = getCoeffsForAnalysis(obj);
        [Phi,W] = phasez(c{:},varargin{:});
      else
        phasez(obj.DfiltObject,varargin{:});
      end
    catch ME
        throwAsCaller(ME);
    end    
  end
  %------------------------------------------------------------------------
  function [H, W] = zerophase(obj,varargin)
  %ZEROPHASE Zero-phase response of digital filter
  %   [H,W] = zerophase(D,N) returns length N vectors H and W containing
  %   the zero-phase response of the digital filter, D, and the frequencies
  %   (in radians) at which it is evaluated. The zero-phase response is
  %   evaluated at N points equally spaced around the upper half of the
  %   unit circle. If you don't specify N, it defaults to 8192.
  %
  %   zerophase(D) with no output argument launches FVTool and shows the
  %   zero-phase response of filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/zerophase">signal/zerophase</a>.
  %
  %   See also FVTOOL.
  
    try      
      if nargout
        c = getCoeffsForAnalysis(obj);
        [H,W] = zerophase(c{:},varargin{:});
      else
        zerophase(obj.DfiltObject,varargin{:});
      end
    catch ME
        throwAsCaller(ME);
    end    
  end
  %------------------------------------------------------------------------
  function [Gd, W] = grpdelay(obj,varargin)
  %GRPDELAY Group delay response of digital filter
  %   [Gd,W] = grpdelay(D,N) returns length N vectors Gd and W containing
  %   the group delay of the digital filter, D, and the frequencies (in
  %   radians) at which it is evaluated. The group delay is defined as
  %   -d{angle(w)}/dw. The frequency response is evaluated at N points
  %   equally spaced around the upper half of the unit circle. If you do
  %   not specify N, it defaults to 8192.
  %
  %   grpdelay(D) with no output argument launches FVTool and shows the
  %   group delay response of the filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/grpdelay">signal/grpdelay</a>.
  %
  %   See also FVTOOL.

    try
     if nargout
       c = getCoeffsForAnalysis(obj);
       [Gd,W] = grpdelay(c{:},varargin{:});
     else
       grpdelay(obj.DfiltObject,varargin{:});
      end   
    catch ME
       throwAsCaller(ME);
    end    
  end
  %------------------------------------------------------------------------
  function [Phi, W] = phasedelay(obj,varargin)
  %PHASEDELAY Phase delay response of digital filter
  %   [PHI,W] = phasedelay(D,N) returns the N-point phase delay response
  %   vector PHI and the N-point frequency vector W in radians/sample of
  %   the digital filter, D. The phase response is evaluated at N points
  %   equally spaced around the upper half of the unit circle. If N is not
  %   specified, it defaults to 8192.
  %
  %   phasedelay(D) with no output argument launches FVTool and shows the
  %   phase delay response of filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/phasedelay">signal/phasedelay</a>.
  %
  %   See also FVTOOL.

    try
     if nargout
       c = getCoeffsForAnalysis(obj);
       [Phi,W] = phasedelay(c{:},varargin{:});
     else
        phasedelay(obj.DfiltObject,varargin{:});
      end
    catch ME
        throwAsCaller(ME);
    end    
  end
  %------------------------------------------------------------------------
  function [H, T] = impz(obj,varargin)
  %IMPZ Impulse response of digital filter
  %   [H,T] = impz(D) computes the impulse response of the digital filter,
  %   D, and returns the response in column vector H, and a vector of times
  %   (or sample intervals) in T (T = [0 1 2...]'). The number of samples
  %   in the response is chosen automatically.
  %
  %   impz(D) with no output argument launches FVTool and shows the impulse
  %   response of filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/impz">signal/impz</a>.
  %
  %   See also FVTOOL.
  
    try
      if nargout
       c = getCoeffsForAnalysis(obj);        
       [H,T] = impz(c{:},varargin{:});
      else
        impz(obj.DfiltObject,varargin{:});
      end  
    catch ME
        throwAsCaller(ME);
    end    
  end
  %----------------------------------------------------------------------    
  function [H, W] = freqrespest(obj,varargin)
  %FREQRESPEST  Frequency response estimate of digital filter
  %   [H,W] = freqrespest(D,L) estimates the frequency response estimate,
  %   H, of the digital filter, D, by running input data, made up from
  %   sinusoids with uniformly distributed random frequencies, through the
  %   filter and forming the ratio between output data and input data. W is
  %   the vector of frequencies at which the response is estimated.
  %
  %   L is the number of trials used to compute the estimate. If not
  %   specified, L defaults to 10. In general, the more trials used, the
  %   more accurate estimate obtained at the expense of a longer time
  %   needed to compute the estimate.
  %
  %   [H,W] = freqrespest(D,L,P1,V1,P2,V2,...) specifies optional
  %   parameters via name-value pairs. Valid pairs are:
  %
  %         Parameter           Default        Description/Valid values
  %   ---------------------  -----------  ----------------------------------
  %    'NFFT'                 512          Number of FFT points.
  %    'NormalizedFrequency'  true         {true,false}
  %    'Fs'                   N/A          Sample rate. Only applies when
  %                                        'NormalizedFrequency' is false.
  %    'FrequencyRange'       'onesided'   {'onesided','twosided','centered'}
  %
  %   freqrespest(D,...) with no output argument launches FVTool and shows
  %   the magnitude response estimate of the digital filter, D.            
  %
  %   See also FVTOOL.

  
  % See if we have a FrequencyRange input, parse it and translate it to
  % DFILT syntax
  idx = find(strcmpi(varargin,'FrequencyRange')==true);
  
  if ~isempty(idx)
    if numel(varargin) > idx
      freqRange = validatestring(varargin{idx+1},...
        {'onesided','twosided','centered'},'digitalFilter/freqrespest',...
        'FrequencyRange');
      
      varargin(idx:idx+1) = [];
    else
      error(message('signal:digitalFilter:MustSpecifyFreqRange'))
    end
    
    switch freqRange
      case 'centered'
        varargin = [varargin {'CenterDC',true}];
      case 'onesided'
        varargin = [varargin {'SpectrumRange','Half'}];
      case 'twosided'
        varargin = [varargin {'SpectrumRange','Whole'}];
    end
  end
      
   try
      if nargout
        [H,W] = freqrespest(obj.DfiltObject,varargin{:});
      else
        freqrespest(obj.DfiltObject,varargin{:});
      end  
    catch ME
        throwAsCaller(ME);
   end    
  end  
  %----------------------------------------------------------------------    
  function [H,W] = noisepsd(obj,varargin)       
  %NOISEPSD Power spectral density of filter output due to roundoff noise
  %   [PSD,W] = noisepsd(D,L) computes the power spectral density, PSD, at
  %   the output of the digital filter, D, due to roundoff noise produced
  %   by quantization errors within the filter. The PSD is computed by
  %   averaging L trials. The larger the number of trials, the better the
  %   estimate (at the expense of longer computation time). If L is not
  %   specified, it defaults to 10 trials. W is the vector of frequencies
  %   at which the PSD is estimated.
  %
  %   You can use the bandpower function to compute the average power of
  %   the output noise. 
  %  
  %   [PSD,W] = noisepsd(H,L,P1,V1,P2,V2,...) specifies optional parameters
  %   via name-value pairs. Valid pairs are:
  %
  %        Parameter           Default        Description/Valid values
  %   ---------------------  -----------  ----------------------------------
  %   'NFFT'                 512          Number of FFT points.
  %   'NormalizedFrequency'  true         {true,false}
  %   'Fs'                   N/A          Sample rate. Only applies when
  %                                       'NormalizedFrequency' is false.
  %   'FrequencyRange'       'onesided'   {'onesided','twosided','centered'}
  %
  %   noisepsd(D,...) with no output argument launches FVTool and shows
  %   the noise PSD estimate of the digital filter, D.
  %
  %   See also FVTOOL, BANDPOWER.

    % See if we have a FrequencyRange input, parse it and translate it to
    % DFILT syntax
    idx = find(strcmpi(varargin,'FrequencyRange')==true);
    
    if ~isempty(idx)
      if numel(varargin) > idx
        freqRange = validatestring(varargin{idx+1},...
          {'onesided','twosided','centered'},'digitalFilter/noisepsd',...
          'FrequencyRange');
        
        varargin(idx:idx+1) = [];
      else
        error(message('signal:digitalFilter:MustSpecifyFreqRange'))
      end
      
      switch freqRange
        case 'centered'
          varargin = [varargin {'CenterDC',true}];
        case 'onesided'
          varargin = [varargin {'SpectrumType','Onesided'}];
        case 'twosided'
          varargin = [varargin {'SpectrumType','Twosided'}];
      end
    end
    
    try
      if nargout
        Hpsd = noisepsd(obj.DfiltObject,varargin{:});
        H = Hpsd.Data;
        W = Hpsd.Frequencies;
      else
        noisepsd(obj.DfiltObject,varargin{:});
      end  
    catch ME
        throwAsCaller(ME);
   end    
  end  
  %------------------------------------------------------------------------
  function L = impzlength(obj,varargin)
  %IMPZLENGTH Length of impulse response of digital filter
  %   L = impzlength(D) returns the length, L, of the impulse response of
  %   the digital filter, D.
  %
  %   L = impzlength(D,TOL) specifies the tolerance, TOL, to increase or
  %   decrease the length accuracy. By default, TOL = 5e-5.
  
    try
      c = getCoeffsForAnalysis(obj);
      L = impzlength(c{:},varargin{:});
    catch ME
        throwAsCaller(ME);
    end      
  end
  
  %------------------------------------------------------------------------
  function [H, T] = stepz(obj,varargin)
  %STEPZ Step response of digital filter
  %   [H,T] = stepz(D) returns the step response H of the digital filter,
  %   D. The length of column vector H is computed using IMPZLENGTH. The
  %   vector of time samples at which H is evaluated is returned in vector
  %   T.
  %
  %   stepz(D) with no output argument launches FVTool and shows the step
  %   response of filter D.
  %
  %   For additional parameters, see <a href="matlab:help signal/stepz">signal/stepz</a>.
  %
  %   See also FVTOOL.

    try
      if nargout
        c = getCoeffsForAnalysis(obj);
        [H,T] = stepz(c{:},varargin{:});
      else
        stepz(obj.DfiltObject,varargin{:});
      end
    catch ME
        throwAsCaller(ME);
    end    
  end
  
  %------------------------------------------------------------------------
  function [HZ,HP,Hl] = zplane(obj,varargin)  
  %ZPLANE Z-plane zero-pole plot for digital filter
  %   zplane(D) launches FVTool and plots the zeros and poles of the
  %   digital filter, D, with the unit circle for reference. Each zero is
  %   represented with a 'o' and each pole with a 'x' on the plot. Multiple
  %   zeros and poles are indicated by the multiplicity number shown to the
  %   upper right of the zero or pole.
  %
  %   [HZ,HP,Hl] = zplane(D) returns vectors of handles to the lines and
  %   text objects generated. HZ is a vector of handles to the zeros lines,
  %   HP is a vector of handles to the poles lines, and Hl is a vector of
  %   handles to the axes / unit circle line and to text objects which are
  %   present when there are multiple zeros or poles. In case there are no
  %   zeros or no poles, HZ or HP is set to the empty matrix [].
  %
  %   For additional parameters, see <a href="matlab:help signal/zplane">signal/zplane</a>.
  %
  %   See also FVTOOL.
  
    try
      if nargout > 0
        [HZ,HP,Hl] = zplane(obj.DfiltObject,varargin{:});
      else
        zplane(obj.DfiltObject,varargin{:})
      end
    catch ME
        throwAsCaller(ME);
    end    
  end  
  %------------------------------------------------------------------------
  function N = filtord(obj,varargin)
  %FILTORD  Order of digital filter
  %   N = filtord(D) returns the order, N, of the digital filter, D.

    try
      c = getCoeffsForAnalysis(obj);
      N = filtord(c{:},varargin{:});
    catch ME
        throwAsCaller(ME);
    end      
  end  
  %------------------------------------------------------------------------
  function t = firtype(obj,varargin)      
  %FIRTYPE Determine the type (1-4) of a linear phase FIR digital filter
  %   T = firtype(D) determines the type (1 through 4) of the FIR digital
  %   filter, D. The filter must be real and have linear phase.
  %
  %   Types 1 through 4 are defined as follows:
  %       
  %   Type 1: Even order symmetric coefficients.
  %   Type 2: Odd order symmetric coefficients.
  %   Type 3: Even order antisymmetric coefficients.
  %   Type 4: Odd order antisymmetric coefficients.

    if ~isfir(obj)
      error(message('signal:digitalFilter:MethodOnlyValidForFIR','firtype'));
    end

    try      
      t = firtype(obj.Coefficients,varargin{:});        
    catch ME
      throwAsCaller(ME);
    end      
  end   
  %------------------------------------------------------------------------
  function [b,a] = tf(obj,varargin)      
  %TF Convert digital filter to transfer function
  %   [NUM,DEN] = tf(D) converts digital filter, D, to numerator and
  %   denominator vectors.
  
    try      
      [b,a] = tf(obj.DfiltObject,varargin{:});                
    catch ME
      throwAsCaller(ME);
    end  
  end
  %------------------------------------------------------------------------
  function [z,p,k] = zpk(obj,varargin)      
  %ZPK Convert digital filter to zero-pole-gain representation
  %   [Z,P,K] = zpk(D) returns the zeros, poles, and gain corresponding to
  %   the digital filter, D, in vectors Z, P, and scalar K respectively.
    
    try
      [z,p,k] = zpk(obj.DfiltObject,varargin{:});   
    catch ME
      throwAsCaller(ME);
    end
  end  
  %------------------------------------------------------------------------
  function [A,B,C,D] = ss(obj,varargin)    
  %SS Convert digital filter to state-space representation
  %   [A,B,C,D] = ss(DF) converts digital filter, DF, to state-space
  %   representation given by:
  %   x(k+1) = A*x(k) + B*u(k)
  %   y(k)   = C*x(k) + D*u(k)
  %   where x is the state vector, u is the input vector, and y is the 
  %   output vector. 
          
    try
      [A,B,C,D] = ss(obj.DfiltObject,varargin{:});             
    catch ME
      throwAsCaller(ME);
    end
  end  
  %------------------------------------------------------------------------
  function flag = isallpass(obj,varargin)
  %ISALLPASS True for all-pass digital filter
  %   FLAG = isallpass(D) returns true if digital filter, D, is all-pass.
  % 
  %   FLAG = isallpass(D,TOL) uses tolerance TOL to determine when two
  %   numbers are close enough to be considered equal.
       
    try
      c = getCoeffsForAnalysis(obj);
      flag = isallpass(c{:},varargin{:});                 
    catch ME
      throwAsCaller(ME);
    end
  end  
  %------------------------------------------------------------------------
  function flag = islinphase(obj,varargin)
  %ISLINPHASE True for linear phase digital filter
  %   FLAG = islinphase(D) returns true if digital filter, D, is linear
  %   phase.
  %
  %   FLAG = islinphase(D,TOL) uses tolerance TOL to determine when two
  %   numbers are close enough to be considered equal.
    
    try
      c = getCoeffsForAnalysis(obj);
      flag = islinphase(c{:},varargin{:});                         
    catch ME
      throwAsCaller(ME);
    end
  end        
  %------------------------------------------------------------------------
  function flag = ismaxphase(obj,varargin)
  %ISMAXPHASE True for maximum phase digital filter
  %   FLAG = ismaxphase(D) returns true if digital filter, D, is maximum
  %   phase.
  %
  %   FLAG = ismaxphase(D,TOL) uses tolerance TOL to determine when two
  %   numbers are close enough to be considered equal.
    
    try
      c = getCoeffsForAnalysis(obj);
      flag = ismaxphase(c{:},varargin{:});                         
    catch ME
      throwAsCaller(ME);
    end
  end                    
  %------------------------------------------------------------------------
  function flag = isminphase(obj,varargin)
  %ISMINPHASE  True for minimum phase digital filter
  %   FLAG = isminphase(D) returns true if digital filter, D, is minimum
  %   phase.
  %
  %   FLAG = isminphase(H,TOL) uses tolerance TOL to determine when two
  %   numbers are close enough to be considered equal.
  
    try
      c = getCoeffsForAnalysis(obj);
      flag = isminphase(c{:},varargin{:});                         
    catch ME
      throwAsCaller(ME);
    end
  end                    
  %------------------------------------------------------------------------
  function flag = isstable(obj,varargin)
  %ISSTABLE True for stable digital filter
  %   FLAG = isstable(D) returns true if digital filter, D, is stable.
  
    try
      c = getCoeffsForAnalysis(obj);
      flag = isstable(c{:},varargin{:});                         
    catch ME
      throwAsCaller(ME);
    end
  end                      
  %------------------------------------------------------------------------
  function s = info(obj,varargin)
  %INFO Information about digital filter
  %   S = info(D) returns a string matrix with information about the
  %   digital filter, D.
  
    try
      s = info(obj.DfiltObject,varargin{:},'long');
    catch ME
      throwAsCaller(ME);
    end
  end
  %------------------------------------------------------------------------
  function flag = isfir(obj)
  %ISFIR  True if digital filter is FIR
  %   FLAG = isfir(D) returns true if filter, D, is FIR.
  
    resp = getRespType(obj);
    flag = strcmp(resp(end-2:end),'fir');
  end  
  %------------------------------------------------------------------------ 
  function filt2block(obj,varargin)
  %FILT2BLOCK Generate Simulink filter block
  %   FILT2BLOCK(D) generates a filter block with coefficients equal to
  %   those available in the digital filter, D. A Discrete FIR Filter block
  %   is generated  when D is an FIR filter. A Biquad Filter block is
  %   generated when D is an IIR filter (*).
  %
  %   FILT2BLOCK(D,'subsystem') generates a Simulink subsystem block that
  %   implements the digital filter, D, using sum, gain, and delay blocks.
  %
  %   FILT2BLOCK(D,...,'FilterStructure',STRUCT) for FIR filters, specify
  %   the filter structure, STRUCT, as one of 'directForm' |
  %   'directFormTransposed' | 'directFormSymmetric' |
  %   'directFormAntiSymmetric' | 'overlapAdd'. The default is
  %   'directForm'. 'overlapAdd' is only valid when you omit the
  %   'subsystem' flag (*). 
  %   For IIR filters, specify the filter structure, STRUCT, as one of
  %   'directForm1' | 'directForm1Transposed' | 'directForm2' |
  %   'directForm2Transposed'. The default is 'directForm2Transposed'.
  %
  %   (*) A DSP System Toolbox license is required for this syntax.
  %
  %   For additional parameters, see <a href="matlab:help signal/filt2block">signal/filt2block</a>.
  
    try
      if isfir(obj)
        filt2block(obj.Coefficients,varargin{:});
      else
        SOS = verifySOSMatrix(obj);
        filt2block(SOS,varargin{:});
      end
    catch ME
       throwAsCaller(ME);
    end
  end
  %------------------------------------------------------------------------
  function objOut = double(obj)
  %DOUBLE Cast coefficients of digital filter to double precision
  %   F2 = double(F1) casts coefficients in digital filter, F1, to double
  %   precision and returns a new digital filter, F2, that contains these
  %   coefficients.
  
    if isdouble(obj)
      objOut = obj;
      return;
    end
    
    % Going from single to double
    df = copy(obj.DfiltObject);
    
    % Cache metadata since it will be list when we write to the df
    % coefficients
    fdes = getfdesign(df);
    fmthd = getfmethod(df);    
    desmthd = getdesignmethod(df);
    
    df.Arithmetic = 'double';
    % Make sure we store coefficients in their quantized form. Remember
    % that DFILT keeps reference double precision coefficients internally
    % and quantized coefficients on the coefficients properties. The
    % quantized cofficients are always stored in double precision. 
    %
    % The reference coefficients should be the quantized coefficients. 
    if isfir(obj)
      df.Numerator = obj.DfiltObject.Numerator;        
    else
      df.sosMatrix = obj.DfiltObject.sosMatrix;
      df.ScaleValues = obj.DfiltObject.ScaleValues;      
    end
    
    % Restore metadata lost while setting the coefficients of df
    setfdesign(df,fdes)
    setfmethod(df,fmthd)
    setdesignmethod(df,desmthd)
    
    objOut = digitalFilter(df);
  end
  %------------------------------------------------------------------------
  function objOut = single(obj)
  %SINGLE Cast coefficients of digital filter to single precision
  %   F2 = single(F1) casts coefficients in digital filter, F1, to single
  %   precision and returns a new digital filter, F2, that contains these
  %   coefficients.
  
    if issingle(obj)
      objOut = obj;
      return;
    end
    
    % Going from double to single
    df = copy(obj.DfiltObject);
    df.Arithmetic = 'single';
    objOut = digitalFilter(df);
  end  
  %------------------------------------------------------------------------
  function flag = isdouble(obj)
  %ISDOUBLE True if digital filter coefficients are double precision
  %   FLAG = isdouble(D) returns true if the coefficients in digital
  %   filter, D, are double precision.
  
    flag = isa(obj.Coefficients,'double');
  end
  %------------------------------------------------------------------------
  function flag = issingle(obj)
  %ISSINGLE True if digital filter coefficients are single precision
  %   FLAG = issingle(D) returns true if the coefficients in digital
  %   filter, D, are single precision.
  
    flag = isa(obj.Coefficients,'single');
  end
  %------------------------------------------------------------------------
  function y = filter(obj,x)
  %FILTER Filter data with digital filter
  %   y = FILTER(D,x) filters input data x with digital filter D and
  %   produces output data Y. Input x can be a double or single precision
  %   vector or matrix with as many columns as input channels.
  %
  %   You design a <a href="matlab:help digitalFilter">digital filter</a> D using the <a href="matlab:help designfilt">designfilt</a> function and then 
  %   apply the filter to your data using the filter function.
  %
  %   % Example 1:
  %   %   Design a lowpass FIR filter with passband frequency of 300 Hz,
  %   %   stopband frequency of 350 Hz, passband ripple of 0.5 dB, and 
  %   %   stopband attenuation of 65 dB. The sample rate is 1 KHz. 
  %   %   Apply the filter to a vector of random data. 
  %   
  %   lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 300,...
  %            'StopbandFrequency', 350, 'PassbandRipple', 0.5, ...
  %            'StopbandAttenuation', 65, 'SampleRate', 1e3);
  %
  %   data = randn(1000,1);
  %   y = filter(lpFilt,data);
  % 
  %   % Example 2:
  %   %   Design a highpass IIR filter with order 25, passband frequency of 
  %   %   75 KHz, and a passband ripple of 0.2 dB. Sample rate is 200 KHz.
  %   %   Apply the filter to a vector of random data.
  %
  %   hpFilt = designfilt('highpassiir', 'FilterOrder', 25, ...
  %            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
  %            'SampleRate', 200e3);
  %
  %   data = randn(1000,1);
  %   y = filter(hpFilt,data);
  %
  %   See <a href="matlab:help designfilt">designfilt</a>, and <a href="matlab:help digitalFilter">digitalFilter</a>  help for more information.
    narginchk(2,2)
  
    if ~isfloat(x)
      error(message('signal:digitalFilter:InputMustBeFloat'));
    end
  
    if isfir(obj)
      y = filter(obj.Coefficients,1,x);
    else
      y = sosfilt(obj.Coefficients,x);
    end      
  end
  %------------------------------------------------------------------------
  function y = fftfilt(obj,x,varargin)
  %FFTFILT Filter data with overlap-add method using FFT
  %   Y = FFTFILT(D,X) filters X, with the FIR digital filter D using the
  %   overlap/add method, and internal parameters (FFT size and block
  %   length) that guarantee efficient execution.
  %   
  %   Y = FFTFILT(D,X,N) allows you to have some control over the internal
  %   parameters, by using an FFT of at least N points.
  %
  %   % Example:
  %   %   Design a lowpass FIR filter with order 350 and cutoff frequency
  %   %   of 150 Hz. The sample rate is 1.5 KHz. Filter a long vector of 
  %   %   data using the overlap-add method to increase speed.
  %
  %   D = designfilt('lowpassfir', 'FilterOrder', 350, ...
  %    'CutoffFrequency', 150, 'SampleRate', 1500);
  %
  %   data = randn(10e6,1);  
  %   y = fftfilt(D,data);
   
    narginchk(2,3)
  
    if ~isfir(obj)
      error(message('signal:digitalFilter:MethodOnlyValidForFIR','fftfilt'));
    end
    
    if issingle(obj)
      error(message('signal:digitalFilter:FilterMustBeDouble','fftfilt'));
    end
    
    try
      y = fftfilt(obj.Coefficients,x,varargin{:});
    catch ME
      throwAsCaller(ME);
    end
  end
  %------------------------------------------------------------------------
  function y = filtfilt(obj,x)
  %FILTFILT Zero-phase forward and reverse filtering
  %   Y = FILTFILT(D, X) filters the data in vector X with the digital
  %   filter D. The length of the input X must be more than three times the
  %   filter order, You can use filtord(D) to get the order of the digital
  %   filter D.
  %
  %   After filtering in the forward direction, the filtered sequence is
  %   then reversed and run back through the filter; Y is the time reverse
  %   of the output of the second filtering operation.  The result has
  %   precisely zero phase distortion, and magnitude modified by the square
  %   of the filter's magnitude response.
  %
  %   Note that FILTFILT should not be used when the intent of a filter is
  %   to modify signal phase, such as differentiators and Hilbert filters.
  %
  %   % Example:
  %   %   Design a highpass IIR digital filter 
  %   %   with order 4, passband frequency of 75 KHz, and a passband ripple 
  %   %   of 0.2 dB. Sample rate is 200 KHz. Apply zero-phase filtering to 
  %   %   a vector of data. 
  %  
  %   D = designfilt('highpassiir', 'FilterOrder', 4, ...
  %            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
  %            'SampleRate', 200e3);
  %
  %   x = rand(1000,1);
  %   y = filtfilt(D,x);
  
    narginchk(2,2)
    
    if issingle(obj)
      error(message('signal:digitalFilter:FilterMustBeDouble','filtfilt'));
    end
    
    try
      c = getCoeffsForAnalysis(obj);   
      if isfir(obj)
        y = filtfilt(c{:},x);
      else
        y = filtfilt(c{:},1,x);
      end
    catch ME
      throwAsCaller(ME);
    end    
  end

end %(end public methods)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
methods (Access = protected)
  %------------------------------------------------------------------------
  function s = getHeader(obj)
    %getHeader Get header for object display    
        
    helpStr = getHelpHyperlink('digitalFilter');  
    s = sprintf(' %s with properties:\n',helpStr);
  end
  %------------------------------------------------------------------------
  function s = getFooter(obj)
    %getFooter Get footer for object display      
        
    helpStr1 = getHelpHyperlink('fvtool');
    helpStr2 = getHelpHyperlink('designfilt');    
    if isMATLABInNormalMode()
      helpStr3 = sprintf('<a href="matlab:help(''%s'')">%s</a>','digitalFilter/filter','filter');
    else
      helpStr3 = 'filter';
    end
    
    str1 = getString(message('signal:digitalFilter:UseFvtool',helpStr1));      
    str2 = getString(message('signal:digitalFilter:UseDesignFilt',helpStr2));    
    str3 = getString(message('signal:digitalFilter:UseFilterFunction',helpStr3));
            
    s = sprintf(' %s\n %s\n %s\n',str1,str2,str3);
  end
  %------------------------------------------------------------------------
  function propgrp = getPropertyGroups(obj)
    %getPropertyGroups Group properties in order for object display
    
    % List 1 - coefficients
    propList1 = getCoefficientPropNamesValues(obj);
    
    % List 2 - specifications     
    propList2 = getSpecProperties(obj);        
    
    propList2 = [{'FrequencyResponse'; 'ImpulseResponse'; 'SampleRate'}; ...
      propList2; {'DesignMethod'}];
    
    propgrp1 = matlab.mixin.util.PropertyGroup(propList1);
    propgrp2 = matlab.mixin.util.PropertyGroup(propList2,'Specifications:');
    
    propgrp = [propgrp1 propgrp2];
  end
end %(end protected methods)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
methods (Access = private)
  function desObj = getDesignObj(obj)
    %getDesignObj Get design object
    desObj = getfdesign(obj.DfiltObject);
  end
  %------------------------------------------------------------------------
  function resp = getRespType(obj)
    %getRespType Get response type
    
    fdes = getDesignObj(obj);
    
    resp = class(fdes);
    idx = strfind(resp,'.');
    if isa(obj.DfiltObject,'dfilt.dffir')
      resp = [resp(idx+1:end) 'fir'];
    else
      resp = [resp(idx+1:end) 'iir'];
    end
  end
  %------------------------------------------------------------------------
  function [coeffPropNames, coeffPropValues] = getCoefficientPropNamesValues(obj)
    %getCoefficientPropNamesValues Get coefficient property names and values
    
    if isfir(obj)
      coeffPropValues = obj.DfiltObject.Numerator;
    else
      coeffPropValues = getSOSMatrix(obj);
    end
    
    if strcmpi(obj.DfiltObject.Arithmetic,'single')
      coeffPropValues = single(coeffPropValues);
    end    
    
    coeffPropNames = {'Coefficients'};            
  end  
  %------------------------------------------------------------------------
  function sos = getSOSMatrix(obj)
    %getSOSMatrix Get SOS matrix with embedded scale values
    
    if isfir(obj)
      error(message('signal:digitalFilter:MethodOnlyValidForSOS', ...
        'getSOSMatrix'))
    end
    
    sos = obj.DfiltObject.sosMatrix;
    g = obj.DfiltObject.ScaleValues;
    
    numSecs = size(sos,1);
    
    if length(g) == 1
      sos(1,1:3) = sos(1,1:3).*g;
    else
      if (length(g) == numSecs+1)
        lastGain = 1;
        if g(end) ~= 1
          % Spread the last non-unity gain across all sections
          lastGain = g(end)^(1/numSecs);
        end
        g = g(1:end-1)*lastGain;
      end
      g = g(:,[1 1 1]);
      sos(:,1:3) = sos(:,1:3).*g;
    end
  end  
  %------------------------------------------------------------------------
  function [designfiltProps,designfiltPropValues] = getSpecProperties(obj)
    %getSpecProperties Get spec property list from fdesign object

    fdesObj = getDesignObj(obj);
    [designfiltProps,designfiltPropValues] = ...
      signal.internal.DesignfiltProcessCheck.getSpecPropertiesFromFdesign(fdesObj);
  end  
  %------------------------------------------------------------------------    
  function SOS = verifySOSMatrix(obj)
    
    if isfir(obj)
      error(message('signal:digitalFilter:MethodOnlyValidForSOS', ...
        'verifySOSMatrix'))
    end
    
    SOS = obj.Coefficients;
    
    if size(obj.Coefficients,1) == 1
      SOS = [SOS ; [1 0 0 1 0 0]];      
    end    
  end
  %------------------------------------------------------------------------    
  function c = getCoeffsForAnalysis(obj)
  %getCoeffsForAnalysis Return a cell array with coefficients b,a or SOS
    if isfir(obj)
      c = {obj.Coefficients,1};
    else    
      SOS = obj.Coefficients;
    
      % Add a line with an all pass filter if filter is second order to
      % avoid ambiguities when calling analysis function.
      if size(obj.Coefficients,1) == 1
        SOS = [SOS ; [1 0 0 1 0 0]];      
      end 
      c = {SOS};
      
    end
  end
end %(end private methods)

end %(end classdef)

%--------------------------------------------------------------------------
function str = getHelpHyperlink(fcnName)
%getHelpHyperlink Get command line help hyperlink
  if isMATLABInNormalMode()
    str = sprintf('<a href="matlab:help(''%s'')">%s</a>',fcnName,fcnName);
  else
    str = fcnName;
  end
end
%--------------------------------------------------------------------------
function setProp(~,~)
%setProp Set method for all dynamic properties

% Properties are read-only. Error out and suggest the use of designfilt to
% design a filter with new specifications. 
  linkStr = getHelpHyperlink('designfilt');
  error(message('signal:digitalFilter:ReadOnlyProp',linkStr))
end
%------------------------------------------------------------------------
function flag = isMATLABInNormalMode()
% isMATLABInNormalMode True if we have java and desktop
flag = all([usejava('jvm') usejava('swing') usejava('awt') ...
  usejava('desktop')]);
end
