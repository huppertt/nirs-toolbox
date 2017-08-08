classdef eigenvector < spectrum.music
  %EIGENVECTOR   Eigenvector pseudospectrum estimator.
  %
  %   SPECTRUM.EIGENVECTOR is not recommended.  Use <a href="matlab:help peig">peig</a> and <a href="matlab:help rooteig">rooteig</a> instead.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS) returns an eigenvector
  %   pseudospectrum estimator in H with the number of complex sinusoids set
  %   to the numeric value specified by NSINUSOIDS.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS,SEGMENTLENGTH) returns an
  %   eigenvector pseudospectrum estimator with the number of samples in each
  %   segment set to the value specified by SEGMENTLENGTH.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT) sets
  %   the percentage of overlap between segments to the value specified by
  %   OVERLAPPERCENT.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,...
  %   WINNAME) specifies the window as a string. Use set(H,'WindowName') to
  %   get a list of valid <a href="matlab:set(spectrum.eigenvector,'WindowName')">windows</a>.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,...
  %   {WINNAME,WINPARAMETER}) specifies the window in WINNAME and the window
  %   parameter value in WINPARAMETER both in a cell array.
  %
  %   NOTE: Depending on the window specified by WINNAME a window parameter
  %   property will be dynamically added to the eigenvector estimator H. Type
  %   "help <WINNAME>" for more details.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,...
  %   WINNAME,THRESHOLD) specifies the THRESHOLD as the cutoff for the signal
  %   and noise subspace separation.
  %
  %   H = SPECTRUM.EIGENVECTOR(NSINUSOIDS,SEGMENTLENGTH,OVERLAPPERCENT,...
  %   WINNAME,THRESHOLD,INPUTTYPE) specifies the type of input the
  %   eigenvector spectral estimator accepts. INPUTTYPE can be one of the
  %   following strings:
  %       'Vector'  (default)
  %       'DataMatrix'
  %       'CorrelationMatrix'
  %
  %   Eigenvector pseudospectrum estimators can be passed to the following
  %   functions along with the data to perform that function:
  %       <a href="matlab:help spectrum/powerest">powerest</a>           - computes the powers and frequencies of sinusoids
  %       <a href="matlab:help spectrum/pseudospectrum">pseudospectrum</a>     - calculates the pseudospectrum
  %       <a href="matlab:help spectrum/pseudospectrumopts.">pseudospectrumopts</a> - returns options to calculate the pseudospectrum
  %
  %   EXAMPLE: Spectral analysis of a signal containing complex sinusoids and
  %            % noise.
  %            n = 0:99;
  %            s = exp(1i*pi/2*n)+2*exp(1i*pi/4*n)+exp(1i*pi/3*n)+randn(1,100);
  %            h = spectrum.eigenvector(3);   % Create an eigenvector spectral estimator.
  %            pseudospectrum(h,s);           % Calculate and plot the pseudospectrum.
  
  %   Author(s): P. Pacheco
  
  
  
  methods  % constructor block
    function this = eigenvector(varargin)
      
      
      narginchk(0,7);
      
      % Set the properties of the object.
      % this = spectrum.eigenvector;
      this.EstimationMethod = 'Eigenvector';
      initialize(this,varargin{:});
      
      
    end  % eigenvector
    
    
    function varargout = powerest(this,x,Fs)
      %POWEREST   Computes the powers and frequencies of sinusoids.
      %   POW = POWEREST(H,X) returns the vector POW containing the estimates
      %   of the powers of the complex sinusoids contained in the data
      %   represented by X.  H must be a EIGENVECTOR object.
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
      %      H = spectrum.eigenvector(3);
      %      [P,W] = powerest(H,s);
      %
      %   See also ROOTMUSIC, PMUSIC, PEIG, PMTM, PBURG, PWELCH and CORRMTX.
      
      
      narginchk(2,3);
      
      if nargin < 3,
        Fs = 1;
      end
      
      P = [this.NSinusoids, this.SubspaceThreshold];
      
      if strcmpi(this.InputType,'CorrelationMatrix'),
        [w,pow] = rooteig(x,P,'corr',Fs);
      else
        [w,pow] = rooteig(x,P,Fs);
      end
      
      % If input is single, W will also be single so cast it to double
      w = double(w);
      varargout = {pow,w};
      
    end
    
    
    function [Sxx,W] = thispseudospectrum(this,x,opts,P)
      %THISPSEUDOSPECTRUM   Calculate the pseudospectrum via Eigenvector analysis.
      %
      % This is a private method.
      
      
      narginchk(2,4);
      
      if strcmpi(this.InputType,'CorrelationMatrix'),
        [Sxx, W] = peig(x,P,'corr',opts{:});
      else
        [Sxx, W] = peig(x,P,opts{:});
      end
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
    end
        
  end  %% public methods
  
end  % classdef

