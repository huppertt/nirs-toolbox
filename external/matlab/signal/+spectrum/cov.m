classdef cov < spectrum.abstractar
  %COV   Covariance power spectral density (PSD) estimator.
  %
  %   SPECTRUM.COV is not recommended.  Use <a href="matlab:help pcov">pcov</a> instead.
  %
  %   H = SPECTRUM.COV returns a covariance (cov) PSD estimator in H.
  %
  %   H = SPECTRUM.COV(ORDER) returns a covariance spectral estimator with
  %   the order of the autoregressive (AR) model set to the numeric value
  %   specified by ORDER.
  %
  %   COV PSD estimators can be passed to the following functions along with
  %   the data to perform that function:
  %       <a href="matlab:help spectrum/psd">psd</a>     - calculates the PSD
  %       <a href="matlab:help spectrum/psdopts">psdopts</a> - returns options to calculate the PSD
  %
  %   EXAMPLE: Spectral analysis of a 4th order autoregressive (AR) process.
  %      x = randn(100,1);
  %      y = filter(1,[1 1/2 1/3 1/4 1/5],x);
  %      h = spectrum.cov(4);                 % Create a cov estimator.
  %      psd(h,y,'Fs',1000);                  % Calculate and plot the PSD.
  
  
  methods  % constructor block
    function this = cov(varargin)
     
      narginchk(0,2);
      
      % Set the properties of the object.
      % this = spectrum.cov;
      this.EstimationMethod = 'Covariance';
      initialize(this,varargin{:});  % Sets NFFT and ORDER
      
      
    end  % cov
    
    
    function [Pxx,W] = thispsd(this,x,opts)
      %THISPSD Calculate the power spectral density (PSD) via Covariance Method.
      %
      % This is a private method.

      narginchk(2,3);
      
      [Pxx, W] = pcov(x,...
        this.Order,...
        opts{:});  % NFFT, Fs(?), and SpectrumType
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
    end
    
  end  %% public methods
  
end  % classdef

