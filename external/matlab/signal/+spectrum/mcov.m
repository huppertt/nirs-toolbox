classdef mcov < spectrum.abstractar
  %MCOV   Modified covariance power spectral density (PSD) estimator.
  %
  %   SPECTRUM.MCOV is not recommended.  Use <a href="matlab:help pmcov">pmcov</a> instead.
  %
  %   H = SPECTRUM.MCOV returns a modified covariance (mcov) PSD estimator in
  %   H.
  %
  %   H = SPECTRUM.MCOV(ORDER) returns a modified covariance spectral
  %   estimator with the order of the autoregressive (AR) model set to the
  %   numeric value specified by ORDER.
  %
  %   MCOV PSD estimators can be passed to the following functions along with
  %   the data to perform that function:
  %       <a href="matlab:help spectrum/psd">psd</a>     - calculates the PSD
  %       <a href="matlab:help spectrum/psdopts">psdopts</a> - returns options to calculate the PSD
  %
  %   EXAMPLE: Spectral analysis of a 4th order autoregressive (AR) process.
  %      x = randn(100,1);
  %      y = filter(1,[1 1/2 1/3 1/4 1/5],x);
  %      h = spectrum.mcov(4);                 % Create an mcov estimator.
  %      psd(h,y,'Fs',1000);                   % Calculate and plot the PSD.
  
  %   Author(s): P. Pacheco
  
  
  
  methods  % constructor block
    function this = mcov(varargin)
      
      
      narginchk(0,2);
      
      % Set the properties of the object.
      % this = spectrum.mcov;
      this.EstimationMethod = 'Modified Covariance';
      initialize(this,varargin{:});  % Set Nfft and Order
      
      
    end  % mcov
    
    function [Pxx,W] = thispsd(this,x,opts)
      %THISPSD Calculate the power spectral density (PSD) via Modified Covariance Method.
      %
      % This is a private method.
      
      narginchk(2,3);
      
      [Pxx, W] = pmcov(x,...
        this.Order,...
        opts{:});  % NFFT, Fs(?), and SpectrumType
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
    end
    
  end  %% public methods
  
end  % classdef

