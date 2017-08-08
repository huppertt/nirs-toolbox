classdef burg < spectrum.abstractar 
  %   BURG   Burg power spectral density (PSD) estimator.
  %
  %   SPECTRUM.BURG is not recommended.  Use <a href="matlab:help pburg">pburg</a> instead.
  %
  %   H = SPECTRUM.BURG returns a Burg PSD estimator in H.
  %
  %   H = SPECTRUM.BURG(ORDER) returns a Burg spectral estimator with the
  %   order of the autoregressive (AR) model set to the numeric value
  %   specified by ORDER.
  %
  %   Burg PSD estimators can be passed to the following functions along with
  %   the data to perform that function:
  %       <a href="matlab:help spectrum/psd">psd</a>     - calculates the PSD
  %       <a href="matlab:help spectrum/psdopts">psdopts</a> - returns options to calculate the PSD
  %
  %   EXAMPLE: Spectral analysis of a 4th order autoregressive (AR) process.
  %      x = randn(100,1);
  %      y = filter(1,[1 1/2 1/3 1/4 1/5],x);
  %      h = spectrum.burg(4);                % Create a burg estimator.
  %      psd(h,y,'Fs',1000);                  % Calculate and plot the PSD.
  
  
  
  
  methods  % constructor block
    function this = burg(varargin)
     
      narginchk(0,2);
      
      % Set the properties of the object.
      % this = spectrum.burg;
      this.EstimationMethod = 'Burg';
      initialize(this,varargin{:});  % Set Nfft and Order
      
      
      
    end  % burg
    
    function [Pxx,W] = thispsd(this,x,opts)
      %THISPSD Calculate the power spectral density (PSD) via Burg.
      %
      % This is a private method.

      narginchk(2,3);
      
      [Pxx, W] = pburg(x,...
        this.Order,...
        opts{:});  % NFFT, Fs(?), and SpectrumType
      
      % If input is single, W will also be single so cast it to double
      W = double(W);
      
    end
    
  end  %% public methods
  
end  % classdef

