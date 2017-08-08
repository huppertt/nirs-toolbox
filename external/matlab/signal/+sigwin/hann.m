classdef hann < sigwin.samplingflagwin
  %HANN Construct a Hann window object
  %
  %   SIGWIN.HANN is not recommended.  Use <a href="matlab:help hann">hann</a> instead.
  %
  %   H = SIGWIN.HANN(N, S) constructs a Hann window object with length N and
  %   sampling flag S.  If N or S is not specified, they default to 64 and
  %   'symmetric' respectively.  The sampling flag can also be 'periodic'.
  
  
  
  methods  % constructor block
    function hWIN = hann(n, sflag)
      
      
      % hWIN = sigwin.hann;
      hWIN.Name = 'Hann';
      
      if nargin > 0 && isnumeric(n),
        hWIN.Length = n;
      end
      if nargin > 1, hWIN.SamplingFlag = sflag; end
      
      
    end  % hann
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Hann window
      %
      %   sigwin.hann is not recommended.
      %   Use <a href="matlab:help hann">hann</a> instead.
      
      
      data = hann(hWIN.Length, hWIN.SamplingFlag);
      
    end
    
  end  %% public methods
  
end  % classdef

