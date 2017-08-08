classdef hamming < sigwin.samplingflagwin
  %HAMMING Construct a Hamming window object
  %
  %   SIGWIN.HAMMING is not recommended.  Use <a href="matlab:help hamming">hamming</a> instead.
  %
  %   H = SIGWIN.HAMMING(N, S) constructs a Hamming window object with length
  %   N and sampling flag S.  If N or S is not specified, they default to 64
  %   and 'symmetric' respectively.  The sampling flag can also be
  %   'periodic'.
  
  
  
  methods  % constructor block
    function hWIN = hamming(n, sflag)
      
      hWIN.Name = 'Hamming';
      
      if nargin > 0 && isnumeric(n),
        hWIN.Length       = n;
      end
      
      if nargin > 1, hWIN.SamplingFlag = sflag; end
      
      
    end  % hamming
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Hamming window
      %
      %   sigwin.hamming is not recommended.
      %   Use <a href="matlab:help hamming">hamming</a> instead.
      
      
      data = hamming(hWIN.Length, hWIN.SamplingFlag);
    end
    
  end  %% public methods
  
end  % classdef

