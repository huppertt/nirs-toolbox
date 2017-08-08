classdef blackman < sigwin.samplingflagwin
  %BLACKMAN Construct a Blackman window object
  %
  %   SIGWIN.BLACKMAN is not recommended.  Use <a href="matlab:help blackman">blackman</a> instead.
  %
  %   H = SIGWIN.BLACKMAN(N, S) constructs a Blackman window object with
  %   length N and sampling flag S.  If N or S is not specified, they default
  %   to 64 and 'symmetric' respectively.  The sampling flag can also be
  %   'periodic'.
  
    
  methods  % constructor block
    function hWIN = blackman(n, sflag)
      
      
      % hWIN = sigwin.blackman;
      hWIN.Name = 'Blackman';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      if nargin>1,
        hWIN.SamplingFlag = sflag;
      end
      
      
    end  % blackman
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Blackman window
      %
      %   sigwin.blackman is not recommended.
      %   Use <a href="matlab:help blackman">blackman</a> instead.
      
      data = blackman(hWIN.Length, hWIN.SamplingFlag);
    end
    
  end  %% public methods
  
end  % classdef

