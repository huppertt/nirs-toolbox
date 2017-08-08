classdef blackmanharris < sigwin.samplingflagwin
  %BLACKMANHARRIS Construct a Blackman-Harris window object
  %
  %   SIGWIN.BLACKMANHARRIS is not recommended.  Use <a href="matlab:help blackmanharris">blackmanharris</a> instead.
  %
  %   H = SIGWIN.BLACKMANHARRIS(N) constructs a Blackman-Harris window object
  %   with length N.  If N is not specified, it defaults to 64.
  
  
  
  methods  % constructor block
    function hWIN = blackmanharris(n)
      
      
      % hWIN = sigwin.blackmanharris;
      hWIN.Name = 'Blackman-Harris';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % blackmanharris
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Blackman-Harris window
      %
      %   sigwin.blackmanharris is not recommended.
      %   Use <a href="matlab:help blackmanharris">blackmanharris</a> instead.
      
      
      data = blackmanharris(hWIN.Length, hWIN.SamplingFlag);
    end
    
  end  %% public methods
  
end  % classdef

