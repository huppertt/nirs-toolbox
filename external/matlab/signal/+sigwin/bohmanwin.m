classdef bohmanwin < sigwin.simplewin
  %BOHMANWIN Construct a Bohman window object
  %
  %   SIGWIN.BOHMANWIN is not recommended.  Use <a href="matlab:help bohmanwin">bohmanwin</a> instead.
  %
  %   H = SIGWIN.BOHMANWIN(N) constructs a Bohman window object with length
  %   N. If N is not specified, it defaults to 64.
  
  
  
  methods  % constructor block
    function hWIN = bohmanwin(n)
      
      % hWIN = sigwin.bohmanwin;
      hWIN.Name = 'Bohman';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % bohmanwin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Bohman window
      %
      %   sigwin.bohmanwin is not recommended.
      %   Use <a href="matlab:help bohmanwin">bohmanwin</a> instead.
      
      
      data = bohmanwin(hWIN.Length);
    end
    
  end  %% public methods
  
end  % classdef

