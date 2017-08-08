classdef parzenwin < sigwin.simplewin
  %PARZENWIN Construct a Parzen window object
  %
  %   SIGWIN.PARZENWIN is not recommended.  Use <a href="matlab:help parzenwin">parzenwin</a> instead.
  %
  %   H = SIGWIN.PARZENWIN(N) constructs a Parzen window object with length
  %   N.  If N is not specified, it defaults to 64.
  
  
  methods  % constructor block
    function hWIN = parzenwin(n)
      
      hWIN.Name = 'Parzen';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % parzenwin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Parzen window
      %
      %   sigwin.parzenwin is not recommended.
      %   Use <a href="matlab:help parzenwin">parzenwin</a> instead.
      
      data = parzenwin(hWIN.Length);
    end
    
  end  %% public methods
  
end  % classdef

