classdef rectwin < sigwin.simplewin
  %RECTWIN Construct a Rectangular window object
  %
  %   SIGWIN.RECTWIN is not recommended.  Use <a href="matlab:help rectwin">rectwin</a> instead.
  %
  %   H = SIGWIN.RECTWIN(N) constructs a Rectangular window object with length
  %   N.  If N is not specified, it defaults to 64.
  
  
  
  methods  % constructor block
    function hWIN = rectwin(n)
      
      hWIN.Name = 'Rectangular';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % rectwin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the rectangular window
      %
      %   sigwin.rectwin is not recommended.
      %   Use <a href="matlab:help rectwin">rectwin</a> instead.
      
      
      data = rectwin(hWIN.Length);
      
    end
    
  end  %% public methods
  
end  % classdef

