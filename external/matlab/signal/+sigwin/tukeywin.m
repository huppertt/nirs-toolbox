classdef tukeywin < sigwin.parameterizewin
  %TUKEYWIN Construct a Tukey object
  %
  %   SIGWIN.TUKEYWIN is not recommended.  Use <a href="matlab:help tukeywin">tukeywin</a> instead.
  %
  %   H = SIGWIN.TUKEYWIN(N, A) constructs a Tukey window object with length
  %   N and Alpha A.  If N or A is not specified, they default to 64 and .5
  %   respectively.
  
  
  
  methods  % constructor block
    function hWIN = tukeywin(n, param)
      
      hWIN.Name = 'Tukey';
      createdynamicprops(hWIN, 'Alpha', 'double','Alpha');
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      if nargin>1,
        hWIN.Alpha = param;
      else
        hWIN.Alpha = 0.5;
      end
      
      
    end  % tukeywin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Tukey window
      %
      %   sigwin.tukeywin is not recommended.
      %   Use <a href="matlab:help tukeywin">tukeywin</a> instead.
      
      
      data = tukeywin(hWIN.Length, hWIN.Alpha);
    end
    
  end  %% public methods
  
end  % classdef

