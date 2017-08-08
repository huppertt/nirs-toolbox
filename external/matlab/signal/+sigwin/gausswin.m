classdef gausswin < sigwin.parameterizewin
  %GAUSSWIN Construct a Gaussian object
  %
  %   SIGWIN.GAUSSWIN is not recommended.  Use <a href="matlab:help gausswin">gausswin</a> instead.
  %
  %   H = SIGWIN.GAUSSWIN(N, A) constructs a Gaussian window object with
  %   length N and Alpha A.  If N or A is not specified, they default to 64
  %   and 2.5 respectively.
  
  
  
  methods  % constructor block
    function hWIN = gausswin(n, param)
      
      hWIN.Name = 'Gaussian';
      createdynamicprops(hWIN, 'Alpha', 'double','Alpha');
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      if nargin>1,
        hWIN.Alpha = param;
      else
        hWIN.Alpha = 2.5;
      end
      
      
    end  % gausswin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Gaussian window
      %
      %   sigwin.gausswin is not recommended.
      %   Use <a href="matlab:help gausswin">gausswin</a> instead.
      
      
      data = gausswin(hWIN.Length, hWIN.Alpha);
    end
    
  end  %% public methods
  
end  % classdef

