classdef bartlett < sigwin.simplewin
  %BARTLETT Bartlett window.
  %
  %   SIGWIN.BARTLETT is not recommended.  Use <a href="matlab:help bartlett">bartlett</a> instead.
  %
  %   H = SIGWIN.BARTLETT(N) returns a N-point Bartlett window object H.
  %
  %   EXAMPLE:
  %     N = 64;
  %     h = sigwin.bartlett(N);
  %     w = generate(h);
  %     stem(w); title('64-point Bartlett Window');
  
  
  
  
  methods  % constructor block
    function hWIN = bartlett(n)
      
      
      % hWIN = sigwin.bartlett;
      hWIN.Name = 'Bartlett';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % bartlett
    
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Bartlett window
      %
      %   sigwin.bartlett is not recommended.
      %   Use <a href="matlab:help bartlett">bartlett</a> instead.
      
      
      data = bartlett(hWIN.Length);
    end
  end  %% public methods
  
end  % classdef

