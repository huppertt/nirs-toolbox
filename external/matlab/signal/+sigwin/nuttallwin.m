classdef nuttallwin < sigwin.samplingflagwin
  %NUTTALLWIN Construct a Nuttall defined minimum 4-term Blackman-Harris window object
  %
  %   SIGWIN.NUTTALLWIN is not recommended.  Use <a href="matlab:help nuttallwin">nuttallwin</a> instead.
  %
  %   H = SIGWIN.NUTTALLWIN(N) constructs a Nuttall defined minimum 4-term
  %   Blackman-Harris window object with length N.  If N is not specified, it
  %   defaults to 64.
  
  
  
  methods  % constructor block
    function hWIN = nuttallwin(n)
      
      hWIN.Name = 'Nuttall';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % nuttallwin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Nuttall window
      %
      %   sigwin.nuttallwin is not recommended.
      %   Use <a href="matlab:help nuttallwin">nuttallwin</a> instead.
      
      
      data = nuttallwin(hWIN.Length, hWIN.SamplingFlag);
    end
    
  end  %% public methods
  
end  % classdef

