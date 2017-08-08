classdef taylorwin < sigwin.parameterizewin
  %TAYLORWIN Construct a TAYLORWIN object
  %
  %   SIGWIN.TAYLORWIN is not recommended.  Use <a href="matlab:help taylorwin">taylorwin</a> instead.
  %
  %   H = SIGWIN.TAYLORWIN(N, NBAR, SLL) constructs a Taylor window object
  %   with length N, number of nearly constant-level sidelobes adjacent to
  %   the mainlobe NBAR, and maximum sidelobe level SLL (in dB).  If N, NBAR,
  %   and SLL are not specified, their default values are 64, 4, and 30
  %   respectively.
  
  
  
  methods  % constructor block
    function hWin = taylorwin(n, nbar, sll)
      
      
      
      narginchk(0,3);
      
      % hWin = sigwin.taylorwin;
      hWin.Name = 'Taylor';
      
      createdynamicprops(hWin, 'Nbar', ...
        'posint', 'Number of nearly constant-level sidelobes');
      
      createdynamicprops(hWin, 'SidelobeLevel', ...
        'udouble', 'Maximum sidelobe level');
      
      if nargin>0 && isnumeric(n),
        hWin.Length = n;
      end
      
      if nargin>1,
        hWin.Nbar = nbar;
      else
        hWin.Nbar = 4;
      end
      
      if nargin>2,
        hWin.SidelobeLevel = sll;
      else
        hWin.SidelobeLevel = 30;
      end
      
      
    end  % taylorwin
    
    function data = generate(hWin)
      %Generate(hWin) Generates the Taylor Window.
      %
      %   sigwin.taylorwin is not recommended.
      %   Use <a href="matlab:help taylorwin">taylorwin</a> instead.
      
      data = taylorwin(hWin.Length, hWin.Nbar, -hWin.SidelobeLevel);
    end
    
  end  %% public methods
  
end  % classdef

