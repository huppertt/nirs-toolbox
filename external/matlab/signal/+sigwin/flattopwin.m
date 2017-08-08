classdef flattopwin < sigwin.samplingflagwin
  %FLATTOPWIN Construct a Flat Top window object
  %
  %   SIGWIN.FLATTOPWIN is not recommended.  Use <a href="matlab:help flattopwin">flattopwin</a> instead.
  %
  %   H = SIGWIN.FLATTOPWIN(N, S) constructs a Flat Top window object with length N
  %   and sampling flag S.  If N or S is not specified, they default to 64 and
  %   'symmetric' respectively.  The sampling flag can also be 'periodic'.
  
  
  
  methods  % constructor block
    function this = flattopwin(n, sflag)
      
      
      % this = sigwin.flattopwin;
      this.Name = 'Flat Top';
      
      if nargin > 0 && isnumeric(n),
        this.Length = n;
      end
      
      if nargin > 1,
        this.SamplingFlag = sflag;
      end
      
      
    end  % flattopwin
    
    
    function data = generate(hWIN)
      %GENERATE Generates the Flat Top window
      %
      %   sigwin.flattopwin is not recommended.
      %   Use <a href="matlab:help flattopwin">flattopwin</a> instead.
      
      
      data = flattopwin(hWIN.Length, hWIN.SamplingFlag);
      
    end
  end  %% public methods
  
end  % classdef

