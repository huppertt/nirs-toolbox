classdef triang < sigwin.simplewin 
  %sigwin.triang class
  %   sigwin.triang extends sigwin.simplewin.
  %
  %    sigwin.triang properties:
  %       Name - Property is of type 'string' (read only)
  %       Length - Property is of type 'spt_uint32 user-defined'
  %
  %    sigwin.triang methods:
  %       generate - hWIN) Generates the triangular window
  
  
  
  methods  % constructor block
    function hWIN = triang(n)
      %TRIANG Construct a Triangular window object
      %
      %   SIGWIN.TRIANG is not recommended.  Use <a href="matlab:help triang">triang</a> instead.
      %
      %   H = SIGWIN.TRIANG(N) constructs a Triangular window object with length
      %   N.  If N is not specified, it defaults to 64.
      
      %   Author(s): V.Pellissier
      
      % hWIN = sigwin.triang;
      hWIN.Name = 'Triangular';
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % triang
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the triangular window
      %
      %   sigwin.triang is not recommended.
      %   Use <a href="matlab:help triang">triang</a> instead.
      
      %   Author(s): V.Pellissier
      %   Copyright 1988-2012 The MathWorks, Inc.
      
      data = triang(hWIN.Length);
      
    end
    
  end  %% public methods
  
end  % classdef

