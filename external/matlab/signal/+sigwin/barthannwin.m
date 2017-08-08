classdef barthannwin < sigwin.simplewin
  %BARTHANNWIN Modified Bartlett-Hanning window.
  %
  %   SIGWIN.BARTHANNWIN is not recommended.  Use <a href="matlab:help barthannwin">barthannwin</a> instead.
  %
  %   H = SIGWIN.BARTHANNWIN(N) returns an N-point Modified Bartlett-Hanning
  %   window object H.
  %
  %   EXAMPLE:
  %      N = 64;
  %      h = sigwin.barthannwin(N);
  %      w = generate(h);
  %      stem(w); title('64-point Modified Bartlett-Hanning window');
  
  %   Reference:
  %     [1] Yeong Ho Ha and John A. Pearce, A New Window and Comparison
  %         to Standard Windows, IEEE Transactions on Acoustics, Speech,
  %         and Signal Processing, Vol. 37, No. 2, February 1999
  
  
  
  methods  % constructor block
    function hWIN = barthannwin(n)
      
      
      % hWIN = sigwin.barthannwin;
      hWIN.Name = 'Bartlett-Hanning';
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      
    end  % barthannwin
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Bartlett-Hanning window
      %
      %   sigwin.barthannwin is not recommended.
      %   Use <a href="matlab:help barthannwin">barthannwin</a> instead.
      
      
      data = barthannwin(hWIN.Length);
    end
    
  end  %% public methods
  
end  % classdef

