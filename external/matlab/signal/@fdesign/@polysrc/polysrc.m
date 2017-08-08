function this = polysrc(L,M,response,varargin)
%POLYSRC Construct a polynomial sample-rate converter (POLYSRC) filter designer.
%   D = FDESIGN.POLYSRC(L,M) constructs a polynomial sample-rate converter
%   filter designer D with an interpolation factor L and a decimation factor
%   M.  If L is not specified, it defaults to 3.  If M is not specified it
%   defaults to 2. Notice that L and M can be arbitrary positive numbers.
%
%   D = FDESIGN.POLYSRC(L,M,'Fractional Delay') initializes the filter
%   designer 'Response' property with 'Fractional Delay'.
%
%   D = FDESIGN.POLYSRC(L,M,'Fractional Delay','Np',Np) initializes the
%   filter designer specification with 'Np'and sets the polynomial order to
%   the value Np. If omitted Np defaults to 3.
%
%   D = FDESIGN.POLYSRC(...,Fs) specifies the sampling frequency (in Hz).
%
%   Example
%      %Design sample-rate converter that uses a 3rd order Lagrange 
%      %interpolation filter to convert from 44.1kHz to 48kHz.
%      [L,M] = rat(48/44.1);
%      f = fdesign.polysrc(L,M,'Fractional Delay','Np',3);
%      Hm = design(f,'lagrange');
%      Fs = 44.1e3;                         % Original sampling frequency
%      n = 0:9407;                          % 9408 samples, 0.213 seconds long
%      x  = sin(2*pi*1e3/Fs*n);             % Original signal, sinusoid at 1kHz
%      y = filter(Hm,x);                    % 10241 samples, still 0.213 seconds
%      stem(n(1:45)/Fs,x(1:45))             % Plot original sampled at 44.1kHz
%      hold on
%      % Plot fractionally interpolated signal (48kHz) in red
%      stem((n(3:51)-2)/(Fs*L/M),y(3:51),'r','filled') 
%      xlabel('Time (sec)');ylabel('Signal value')
%      legend('44.1 kHz sample rate','48 kHz sample rate')
%
%   For more information about Farrow SRCs, see the
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\efficientsrcdemo.html'])">Efficient Sample Rate Conversion between Arbitrary Factors</a> demo. 
%

%   Copyright 2007-2010 The MathWorks, Inc.

this = fdesign.polysrc;
set(this, 'MultirateType', 'Polynomial Sample Rate Converter');
set(this,'Response','Fractional Delay')

if nargin > 0
    set(this, 'InterpolationFactor', L);
    if nargin > 1
        set(this, 'DecimationFactor', M);
        if nargin > 2
            if strcmpi(lower(response),'Fractional Delay'),
                set(this,'Response','Fractional Delay')
            else
                error(message('signal:fdesign:polysrc:polysrc:InvalidResponse', response));
            end
        end
    end
end
setspecs(this, varargin{:});


% [EOF]
