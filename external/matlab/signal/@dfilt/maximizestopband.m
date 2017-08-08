%MAXIMIZESTOPBAND Maximize stopband attenuation.
%   Hq = MAXIMIZESTOPBAND(Hd,WL) quantizes the filter Hd and returns the
%   fixed-point filter Hq using wordlength WL to represent the filter 
%   coefficients. The stopband attenuation of Hq is maximized for the given
%   wordlength WL.
%
%   For multistage filters, WL can either be a scalar or vector. If WL is a
%   scalar, the same wordlength is used for all stages. If WL is a vector,
%   each stage uses the corresponding element in the vector. The vector
%   length must equal the number of stages.
%
%   Noise shaping is used to maximize the stopband attenuation of the
%   filter Hq. As a trade-off, the passband ripple usually increases
%   slightly.
%
%   The noise shaping procedure is stochastic. To obtain repeatable results
%   on successive function calls, the uniform random number generator RAND
%   should be initialized prior to calling MAXIMIZESTOPBAND.
%
%   Hd must be generated with FDESIGN and contain its design specifications.
%   Use getfdesign(Hd) to verify the design specifications.
%
%   Hq = MAXIMIZESTOPBAND(Hd,WL,'NTrials',N) specifies the number of Monte
%   Carlo trials. Hq is the filter with the largest stopband attenuation
%   among the trials. The default number of trials is one.
%
%   Note: due to the stochastic aspect of noise shaping, the examples below
%   may not always produce the results indicated. You may want to
%   experiment with the 'NTrials' option and/or different seeds to
%   initialize RAND in order to reproduce the results. 
%
%   % Example #1: Compare the stopband attenuation of a quantized filter 
%   % with and without stopband maximization. Results for Hq may vary.
%       Hf = fdesign.lowpass('Fp,Fst,Ap,Ast',0.4,0.45,0.5,60);
%       Hd = design(Hf,'equiripple');
%       Hd.Arithmetic = 'fixed';
%       WL = 16; % use 16 bits to represent coefficients
%       Hd.CoeffWordLength = WL;
%       Hq = maximizestopband(Hd,WL);
%       md = measure(Hd); % Apass = 0.494 dB; Astop = 59.1436 dB
%       mq = measure(Hq); % Apass = 0.5 dB;   Astop = 60 dB 
%
%   % Example #2: Obtain close to 70 dB attenuation with a 70th order
%   % halfband decimator while using 14 bits to represent the coefficients.
%       Hf = fdesign.decimator(2,'halfband','N,Ast',70,80);
%       Hd = design(Hf,'equiripple');
%       Hq = maximizestopband(Hd,14);
%       mq = measure(Hq); % Astop = 69.4 dB (compare to 66.2 dB)
%
%   % Example #3: Maximize stopband for a 4-stage decimator. Use a
%   % different wordlength for each stage.
%       Hf = fdesign.decimator(16,'lowpass',0.05,0.06,0.75,65);
%       Hd = design(Hf,'multistage');
%       WL = [12 11 13 17];
%       Hq  = maximizestopband(Hd,WL);
%
%   See also CONSTRAINCOEFFWL, MINIMIZECOEFFWL, RAND.

%   Copyright 2009 The MathWorks, Inc.

% [EOF]
