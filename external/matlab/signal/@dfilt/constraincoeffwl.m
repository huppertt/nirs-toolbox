%CONSTRAINCOEFFWL Constrain coefficient wordlength.
%   Hq = CONSTRAINCOEFFWL(Hd,WL) returns the fixed-point filter Hq that
%   meets the design specifications of the single or multistage FIR filter
%   Hd using at most a wordlength WL to represent the filter coefficients.
%
%   For multistage designs, WL can either be a scalar or vector. If WL is a
%   scalar, the same wordlength constraint is applied to all stages. If WL
%   is a vector, each stage is constrained by the corresponding element in
%   the vector. The vector length must equal the number of stages.
%
%   Hd must be generated with FDESIGN and contain its design specifications.
%   Use getfdesign(Hd) to verify the design specifications.
%
%   The wordlength constraint may result in a fixed-point filter with a 
%   different order than the input filter object. 
%
%   Noise shaping is used by default in order to satisfy the wordlength
%   constraint. As a trade-off, the passband ripple usually increases
%   slightly.
%
%   The noise shaping procedure is stochastic. To obtain repeatable results
%   on successive function calls, the uniform random number generator RAND
%   should be initialized prior to calling CONSTRAINCOEFFWL.
%
%   Hq = CONSTRAINCOEFFWL(Hd,WL,...,'NTrials',N) specifies the number of
%   Monte Carlo trials. Hq is the first filter among the trials that meets
%   the specifications in Hd while using at most a wordlength WL. The
%   default number of trials is one.
%
%   Hq = CONSTRAINCOEFFWL(Hd,WL,...,'NoiseShaping',NSFlag) specifies
%   whether to enable or disable noise shaping. NSFlag can be either true
%   or false. If unspecified, 'noiseShaping' defaults to true.
%
%   Hq = CONSTRAINCOEFFWL(...,'Apasstol',Apasstol,'Astoptol',Astoptol)
%   specifies the passband and stopband tolerances. The default Apasstol is
%   1e-4, while the default Astoptol is 1e-2. Both are specified in dB.
%
%   Note: due to the stochastic aspect of noise shaping, the examples below
%   may not always produce the results indicated. You may want to
%   experiment with the 'NTrials' option and/or different seeds to
%   initialize RAND in order to reproduce the results. 
%
%   % Example #1: Constrain filter to meet specs using at most 11-bit
%   % coefficients
%       Hf = fdesign.lowpass('Fp,Fst,Ap,Ast',.4,.5,1,60);
%       Hd = design(Hf,'equiripple'); % 43 coefficients
%       Hq = constraincoeffwl(Hd,11); % 45 11-bit coefficients
%
%   % Example #2: Constrain halfband filter to meet specs using at most 
%   % 16-bit coefficients
%       Hf = fdesign.halfband('TW,Ast',0.2,80);
%       Hd = design(Hf,'kaiserwin');  % 53 coefficients
%       Hq = constraincoeffwl(Hd,16); % 57 16-bit coefficients
%
%   % Example #3: Constrain each stage of a multistage filter to use 9, 13,
%   % and 14 bits respectively. Overall, the multistage filter meets specs.
%       Hf = fdesign.decimator(8,'lowpass','Fp,Fst,Ap,Ast',0.1,0.12,1,70);
%       Hd  = design(Hf,'multistage','nstages',3);
%       WL = [9 13 14];
%       Hq = constraincoeffwl(Hd,WL); % Try: measure(Hq)
%
%   See also MINIMIZECOEFFWL, MAXIMIZESTOPBAND, RAND.

%   Copyright 2009 The MathWorks, Inc.

% [EOF]
