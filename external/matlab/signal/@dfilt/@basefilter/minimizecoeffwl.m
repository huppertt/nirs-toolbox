function Hbest = minimizecoeffwl(this,varargin)
%MINIMIZECOEFFWL Minimize coefficient wordlength.
%   Hq = MINIMIZECOEFFWL(Hd) returns the minimum wordlength fixed-point 
%   filter that meets the design specifications of the single or 
%   multistage FIR filter Hd. 
%
%   For multistage designs, the wordlength for each stage is minimized
%   individually.
%
%   Noise shaping is used by default in order to reduce the wordlength
%   required.
%
%   The noise shaping procedure is stochastic. To obtain repeatable results
%   on successive function calls, the uniform random number generator RAND
%   should be initialized prior to calling MINIMIZECOEFFWL.
%
%   Hd must be generated with FDESIGN and contain its design specifications.
%   Use getfdesign(Hd) to verify the design specifications.
%
%   Hq = MINIMIZECOEFFWL(Hd,...,'NTrials',N) specifies the number of Monte
%   Carlo trials. Hq is the filter with the smallest wordlength from all
%   trials that meets the design specifications for Hd. The default number
%   of trials is one.
%
%   Hq = MINIMIZECOEFFWL(Hd,...,'NoiseShaping',NSFlag) specifies whether to
%   enable or disable noise shaping. NSFlag can be either true or false.
%   If unspecified, 'NoiseShaping' defaults to true.
%
%   Hq = MINIMIZECOEFFWL(Hd,...,'MatchRefFilter',MRFFlag) specifies whether
%   the wordlength minimization only takes into account the specs used to
%   obtain Hd or if other aspects of Hd should be matched as well. The
%   filter order or the filter transition width will be matched if MRFFlag
%   is true. If unspecified, 'MatchRefFilter' defaults to false.
%
%   Hq = MINIMIZECOEFFWL(Hd,...,'Apasstol',Apasstol,'Astoptol',Astoptol)
%   specifies the passband and stopband tolerances. The default Apasstol is
%   1e-4, while the default Astoptol is 1e-2. Both are specified in dB.
%
%   Note: due to the stochastic aspect of noise shaping, the examples below
%   may not always produce the results indicated. You may want to
%   experiment with the 'NTrials' option and/or different seeds to
%   initialize RAND in order to reproduce the results. 
%
%   % Example #1: Minimize coefficient wordlength for a lowpass design.
%       Hf = fdesign.lowpass('Fp,Fst,Ap,Ast',.3,.42,.5,80);
%       Hd = design(Hf,'equiripple');
%       Hq = minimizecoeffwl(Hd); 
%
%   % Example #2: Minimize wordlength on each of 4 stages of a design.
%   % Compare with and without the use of noise shaping.
%       Hf  = fdesign.decimator(16,'lowpass','Fp,Fst,Ap,Ast',.05,.06,.75,65);
%       Hd  = design(Hf,'multistage');
%       Hq1  = minimizecoeffwl(Hd); % Wordlengths: 10,9,11,13
%       Hq2  = minimizecoeffwl(Hd,'NoiseShaping',false); % WL: 10,11,13,16
%
%   % Example #3: Save 1 bit and 4 coefficients by using noise shaping.
%       Hf  = fdesign.halfband('TW,Ast',.08,59);
%       Hd  = design(Hf,'kaiserwin');
%       Hq1 = minimizecoeffwl(Hd,'noiseShaping',true);  % 91 13-bit coeffs.
%       Hq2 = minimizecoeffwl(Hd,'noiseShaping',false); % 95 14-bit coeffs.
%
%   % Example #4: Save one bit by increasing tolerance.
%       Hf = fdesign.lowpass('Fp,Fst,Ap,Ast',.1,.12,.5,80);
%       Hd = design(Hf,'equiripple','UniformGrid',false);
%       Hq1 = minimizecoeffwl(Hd,'Astoptol',0);    % 21 bits
%       Hq2 = minimizecoeffwl(Hd,'Astoptol',7e-2); % 20 bits
%
%   See also CONSTRAINCOEFFWL, MAXIMIZESTOPBAND, RAND.

%   Copyright 2009 The MathWorks, Inc.

% Make sure to work with reference filter in case filter has been quantized
Href = reffilter(this);

try
    Hbest = minimizecoeffwlfir(this,Href,varargin{:});
catch ME
    throw(ME);
end




% [EOF]
