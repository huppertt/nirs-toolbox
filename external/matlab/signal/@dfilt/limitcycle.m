%LIMITCYLE  Detect zero-input limit cycles in IIR quantized filters.
%   LIMITCYCLE(Hd) runs 20 Monte Carlo trials with random initial states
%   and zero input vector of twice the length of the impulse response. This
%   function stops if a zero-input limit cycle is detected in the quantized
%   filter Hd. A limit cycle type is returned which is one of 'granular' to
%   indicate that a granular overflow occurred; 'overflow' to indicate that
%   an overflow limit cycle occurred; or 'none' to indicate that no limit
%   cycles were detected during the Monte Carlo trials.
%
%   LIMITCYCLE(Hd, NTRIALS, INPUTLENGTHFACTOR, STOPCRITERION) allows
%   you to set
%
%     NTRIALS, the number of Monte Carlo trials (default is 20).
%
%     INPUTLENGTHFACTOR, the length of the zero vector as a multiple of the
%     length of the impulse response of the filter (default is 2).
%
%     STOPCRITERION, the stop criterion, a string containing one of
%     'either' (the default), 'granular' or 'overflow'.  If STOPCRITERION
%     is 'either', then the Monte Carlo trials will stop if either a
%     granular or overflow limit cycle is detected; 'granular', stop only
%     if a granular limit cycle was detected; 'overflow', stop only if an
%     overflow limit cycle was detected.
%
%   If any of the input values are empty ([]), then the default
%   values are used.
%
%   REPORT = LIMITCYCLE(Hd, ...) also returns a REPORT object with fields:
%
%     LIMITCYCLETYPE, one of 'granular' to indicate that a granular
%     limit cycle occurred; 'overflow' to indicate that an overflow limit
%     cycle occurred; or 'none' to indicate that no limit cycles were
%     detected during the Monte Carlo trials.
%
%     ZI the initial conditions that caused the limit cycle.
%
%     OUTPUT the output of the filter in the limit cycle (steady state).
%
%     TRIAL the number of the Monte Carlo trial that was stopped on.
%
%   If no limit cycles are detected, then the parameters of the last Monte
%   Carlo trial are returned.
%
%   Example:
%     s = [1 0 0 1 0.9606 0.9849];
%     Hd = dfilt.df2sos(s);
%     Hd.Arithmetic = 'fixed';
%     greport = limitcycle(Hd,20,2,'granular')
%     oreport = limitcycle(Hd,20,2,'overflow')
%     figure,
%     subplot(211),plot(greport.Output(1:20)), title('Granular Limit Cycle');
%     subplot(212),plot(oreport.Output(1:20)), title('Overflow Limit Cycle');

%   Copyright 2005 The MathWorks, Inc.

% Help for the p-coded LIMITCYCLE method of DFILT classes.
