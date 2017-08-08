function inputdata = generatetbstimulus(filterobj, varargin)
%GENERATETBSTIMULUS Generates and returns HDL Test Bench Stimulus
%   GENERATETBSTIMULUS(Hd) automatically generates the filter
%   input stimulus based on the settings for the current filter.
%   The stimulus consists of any or all of impulse, step, ramp,
%   chirp, noise, or user-defined stimulus.  Note that the results
%   are quantized using the input quantizer settings of Hd. If the
%   arithmetic property of the filter Hd is set to 'double', then
%   double-precision stimulus is return. If the arithmetic property
%   is set to 'fixed', then the stimulus is return as a fixed-point
%   object.
%
%   GENERATETBSTIMULUS(Hd, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) 
%   generates the test bench with parameter/value pairs. Valid 
%   properties and values for GENERATETBSTIMULUS are listed in 
%   the HDL Filter Designer documentation section "Property Reference."
%
%   Y = GENERATETBSTIMULUS(Hd,...) returns the stimulus to MATLAB
%   variable Y. 
%
%   GENERATETBSTIMULUS(Hd,...) with no output argument plots the 
%   stimulus in the current figure window.
%
%   EXAMPLE:
%   h = firceqrip(30,0.4,[0.05 0.03]);
%   Hb = dfilt.dffir(h);
%   y = generatetbstimulus(Hb, 'TestBenchStimulus',{'ramp','chirp'});
%   generatetbstimulus(Hb, 'TestBenchStimulus','noise');
%
%   See also GENERATEHDL, GENERATETB, ISHDLABLE.

%   Copyright 2003-2012 The MathWorks, Inc.

  % check for Filter Design HDL Coder
  fdhdlcInstallCheck;
  
  [cando, ~, errObj] = ishdlable(filterobj);
  if ~cando
    error(errObj);
  end

  % Parse TestBenchStimulus and TestBenchUserStimulus parameters
  hprop = PersistentHDLPropSet; % Get persistent HDLPropSet
  if isempty(hprop),
      hprop = hdlcoderprops.HDLProps;
      PersistentHDLPropSet(hprop); % Make hprop persistent
      hdlsetparameter('tbrefsignals', false); % Default specific to FDHC
  end
  
  set(hprop.CLI,'TestbenchUserStimulus', []); % resets to "forget"
  set(hprop.CLI,'TestbenchStimulus', '');
  set(hprop.CLI, 'TestbenchFracDelayStimulus', '');
  set(hprop.CLI, 'TestbenchCoeffStimulus', []);
  
  set(hprop.CLI,varargin{:});
  updateINI(hprop);
  
  hF = createhdlfilter(filterobj);
  inputdata =  hF.maketbstimulus(filterobj);
  

% [EOF]

