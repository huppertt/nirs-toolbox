%GENERATEHDL Generate HDL.
%   GENERATEHDL(Hb) automatically generates VHDL or Verilog code for 
%   the quantized filter, Hb. The default language is VHDL; to generate
%   Verilog, supply the property/value pair 'TargetLanguage','Verilog'. 
%   The default file name is the name of the filter variable, e.g. 
%   Hb.vhd for VHDL and Hb.v for Verilog.  The file is written to 
%   the HDL source directory which defaults to 'hdlsrc' under the 
%   current directory. This directory will be created if necessary.
%
%   GENERATEHDL(Hb, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...)
%   generates the HDL code with parameter/value pairs. Valid properties
%   and values for GENERATEHDL are listed in the Filter Design HDL 
%   Coder documentation section "Property Reference."
%
%   EXAMPLE:
%   filtdes = fdesign.lowpass('N,Fc,Ap,Ast',30,0.4,0.05,0.03,'linear');
%   Hb = design(filtdes,'filterstructure','dffir');
%   Hb.arithmetic = 'fixed';
%   generatehdl(Hb);
%   generatehdl(Hb,'TargetLanguage','Verilog');
%   generatehdl(Hb,'Name','myfiltername','TargetDirectory','mysrcdir');
%   generatehdl(Hb,'InputPort','adc_data','OutputPort','dac_data');
%   generatehdl(Hb,'AddInputRegister','on','AddOutputRegister','off');
%   generatehdl(Hb,'OptimizeForHDL','on','CoeffMultipliers','csd');
%   generatehdl(Hb,'SerialPartition', 31);
%   generatehdl(Hb,'SerialPartition', [12 11 8], 'ReuseAccum', 'on');
%
%   See also GENERATETB, GENERATETBSTIMULUS, FDHDLTOOL.

%   Copyright 2003-2011 The MathWorks, Inc.

% [EOF]

