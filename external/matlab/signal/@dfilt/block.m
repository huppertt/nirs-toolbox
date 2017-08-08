%BLOCK Generate a DSP System Toolbox block equivalent to the filter object.
%   BLOCK(Hd) generates a DSP System Toolbox block equivalent to Hd.
%
%   BLOCK(Hd, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) generates a
%   DSP System Toolbox block using the options specified in the
%   parameter/value pairs. The available parameters are:
%
%     -------------       ---------------      ----------------------------
%     Property Name       Property Values      Description
%     -------------       ---------------      ----------------------------
%     Destination         [{'current'}         Specify whether to add the block
%                          'new'               to your current Simulink model,
%                          <user defined>]     create a new model to contain the
%                                              block, or specify the name of the
%                                              target subsystem. 
%
%     Blockname           {'filter'}           Provides the name for the new 
%                                              subsystem block. By default the 
%                                              block is named 'filter'.
%
%     OverwriteBlock      ['on' {'off'}]       Specify whether to overwrite an
%                                              existing block with the same name
%                                              as specified by the Blockname 
%                                              property or create a new block.
%
%     MapStates           ['on' {'off'}]       Specify whether to map the States
%                                              of the filter as initial conditions
%                                              of the block.
%
%     Link2Obj            ['on' {'off'}]       Specify whether to set the
%                                              filter variable in the block
%                                              mask rather than setting the
%                                              coefficient values.
%
%     MapCoeffsToPorts    ['on' {'off'}]       Specify whether to map the 
%                                              coefficients of the filter 
%                                              to the ports of the block.
%
%     CoeffNames          {'Num'}              Specify the coefficients
%                                              variables names. By default 
%                                              the coefficient variables take
%                                              the names of the ports.
%                                              MapCoeffsToPorts must be 'on' 
%                                              for this property to apply.
%
%     InputProcessing   [{'columnsaschannels'} Specify sample-based 
%                        'elementsaschannels'  ('elementsaschannels') vs. 
%                        'inherited']          frame-based ('columnsaschannels') 
%                                              processing.
%
%    EXAMPLES:
%    [b,a] = butter(5,.5);
%    Hd = dfilt.df1(b,a);
% 
%    %#1 Default syntax:
%    block(Hd);
% 
%    %#2 Using parameter/value pairs:
%    block(Hd, 'Blockname', 'DF1');

% Copyright 1988-2010 The MathWorks, Inc.

% Help for the filter's BLOCK method.

% [EOF]
