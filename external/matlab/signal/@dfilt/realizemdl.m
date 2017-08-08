%REALIZEMDL Filter realization (Simulink diagram).
%     REALIZEMDL(Hd) automatically generates architecture model of filter
%     Hd in a Simulink subsystem block using individual sum, gain, and
%     delay blocks, according to user-defined specifications.
%
%     REALIZEMDL(Hd, PARAMETER1, VALUE1, PARAMETER2, VALUE2, ...) generates
%     the model with parameter/value pairs. Valid properties and values for
%     realizemdl are listed in this table, with the default value noted and
%     descriptions of what the properties do.
%
%     -------------       ---------------       ----------------------------
%     Property Name       Property Values       Description
%     -------------       ---------------       ----------------------------
%     Destination         [{'current'}          Specify whether to add the block
%                          'new'                to your current Simulink model,
%                          <user defined>]      create a new model to contain the
%                                               block, or specify the name of the
%                                               target subsystem. 
%
%     Blockname           {'filter'}            Provides the name for the new 
%                                               subsystem block. By default the 
%                                               block is named 'filter'.
%
%     OverwriteBlock      ['on' {'off'}]        Specify whether to overwrite an
%                                               existing block with the same name
%                                               as specified by the Blockname 
%                                               property or create a new block.
%
%     OptimizeZeros       [{'on'} 'off']        Specify whether to remove zero-gain
%                                               blocks.
% 
%     OptimizeOnes        [{'on'} 'off']        Specify whether to replace unity-gain
%                                               blocks with direct connections.
%
%     OptimizeNegOnes     [{'on'} 'off']        Specify whether to replace negative
%                                               unity-gain blocks with a sign
%                                               change at the nearest sum block.
%
%     OptimizeDelayChains [{'on'} 'off']        Specify whether to replace cascaded 
%                                               chains of delay blocks with a
%                                               single integer delay block to 
%                                               provide an equivalent delay.
%
%     MapStates           ['on' {'off'}]        Specify whether to map the States
%                                               of the filter as initial conditions
%                                               of the block.
%
%     MapCoeffsToPorts    ['on' {'off'}]        Specify whether to map the 
%                                               coefficients of the filter
%                                               to the ports of the block.
%                                            
%     CoeffNames          {'Num'}               Specify the coefficient 
%                                               variables names. MapCoeffsToPorts 
%                                               must be 'on' for this property 
%                                               to apply.
%
%     InputProcessing(*) [{'columnsaschannels'} Specify sample-based 
%                        'elementsaschannels'   ('elementsaschannels') vs. 
%                        'inherited']           frame-based ('columnsaschannels') 
%                                               processing.
%
%    (*) A DSP System Toolbox is required when you specify InputProcessing
%    as 'columnsaschannels'. Input processing is set to frame-based
%    ('columnsaschannels') by default. If you do not have a license,
%    explicitly set the input processing to sample-based by setting the
%    'InputProcessing' parameter to 'elementsaschannels'.
%
%    % Examples:
%    [b,a] = butter(4,.4);
%    Hd = dfilt.df1(b,a);
% 
%    %#1 Default syntax:
%    realizemdl(Hd);
% 
%    %#2 Using parameter/value pairs:
%    realizemdl(Hd, 'OverwriteBlock', 'on');
%
%    %#3 Using map coefficients to ports
%    realizemdl(Hd, 'MapCoeffsToPorts', 'on','CoeffNames',{'Num','Den'});
 
%    Copyright 1988-2012 The MathWorks, Inc.


