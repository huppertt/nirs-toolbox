%COEWRITE Write a XILINX CORE Generator(tm) coefficient (.COE) file.
%   COEWRITE(Hd) writes a XILINX Distributed Arithmetic FIR filter
%   coefficient .COE file which can be loaded into the XILINX CORE
%   Generator.  The coefficients  are extracted from the fixed-point DFILT
%   filter object, Hd. A dialog box is displayed  to fill in a file name.
%   The default file name is 'untitled.coe'. 
%
%   COEWRITE(Hd,RADIX) indicates the radix (number base) being used to
%   specify  the FIR filter coefficients.  Valid RADIX values are 2, 10,
%   and 16 (default).
%
%   COEWRITE(...,FILENAME) writes a XILINX .COE file to a disk file  called
%   FILENAME.
%
%   The extension '.coe' will be added to FILENAME if it doesn't already
%   have  an extension.
%
%   EXAMPLE:
%   b = firceqrip(30,0.4,[0.05 0.03]);
%   Hd = dfilt.dffir(b); 
%   Hd.arithmetic = 'fixed'; % Requires the Fixed-Point Designer
%   coewrite(Hd,10,'mycoefile');
%
%   See also COEREAD.

%   Author: P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

% Help for the DFILT/COEWRITE method.

% [EOF]
