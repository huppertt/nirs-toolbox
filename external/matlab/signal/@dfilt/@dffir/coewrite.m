function coewrite(Hq,radix,filename)
%COEWRITE Write a XILINX CORE Generator(tm) coefficient (.COE) file.
%   COEWRITE(Hd) writes a XILINX Distributed Arithmetic FIR filter
%   coefficient .COE file which can be loaded into the XILINX CORE
%   Generator.  The coefficients are extracted from the fixed-point DFILT
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

%   Author(s): P. Costa
%   Copyright 1999-2010 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct'));
if nargin < 2, radix = []; end
if isempty(radix), radix = 16; end
if nargin < 3, filename = []; end

if ~strcmpi(Hq.Arithmetic, 'Fixed'),
    error(message('signal:dfilt:dffir:coewrite:NotFixedPoint'));
end

if ~isreal(Hq),
    error(message('signal:dfilt:dffir:coewrite:NotReal'));
end

file = 'untitled.coe';
ext = 'coe';
dlgStr = fdatoolmessage('ExpQuantizedCoeffsToCOEFile');

if isempty(filename),
    % Put up the file selection dialog
    [filename, pathname] = lcluiputfile(file,dlgStr);
else
    % File will be created in present directory
    s = pwd;
    pathname = [s filesep];
end

if ~isempty(filename),
    if isempty(strfind(filename,'.')), filename=[filename '.' ext]; end
    deffile = [pathname filename];
 
    save2coefile(Hq, radix, deffile);
end


%------------------------------------------------------------------------
function save2coefile(Hq, radix, file)

% Unix returns a path that sometimes includes two paths (the
% current path followed by the path to the file) separated by '//'.
% Remove the first path.
indx = strfind(file,[filesep,filesep]);
if ~isempty(indx)
    file = file(indx+1:end);
end

% Write the coefficients out to a .COE file.
fid = fopen(file,'w');

% Display header information
strheader = sptfileheader(fdatoolmessage('XilinxCoreGenDistributedArithFIRCoeffFile'), ...
    'dsp', ';');
fprintf(fid,'%s\n',strheader);

% Display the Radix
fprintf(fid, 'Radix = %d; \n',radix); % 'Radix' is a .COE file keyword

s = get(Hq);

% Display the coefficient wordlength
fprintf(fid,'Coefficient_Width = %d; \n',s.CoeffWordLength);

% Get the coefficients from the object
qcoeffs = Hq.Numerator;

% Create a quantizer so that we can use the methods NUM2BIN, NUM2INT and
% NUM2HEX.
DataMode = 'fixed';
if ~s.Signed,
	DataMode = 'ufixed';
end
q = quantizer([s.CoeffWordLength, s.NumFracLength],'mode',DataMode);

switch radix,
case 2,   HqCoeff = num2bin(q,qcoeffs);
case 10,  HqCoeff = num2str(num2int(q,qcoeffs')); 
case 16,  HqCoeff = num2hex(q,qcoeffs);
otherwise
    error(message('signal:dfilt:dffir:coewrite:InvalidParam'));
end

% Display the Coefficients
fprintf(fid,'%s','CoefData = '); % 'CoefData' is a .COE file keyword
[r,~] = size(HqCoeff);
for i = 1:r-1;
    fprintf(fid,'%s,\n', HqCoeff(i,:));   
end

% Add the semicolon to the last coefficient
fprintf(fid,'%s;\n',HqCoeff(r,:)); 

fclose(fid);

% Launch the MATLAB editor (to display the generated file)
edit(file);


%------------------------------------------------------------------------
function [filename, pathname] = lcluiputfile(file,dlgStr)
% Local UIPUTFILE: Return an empty string for the "Cancel" case

[filename, pathname] = uiputfile(file,dlgStr);

% filename is 0 if "Cancel" was clicked
if filename == 0, filename = ''; end

% [EOF]
