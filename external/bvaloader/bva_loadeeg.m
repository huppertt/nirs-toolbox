function eeg = bva_loadeeg(hdrFile)

%bva_loadeeg - load continous eeg from  Brain Vision Data Exchange File 
%
% function [eeg] = bva_loadeeg(header)
%
% Load Brain Vision data exchange file data. 
% 
%
% Input:
%	file = Brain Vision Data HeaderFile (.vhdr)
%
% Output:
%	eeg = continuous EEG data
%
% requires: bva_readheader.m
%
% see also: ERRP/io
%

% Copyright (C) 2008-2012 Stefan Schinkel, HU Berlin
% http://people.physik.hu-berlin.de/~schinkel/
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% $Log$
if nargin < 1
	%use open dialog
   [fileName, pathName] = uigetfile( {'*.vhdr'},  'Choose Header File');

	if fileName == 0
		error('ERRP:IO:bva_readheader:NoHeaderFile','No Header File Provided')
	else
		hdrFile = fullfile(pathName, fileName);
	end
else
	%make sure we have full aht
	[pathName fileName fileExt] = fileparts(hdrFile);
	
	if isempty(pathName); pathName = pwd;end
	hdrFile  = fullfile(pathName, [fileName fileExt]);
end

hdrFile;


% read meta data
[fs label meta ] = bva_readheader(hdrFile);
nChans = numel(label);

% assemble eegFile name and check if exists
if exist(fullfile(pathName, meta.eegFile),'file')
	eegFile = fullfile(pathName, meta.eegFile);
elseif exist(fullfile(pathName, lower( meta.eegFile )),'file');
	eegFile = fullfile(pathName, lower( meta.eegFile ));
else
	error('ERRP:io:bva_loadeeg','Couldn''t find find .eeg file');
end


% make sure we read binary data
if strcmpi(meta.DataFormat, 'ascii')
	error('Reading ASCII is not implemented. Yet.')  
end;

% extract data format, determins bytesPerSample
switch lower(meta.DataType)
    case 'int_16',        binformat = 'int16'; bytesPerSample = 2;
    case 'uint_16',       binformat = 'uint16'; bytesPerSample = 2;
    case 'ieee_float_32', binformat = 'float32'; bytesPerSample = 4;
    otherwise, error('Unsupported binary format');
end

% open file
fp = fopen(eegFile,'r');
%seek to file end and get position in byte (total bytes in file)
fseek(fp, 0, 'eof');
totalBytes =  ftell(fp);

% number of Frames 
nFrames =  totalBytes / (bytesPerSample * nChans);

%pre-allocate memory
eeg = single( zeros(nChans,nFrames) );

% Read data
switch lower(meta.DataOrientation)
    case 'multiplexed'
		frewind(fp);                
		eeg = fread(fp, [nChans, nFrames], [binformat '=>float32']);
		fclose(fp);
    case 'vectorized'
		error('Reading vectorized binary is not Not implemented')
    otherwise
		error('Not implemented')
end

%scale data
eeg = eeg .* repmat(meta.scaleFactor',1,size(eeg,2));

