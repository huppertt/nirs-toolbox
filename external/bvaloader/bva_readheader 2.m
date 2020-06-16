function [fs label meta] = bva_readheader(varargin)

%bva_readheader - read BVA Header File and extract sampling rate and channels
%
% function [fs label meta] = bva_readheader(file)
%
% Read Brain Vision Header File and return sampling rate, channel
% labels and various other info to caller. 
%
% Input:
%	file = Brain Vision Data HeaderFile (.vhdr)
%
% Output:
%	fs = sampling rate
%	label = electrode labels
%	meta = other info (required for bva_loadeeg.m)
%
% requires: 
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

%% debug settings
debug = 0;
if debug;warning('on','all');else warning('off','all');end

%% check number of input arguments
error(nargchk(0,1,nargin))

%% check number of out arguments
error(nargoutchk(0,4,nargout))

%% check && assign input
varargin{2} = [];

if ~isempty(varargin{1}),
	file = varargin{1};
else
	%use open dialog
   [filename, pathname] = uigetfile( {'*.vhdr'},  'Choose Header File');

	if filename == 0
		error('ERRP:IO:bva_readheader:NoHeaderFile','No Header File Provided')
	else
		file = fullfile(pathname, filename);
	end
end

%open fp (for reading)
fp = fopen(file);

if fp == -1
	error('ERRP:IO:bva_readheader:FileNotReadable','Couldn''t read header file')
end

% read the whole file as one cell array
raw={};
while ~feof(fp)
    raw = [raw; {fgetl(fp)}];
end
fclose(fp);

% Remove comments and empty lines
raw(strmatch(';', raw)) = [];
raw(cellfun('isempty', raw) == true) = [];

% Find sections
sections = [strmatch('[', raw)' length(raw) + 1];
for section = 1:length(sections) - 1

    % Convert section name
    fieldname = lower(char(strread(raw{sections(section)}, '[%s', 'delimiter', ']')));
    fieldname(isspace(fieldname) == true) = [];

    % Fill structure with parameter value pairs
    switch fieldname
        case {'commoninfos' 'binaryinfos'}
            for line = sections(section) + 1:sections(section + 1) - 1
                [parameter, value] = strread(raw{line}, '%s%s', 'delimiter', '=');
                hdr.(fieldname).(char(parameter)) = char(value);
            end
        case {'channelinfos' 'coordinates' 'markerinfos'}
            for line = sections(section) + 1:sections(section + 1) - 1
                [parameter, value] = strread(raw{line}, '%s%s', 'delimiter', '=');
                hdr.(fieldname)(str2double(parameter{1}(3:end))) = value;
            end
        case 'comment'
            hdr.(fieldname) = raw(sections(section) + 1:sections(section + 1) - 1);
    end
end


% slice out what we need 

% 1st sampling rate
fs = 1000/ str2num(hdr.commoninfos.SamplingInterval) * 1000;

% 2nd label and chan and the scaleFactor for later on
label = {};
scale  = [];
for iChan = 1:size(hdr.channelinfos,2)
	[t scale(iChan)] =  strread(hdr.channelinfos{iChan},'%s%*d%f%*s','delimiter',',');
	label{iChan} = t{1};
end


% info for loading data
meta.eegFile = hdr.commoninfos.DataFile;
meta.DataFormat = hdr.commoninfos.DataFormat;
meta.DataType = hdr.binaryinfos.BinaryFormat;
meta.DataOrientation = hdr.commoninfos.DataOrientation;
meta.scaleFactor  = scale;

