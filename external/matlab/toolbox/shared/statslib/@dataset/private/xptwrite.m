function xptwrite(data,filename,varargin)
%XPTWRITE Writes a dataset array to a SAS XPORT format file.
%
%   XPTWRITE(DATA) writes the data stored in dataset array DATA to a SAS
%   XPORT format file. The function displays a dialog box for selecting the
%   name and location of the file. When DATA has observation names XPTWRITE
%   writes the observation names to a variable called OBSNAMES.
%
%   XPTWRITE(DATA, FILENAME) writes the data stored in dataset array DATA
%   to a SAS XPORT format file, FILENAME.
%
%   XPTWRITE(...,'WriteObsNames',FALSE) does not write the observation
%   names to the text file.
%
%   Note: XPTWRITE only supports scalar valued observations.
%
%   Example:
%
%      % Write data to a SAS XPORT format dataset:
%      patients = dataset('file','hospital.dat',...
%                   'delimiter',',',...
%                   'ReadObsNames',true);
%      xptwrite(patients,'patients.xpt')
%
%   SAS is a registered trademarks of SAS Institute Inc.
%
%   See also DATASET, XPTREAD.

%   Copyright 2009-2013 The MathWorks, Inc.


% The XPORT format specification can be found here:
% http://support.sas.com/techsup/technote/ts140.html

% Thanks to Cleve Moler for the ieee2ibm routines.

% Process input parameter name/value pairs
pnames = {'writeobsnames' 'datasetname'};
dflts =  { true              'MATLABDS'};
[writeObsFlag,dsName] = ...
    dataset.parseArgs(pnames, dflts, varargin{:});

if ~isa(data,'dataset')
    error(message('stats:dataset:xptwrite:DataNotDataset'));
end

% Treat observation names as a variable
if writeObsFlag
    ObsNames = get(data,'ObsNames');
    if ~isempty(ObsNames)
        data = [dataset({get(data,'ObsNames'),'OBSNAMES'}), data];
    end
end

if nargin<2 || isempty(filename)
    [F,P]=uiputfile('*.xpt');
    if isequal(F,0)
        return
    end
    filename = [P,F];
end

% Check that the extension is .xpt and add it if not
[~,~,fx] = fileparts(filename);
if ~strcmpi(fx,'.xpt')
    filename = [filename '.xpt'];
end

blanks8 = blanks(8);
[fid,theMessage] = fopen(filename,'w','b');
if fid == -1
    error(message('stats:dataset:xptwrite:FileOpenError', filename, theMessage));
end

% First line should be a standard library header file line
firstHeaderLine = 'HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!000000000000000000000000000000  ';

fwrite(fid,firstHeaderLine,'char');

sasString = padString('SAS',8);
sasLibString = padString('SASLIB',8);
sasDataString = padString('SASDATA',8);
sasVersion = padString('0.0',8);
OSString = padString(computer,8,true);
createdString = padString(datestr(now,'ddmmmyy:HH:MM:SS'),16);
% The format says that the real header looks like this:

% aaaaaaaabbbbbbbbccccccccddddddddeeeeeeee                 ffffffffffffffff
% In this record:
% -- aaaaaaaa and bbbbbbbb specify 'SAS '
% -- cccccccc specifies 'SASLIB '.
% -- dddddddd specifies the version of the SAS(r) System under which the file was created.
% -- eeeeeeee specifies the operating system that creates the record.
% -- ffffffffffffffff specifies the date and time created, formatted as ddMMMyy:hh:mm:ss. Note
% that only a 2-digit year appears. If any program needs to read in this 2-digit year, be prepared
% to deal with dates in the 1900s or the 2000s.
% Another way to consider this record is as a C structure:
% struct REAL_HEADER {
% char sas_symbol[2][8];
% char saslib[8];
% char sasver[8];
% char sas_os[8];
% char blanks[24];
% char sas_create[16];
% };
nextHeaderLine = [sasString,sasString,sasLibString,sasVersion,OSString,blanks(24),createdString];
fwrite(fid,nextHeaderLine,'char');
% The next line should be the date that the file was modified

nextHeaderLine = [createdString,blanks(64)];
fwrite(fid,nextHeaderLine,'char');
% Next is a MEMBER HEADER
%HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!000000000000000001600000000140

blockSize = 140;
% Note that the 140 at the end corresponds to the number of bytes in each
% VAR block
nextHeaderLine = ...
    sprintf('HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!000000000000000001600000000%s  ',num2str(blockSize));
fwrite(fid,nextHeaderLine,'char');
% Then a DESCRIPTOR HEADER
nextHeaderLine = ...
    'HEADER RECORD*******DSCRPTR HEADER RECORD!!!!!!!000000000000000000000000000000  ';
fwrite(fid,nextHeaderLine,'char');

dataSetName = padString(upper(dsName),8,true);
nextHeaderLine = [sasString,dataSetName,sasDataString,sasVersion,OSString,blanks(24),createdString];
fwrite(fid,nextHeaderLine,'char');

nextHeaderLine = [createdString,blanks(16),padString(dataSetName,40),blanks8];
fwrite(fid,nextHeaderLine,'char');

theVars = data.Properties.VarNames;
theVarDescriptions = data.Properties.VarDescription;

% We need to truncate the variable names to 8 characters without creating
% duplicates names
theXPTVars = truncateVarNames(theVars,8);

numVariables = numel(theVars);

namestrHeader = sprintf('HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!000000%04d00000000000000000000  ',numVariables);

fwrite(fid,namestrHeader,'char');

% Write out thr header information
pad52 = blanks(52);
npos = 0;
isVarChar = false(numVariables,1);
isVarFP = false(numVariables,1);
varTypes = cell(numVariables,1);
varLengths = zeros(numVariables,1);
charLengths = zeros(numVariables,1);

for count = 1:numVariables
    % Next we need to decide what type we have and set up the appropriate sizes
    theData = data.data{count};
    theClass = class(theData);
    varTypes{count} = theClass;
    if isa(theData,'cell')
        ntype = 2;
        if ~iscellstr(theData)
            fclose(fid);
            delete(filename);
            error(message('stats:dataset:xptwrite:CellVar', theVars{ count }))
        end
        nlng = max(cellfun(@length,theData));
        charLengths(count) = nlng;
        varTypes{count} = 'char';
        isVarChar(count) = true;
    elseif isa(theData,'categorical')
        ntype = 2;
        nlng =  max(cellfun(@length,categories(theData)));
        charLengths(count) = nlng;
        varTypes{count} = 'char';
        isVarChar(count) = true;
    elseif isa(theData,'int8') || isa(theData,'uint8')
        ntype = 1;
        nlng = 1;
    elseif isa(theData,'int16') || isa(theData,'uint16')
        ntype = 1;
        nlng = 2;
    elseif isa(theData,'int32') || isa(theData,'uint32')
        ntype = 1;
        nlng = 4;
    else % Assume 8 bytes numeric
        ntype = 1;
        nlng = 8;
        isVarFP(count) = true;
    end
    varLengths(count) = nlng;
    fwrite(fid,ntype,'int16');
    nhfun = 0;
    fwrite(fid,nhfun,'int16');
    fwrite(fid,nlng,'int16');
    fwrite(fid,count,'int16');
    
    % Warn if the variable name was truncated
    if ~strcmp(theXPTVars{count},theVars{count})
        warning(message('stats:dataset:xptwrite:TruncatedVarName', theVars{ count }, theXPTVars{ count }))
    end
    
    nname = padString(theXPTVars{count},8,true);
    fwrite(fid,nname,'char');
    % If there is no varDescription, use the un-truncated variable name.
    if isempty(theVarDescriptions) || isempty(theVarDescriptions{count})
        label = padString(theVars{count},40,true);
    else
        label = padString(theVarDescriptions{count},40,true);
    end
    fwrite(fid,label,'char');
    
    nform = blanks8;
    nfl = 0;
    nfd = 0;
    nfj = 0;
    nfill = blanks(2);
    fwrite(fid,nform,'char');
    fwrite(fid,nfl,'int16');
    fwrite(fid,nfd,'int16');
    fwrite(fid,nfj,'int16');
    fwrite(fid,nfill,'char');
    
    niform =  blanks8;
    nifl = 0;
    nifd = 0;
    fwrite(fid,niform,'char');
    fwrite(fid,nifl,'int16');
    fwrite(fid,nifd,'int16');
    fwrite(fid,npos,'int32');
    % Add some padding
    fwrite(fid,pad52,'char');
    
    npos = npos+nlng;
end

% Pad to 80 char boundary
padFile(fid);

% Now write the observation header
obsHeader = 'HEADER RECORD*******OBS     HEADER RECORD!!!!!!!000000000000000000000000000000  ';
fwrite(fid,obsHeader,'char');

numObs = size(data,1);

% Set the missing type to be . type
nanVal = uint8([hex2dec('2e') 0 0 0 0 0 0 0]);

% Observations must be scalar valued
for outer = 1:numObs
    for inner = 1:numVariables
        if isVarChar(inner)
            fwrite(fid,padString(char(data.data{inner}(outer)),charLengths(inner)),'char');
        else
            theVal = data.data{inner}(outer,:);
            if ~isscalar(theVal)
                fclose(fid);
                delete(filename);
                error(message('stats:dataset:xptwrite:ObsNotScalar', theVars{ inner }))
            end
            if isVarFP(inner) % Convert single to double before writing out as IBM format
                if ~isnan(theVal)
                    theVal = ieee2ibm(double(theVal));
                else
                    theVal = nanVal;
                end
                fwrite(fid,theVal,'uint8');
            else
                fwrite(fid,theVal,varTypes{inner});
            end
        end
    end
end

% Pad to 80 char boundary
padFile(fid);

% Close the file
fclose(fid);

function padFile(fid)
% Pad if necessary so that file contains a multiple of 80 chars

curpos = ftell(fid);
nfrag = rem(curpos,80);
if nfrag>0
    fwrite(fid,blanks(80-nfrag),'char');
end

function str = padString(inStr,paddedLen,truncFlag,fid,filename)
% The format requires strings to be padded with spaces of fixed length.
strLen = numel(inStr);
if nargin >2 && truncFlag
    inStr = inStr(1:min(strLen,paddedLen));
end
str = blanks(paddedLen);
try
    str(1:numel(inStr)) = inStr;
catch theException
    fclose(fid);
    delete(filename);
    if strLen>paddedLen
        error(message('stats:dataset:xptwrite:StrTooLong'));
    else
        rethrow(theException);
    end
end

function truncatedNames = truncateVarNames(names,maxLength,fid,filename)
% Truncates variables names to 8 characters and attempts to make unique
% names in cases where one of more names are truncated to the same
% substring. This is not foolproof so throw an error if we run into
% trouble.
if nargin < 2
    maxLength = 8;
end
numNames = numel(names);
truncatedNames = cellfun(@(x)x(1:min(maxLength,length(x))),names(:),'UniformOutput',false);
% see if they are unique
[sortedNames,~,j] = unique(truncatedNames);

% if there are no duplicates then we are done
if numel(sortedNames) == numNames
    return
end

% Try to get rid of duplicates by putting numerals on the end.
numDups = accumarray(j,1);
for theName = find(numDups>1)
    numericText = int2str((1:numDups(theName))');
    spaceNeeded = size(numericText(1,:),1);
    newNames = repmat(sortedNames(theName),numDups(theName),1);
    for count = 1:numDups(theName)
        newNames{count} = [newNames{count}(1:maxLength-spaceNeeded) strrep(numericText(count,:),' ','0')];
    end
    truncatedNames(j==theName) = newNames;
end

% Check that we didn't introduce a new problem
sortedNames = unique(truncatedNames);
if numel(sortedNames) ~= numNames
    % Maybe instead of erroring we could just spit out var1,var2,...
    fclose(fid);
    delete(filename);
    error(message('stats:dataset:xptwrite:CannotTruncateNames'))
end

function b = ieee2ibm(d)
% b = ieee2ibm(d)
% Convert an IEEE double, as used by MATLAB, to IBM 360 floating point
%    format, as used in SAS XPORT files.
% Input: a solitary double.
% Missing values are represented by denormal doubles with the code for
%    the missing value in the last byte.
% Output: 1-by-8 vector of "flints", integer-valued doubles in the
%    range 0 <= b(j) < 256.
% Ex:  d = 0, b = [0 0 0 0 0 0 0 0]
%      d = 1, b = [65 16 0 0 0 0 0 0]
%      d = 180.3, b = [66 180 76 204 204 204 204 208]
%      d = pow2(real('A'),-1074), b = [65 0 0 0 0 0 0 0]

% Thanks to Clever Moler for this routine.

% d = (+|-)16^(b(1)-64)*(b(2)/16^2 + b(3)/16^4 + ... + b(8)/16^14)

if d == 0
    b = zeros(1,8);
elseif abs(d) < pow2(256,-1074)  % Missing value
    b = [pow2(d,1074) zeros(1,7)];
elseif abs(d) < 2^(-256)   % IBM 360 underflow
    b = zeros(1,8);
elseif abs(d) >= 2^252     % IBM 360 overflow
    b = [Inf zeros(1,7)];
elseif isnan(d)
    b = [d zeros(1,7)];
else
    s = d < 0;              % Sign
    d = abs(d);
    [f,e2] = log2(d);       % Exponent and fraction, base 2,  d = f*2^e2
    e = floor(e2/4)+1;
    f = f*2^(e2-4*e);       % Exponent and fraction, base 16.
    f = 256*f;
    if f < 16
        f = 16*f;
        e = e-1;
    end
    b = zeros(1,8);
    b(1) = e + 64 + 128*s;  % Biased exponent and sign bit.
    for k = 2:8             % Peel off digits.
        b(k) = floor(f);
        f = 256*(f - b(k));
    end
end
