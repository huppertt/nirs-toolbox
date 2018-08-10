function t = readFromFile(filename,args)
% READFROMFILE Create a table by reading from a file.

%   Copyright 2012 The MathWorks, Inc.
try
    pnames = {'FileType'};
    dflts =  {       [] };
    [fileType,supplied,otherArgs] = matlab.internal.table.parseArgs(pnames, dflts, args{:});

    if ~supplied.FileType
        [~,~,fx] = fileparts(filename);
        switch fx
        case {'.txt' '.dat' '.csv'}, fileType = 'text';
        case {'.xls' '.xlsx' '.xlsb' '.xlsm' '.xltm' '.xltx' '.ods'}, fileType = 'spreadsheet';
        case '', fileType = 'text';
        otherwise
            error(message('MATLAB:readtable:UnrecognizedFileExtension',fx));
        end
    else
        fileTypes = {'text' 'spreadsheet'};
        itype = find(strncmpi(fileType,fileTypes,length(fileType)));
        if isempty(itype)
            error(message('MATLAB:readtable:UnrecognizedFileType',fileType));
        elseif ~isscalar(itype)
            error(message('MATLAB:readtable:AmbiguousFileType',fileType));
        end
    end

    % readTextFile and readXLSFile will add an extension if need be, no need to add one here.

    switch lower(fileType)
    case 'text'
        t = table.readTextFile(filename,otherArgs);
    case 'spreadsheet'
        t = table.readXLSFile(filename,otherArgs);
    end
catch ME
    throwAsCaller(ME)
end
