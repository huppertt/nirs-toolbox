function fileProps = parseDesignfiltCallCodeFromFile(callingFileName,callingFileLine)
%parseDesignfiltCallCodeFromFile Parse designfilt calls when they come from
%a file or script.

%   Copyright 2013 The MathWorks, Inc.

strBuff = StringWriter;
strBuff.readfile(callingFileName);
strBuffLines = cellstr(strBuff);

isDesignFiltCall = false;
fileProps.InitialCallingLine = callingFileLine;

% Find beginning and end of call to designfilt

% Look for relevant code below the callingFileLine
for idx = callingFileLine : numel(strBuffLines)
  
  currentLine = strBuffLines{idx};
  
  if ~isempty(strfind(currentLine,'designfilt'))
    isDesignFiltCall = true;
  end
  
  contIdx = strfind(currentLine,'...');
  
  if isempty(contIdx)
    fileProps.FinalCallingLine = idx;
    break;
  end
end

if ~isDesignFiltCall
  error(message('signal:designfilt:InvalidLineOfCode',callingFileName));
end

% Make sure wo do not have relevant code above the callingFileLine
isDone = false;
N = fileProps.InitialCallingLine - 1;
while ~isDone && N > 0
  currentLine = strBuffLines{N};
  contIdx = strfind(currentLine,'...');
  if isempty(contIdx)
    isDone = true;    
  else
    N = N - 1;
  end
  fileProps.InitialCallingLine = N + 1;
end

codeLines = strBuffLines(fileProps.InitialCallingLine:fileProps.FinalCallingLine);
isDone = false;
while ~isDone
  % Parse code using MTREE
  [outputVarName, errFlag, codeCell, fcallLineIdx, indentLength, inputParamInfo] = ...
    signal.internal.parseDesignfiltCallCode(codeLines);
  
  if strcmp(errFlag,'incompleteCode') && numel(strBuffLines) > fileProps.FinalCallingLine
    fileProps.FinalCallingLine = fileProps.FinalCallingLine + 1;
    codeLines = strBuffLines(fileProps.InitialCallingLine:fileProps.FinalCallingLine);
  else
    isDone = true;
    if strcmp(errFlag,'incompleteCode')
      error(message('signal:designfilt:CannotParseCode'));
    end
  end
end

if islogical(errFlag) && errFlag
  error(message('signal:designfilt:InvalidLineOfCode',callingFileName));
end

[~,fName,ext] = fileparts(callingFileName);

fileProps.OutputVarName = outputVarName;
fileProps.HasEqualSign  = ~isempty(outputVarName);
fileProps.ShortFileName = [fName ext];
fileProps.CodeCell = codeCell;
fileProps.FcnCallLineIdx = fcallLineIdx;
fileProps.IndentLength = indentLength;
fileProps.InputParamNames = inputParamInfo.Names;
fileProps.InputParamValues = inputParamInfo.Values;
fileProps.InputParamValuesAreStrings = inputParamInfo.ValuesAreStrs;
fileProps.ErrorFlag = errFlag;

end