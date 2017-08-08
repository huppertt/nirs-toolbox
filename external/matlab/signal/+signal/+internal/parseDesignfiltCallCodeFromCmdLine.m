function fileProps = parseDesignfiltCallCodeFromCmdLine(cmdSet,type)
%parseDesignfiltCallCodeFromCmdLine Parse designfilt code coming from the
%command line.

% cmdSet can be a cell array with code lines, a java.lang.String[]
% object, or a string

% type ca be {'outputVar'} or 'all' - when set to outputVar, only the
% output var name is parsed, otherwise inputs, code, and output var names
% are parsed.

%   Copyright 2013 The MathWorks, Inc.


codeComesFromString = ischar(cmdSet);

fileProps.ErrorFlag = false;

if nargin < 2
  type = 'outputVar';
end

if codeComesFromString
  % Code was read from the comman line directly. Look for the last >>
  % prompt and get code from there and on. 
  
  idx = strfind(cmdSet,'>>');
  if isempty(idx)
    fileProps.ErrorFlag = false;
    return;
  elseif numel(idx) > 1 && idx(end)+2 == numel(cmdSet)
    % When calling this from test cases we have an extra set of prompt
    % symbols >> without any code after them. Remove them.
    idx = idx(1:end-1);
    cmdSet = cmdSet(1:end-3); % remove empty >> line
  end
  
  codeLines = cmdSet(idx(end)+2:end);    
  
  if isempty(codeLines)
    fileProps.ErrorFlag = false;
    return;
  end

  isDone = false;
  idx = 0;
  while ~isDone
    
    % Code must eventually be valid, but we might need to remove some lines
    % from the bottom to top until this happens. 
    T = mtree(codeLines);
    
    if isempty(mtfind(T, 'Kind', 'ERR'))
      isDone = true;
    else      
      strCell = regexp(codeLines, '[\f\n\r]', 'split');
      if numel(strCell) == 1
         error(message('signal:designfilt:CannotParseCode'));
      else
        idx = idx + 1;
        codeLines = getString(strCell,1,numel(strCell)-idx);
      end
    end
  end 
  
else
  N = numel(cmdSet);
  
  designfiltIdx = [];
  desigfiltCallLine = [];
  % Find the place where designfilt is called
  for idx = N:-1:1
    cmdStr = getCodeLine(cmdSet,idx);
    designfiltIdx = strfind(cmdStr,'designfilt');
    if isempty(designfiltIdx)
      if (idx - 1) > 0 && isempty(strfind(getCodeLine(cmdSet,idx-1),'...'))
        break;
      end
    else
      desigfiltCallLine = idx;
      break;
    end
  end
  
  if isempty(designfiltIdx)
    % We did not find a designfilt call so return without parsing parameters
    fileProps.ErrorFlag = true;
    return;
  end
  
  % Need to find begining and end of the command
  isDone = isempty(designfiltIdx);
  N = desigfiltCallLine-1;
  initLineIdx = N+1;
  while ~isDone && (N > 0)
    % Look for cont dots above the call line
    cmdStr = getCodeLine(cmdSet,N);
    if isempty(strfind(cmdStr,'...'))
      initLineIdx = N+1;
      isDone = true;
    else
      N = N-1;
    end
    if N < 1
      initLineIdx = 1;
    end
  end
  
  isDone = isempty(designfiltIdx);
  N = desigfiltCallLine;
  finalLineIdx = N;
  while ~isDone && (N <= numel(cmdSet))
    % Look for cont dots below the call line
    cmdStr = getCodeLine(cmdSet,N);
    if isempty(strfind(cmdStr,'...'))
      finalLineIdx = N;
      isDone = true;
    else
      N = N+1;
    end
  end
  
  [~,codeLines] = getString(cmdSet,initLineIdx,finalLineIdx);
end

isDone = false;
while ~isDone
  % Parse code using MTREE
  
  if strcmp(type,'outputVar')
    [outputVarName,errFlag] = signal.internal.parseDesignfiltCallCode(codeLines);
  else
    [outputVarName,errFlag,codeCell, fcallLineIdx, indentLength, inputParamInfo] = ...
      signal.internal.parseDesignfiltCallCode(codeLines);
  end
  
  if ~codeComesFromString && strcmp(errFlag,'incompleteCode') && numel(cmdSet) > finalLineIdx
    finalLineIdx = finalLineIdx + 1;
    [~,codeLines] = getString(cmdSet,initLineIdx,finalLineIdx);
  else
    isDone = true;
    if strcmp(errFlag,'incompleteCode')
      error(message('signal:designfilt:CannotParseCode'));
    end
  end
end

fileProps.OutputVarName = outputVarName;
fileProps.ErrorFlag = errFlag;
if strcmp(type,'all')
  fileProps.CodeCell = codeCell;
  fileProps.FcnCallLineIdx = fcallLineIdx;
  fileProps.IndentLength = indentLength;
  fileProps.InputParamNames = inputParamInfo.Names;
  fileProps.InputParamValues = inputParamInfo.Values;
  fileProps.InputParamValuesAreStrings = inputParamInfo.ValuesAreStrs;  
end

%--------------------------------------------------------------------------
function str = getCodeLine(cmdSet,N)

if isa(cmdSet,'java.lang.String[]')
  str = toCharArray(cmdSet(N)).';
else
  str = cmdSet{N};
end
%--------------------------------------------------------------------------
function [str,strCell] = getString(cmdSet,initLineIdx,finalLineIdx)

cmdSet = cmdSet(initLineIdx:finalLineIdx);

for idx = 1:numel(cmdSet)
  
  if isa(cmdSet,'java.lang.String[]') || isa(cmdSet,'java.lang.String')
    strCell{idx} = toCharArray(cmdSet(idx)).'; %#ok<*AGROW>
  else
    strCell{idx} = cmdSet{idx};
  end
  
  % Concatenate a string, do not use char since we need a one dimensional
  % string to pass to mtree
  if idx == 1
    str = strCell{idx};
  else
    str = sprintf('%s\n%s',str,strCell{idx});
  end
end


