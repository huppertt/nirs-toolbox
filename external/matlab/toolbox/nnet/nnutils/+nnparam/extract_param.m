function [args,param,nargs] = extract_param(allArgs,defaultParam)

% Copyright 2010-2011 The MathWorks, Inc.

fn = fieldnames(defaultParam);
numArgs = length(allArgs);

% TODO - check parameter values

% Case 0: No arguments or parameters
if numArgs == 0
  args = {};
  param = defaultParam;
  nargs = 0;
  return;
end

% Case 1: Last argument is a parameter structure or object
lastArg = allArgs{end};
if isstruct(lastArg) || isa(lastArg,'nnetParam')
  param = defaultParam;
  foundField = false;
  for i=1:length(fn)
    fni = fn{i};
    if isfield(lastArg,fni)
      foundField = true;
      param.(fni) = lastArg.(fni);
    end
  end
  if foundField || isempty(fieldnames(lastArg))
    args = allArgs(1:(end-1));
    nargs = length(args);
  else
    args = allArgs;
    nargs = length(args);
  end
  
% Case 2: Last arguments variable number of parameter name/value pairs
elseif (numArgs>=2) && nntype.string('isa',allArgs{numArgs-1})
  param = defaultParam;
  fnlower = lower(fn);
  lastarg = length(allArgs);
  for i=(numArgs-1):-2:1
    argi = allArgs{i};
    if nntype.string('isa',argi)
      j = nnstring.first_match(lower(argi),fnlower);
      if ~isempty(j)
        param.(fn{j}) = allArgs{i+1};
        lastarg = i-1;
      else
        args = allArgs(1:lastarg);
        nargs = length(args);
        return
      end
    else
      args = allArgs(1:lastarg);
      nargs = length(args);
      return;
    end
  end
  args = allArgs(1:lastarg);
  nargs = length(args);
  
% Case 3: No parameters provided
else
  args = allArgs;
  param = defaultParam;
  nargs = length(args);
end
