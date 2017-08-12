function helpline = get_mhelp_title(f)

% Copyright 2010-2014 The MathWorks, Inc.

path = which(f);
file = fopen(path,'r','n','UTF-8');

filesepInd = find(f == filesep,1,'last');
if ~isempty(filesepInd)
  f = f((filesepInd+1):end);
end
lenf = length(f);

if (file == -1)
  helpline = '';
  return
end

while (true)
  line = fgetl(file);
  if isnumeric(line) && (numel(line)==1) && (line == -1)
    helpline = ''; break;
  end
  if ischar(line) && (length(line) > 1) && (line(1) == '%')
    helpline = line(2:end);
    if (length(helpline) >= lenf) && all(upper(f)==upper(helpline(1:lenf)))
      helpline = helpline((lenf+1):end);
    end
    while ~isempty(helpline) && (helpline(1) == ' ')
      helpline = helpline(2:end);
    end
    break;
  end
end
fclose(file);
