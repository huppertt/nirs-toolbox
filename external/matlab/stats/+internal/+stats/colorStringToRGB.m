function clrRGB = colorStringToRGB(clrString)
%COLORSTRINGTORGB converts a color string or cell strings to RGB valued
%matrix. If clrString is a numerical matrix, COLORSTRINGTORGB checks if it
%has 3 columns. Each row of the output corresponds to each element of the
%input.
%
%   Examples:
%   >> clrRGB = internal.stats.colorStringToRGB('rgb')
%   clrRGB =
%      1     0     0
%      0     1     0
%      0     0     1
%   >> clrRGB = internal.stats.colorStringToRGB({'r','blue','m'})
%   clrRGB =
%      1     0     0
%      0     0     1
%      1     0     1
%   >> clrRGB = internal.stats.colorStringToRGB({'red'})
%   clrRGB = 1     0     0

%   Copyright 2012 The MathWorks, Inc.



ColorStrShort = {'y','m','c','r','g','b','w','k'};
ColorStrLong = {'yellow','magenta','cyan','red','green','blue',...
    'white','black'};

if ischar(clrString)
    if isrow(clrString)
        clrString = cellstr(clrString')';
    else
        clrString = cellstr(clrString)';
    end
elseif iscellstr(clrString)
    if ~isvector(clrString)
        clrString = clrString';
        clrString = clrString(:);
    end   
elseif isnumeric(clrString)
    if any(clrString(:)<0)||any(clrString(:)>1)
        error(message('stats:internal:colorStringToRGB:BadColorValue'));
    end
    if  size(clrString,2) == 3
        clrRGB = clrString;
        return;
    elseif isequal(size(clrString),[3,1])
        clrRGB = clrString';
        return;
    else
        error(message('stats:internal:colorStringToRGB:ValueMustBe3ElementVector'));
    end   
end

numClrs = numel(clrString);
clrRGB = zeros(numClrs,3);

for i = 1:numClrs
    clr = clrString{i};
    clr = internal.stats.getParamVal(clr,[ColorStrShort,ColorStrLong],...
        '''Color''');
    
    switch clr        
        case {'r' 'red'}
            clrRGB(i,:) = [1 0 0];
        case {'g' 'green'}
            clrRGB(i,:) = [0 1 0];
        case {'b','blue'}
            clrRGB(i,:) = [0 0 1];
        case {'w' 'white'}
            clrRGB(i,:) = [1 1 1];
        case {'k' 'black'}
            clrRGB(i,:) = [0 0 0];
        case {'y' 'yellow'}
            clrRGB(i,:) = [1 1 0];
        case {'m' 'magenta'}
            clrRGB(i,:) = [1 0 1];
        case {'c' 'cyan'}
            clrRGB(i,:) = [0 1 1];
    end
end
 
end