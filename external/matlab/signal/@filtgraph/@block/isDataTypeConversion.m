function flag = isDataTypeConversion(this)
% Check if the block type is data type conversion.  
% Note that
%   cast
%   caststage
%   convert
%   convertio
% are all data type conversion blocks.

%   Copyright 2008 The MathWorks, Inc.

if (strcmpi(this.blocktype,'convert') || strcmpi(this.blocktype,'convertio') ...
        || strcmpi(this.blocktype,'cast') || strcmpi(this.blocktype,'caststage'))
    flag = true;
else
    flag = false;
end
            