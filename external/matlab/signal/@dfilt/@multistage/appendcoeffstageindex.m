function coeff = appendcoeffstageindex(this,coeff,index)
%APPENDCOEFFSTAGEINDEX Append stage index to coefficient names

%   Copyright 2009 The MathWorks, Inc.

if ~isstruct(coeff)
    for k = 1:length(coeff)
        % append stage number after the coefficient names
        coeff{k} = sprintf('%s%s%s',coeff{k},'_',index);
    end
else
    % If coefficient name is a structure, such as multistage or
    % iirmultirate filter.
    fd = fields(coeff);
    for stg = 1:length(fd)
        fieldname = fd{stg};
        stagecoeff = coeff.(fieldname);
        newindex = sprintf('%s%s%d',index,'_',stg);
        stagecoeff = appendcoeffstageindex(this,stagecoeff,newindex);
        coeff.(fieldname) = stagecoeff;
    end
end


% [EOF]
