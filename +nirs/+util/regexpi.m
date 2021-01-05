function START = regexpi(STR,EXPRESSION,forceCellOutput)

%Older matlab versions don't support the 'forceCellOutput' flag

if( ~verLessThan('Matlab','9.6'))
    if(nargin==3)
        START = regexpi(STR,EXPRESSION,forceCellOutput);
    else
        START = regexpi(STR,EXPRESSION);
    end
    return
else
     START = regexpi(STR,EXPRESSION);
     if(nargin==3 & strcmp(lower(forceCellOutput),'forceCellOutput'))
         if(~iscell(START))
             START={START};
         end
     end
end


