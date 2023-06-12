function strnew=fix_strings(strold,extended)
%This function removes any non-alphabetic charectors from a string

if(~ischar(strold))
    strnew=strold;
    return;
end

if(nargin<2)
    extended=true;
end

if(extended)
    allowed='abcdefghijklmnopqrtsuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 /\-_.?!$#*%<>,~[]{}|';
else
    allowed='abcdefghijklmnopqrtsuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 /.\';
end


strnew=strtrim(strold(ismember(double(strold),double(allowed))));
