function mysave(fn,data)

% mysave(fn,data)
%
% save function that will work with matlab compiler
% 
% fn is the filename to be saved
% data is the variable to be saved


dlim = '\t';
data = data';
[nrow,ncol]=size(data);
fid = OpenFile(fn,'w');

str = ['fprintf(fid, ', '''', '%g'];
if nrow>1
    str1 = [];
    for i = 2: nrow
        str1 = [str1, dlim, '%g'];
    end
end
if exist('str1', 'var')
    str = [str, str1, '\n', '''', ',data);'];
else
    str = [str, '\n', '''', ',data);'];
end

eval(str)

fclose(fid);
