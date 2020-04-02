function nirfastroot = get_nirfast_root()
% figure out the root folder of nirfast
% Written by Hamid Ghadyani
% 2012, Dartmouth College

foo = which('nirfast');
if isempty(foo)
    error(' Can''t find nirfast main script.')
end
foo = regexp(fileparts(foo),filesep,'split');
nirfastroot = fullfile(foo{1:end-1});
if ~ispc
    nirfastroot = [filesep nirfastroot];
end