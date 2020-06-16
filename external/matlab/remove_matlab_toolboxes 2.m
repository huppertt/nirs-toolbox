%This script removes all matlab licensed toolboxes 
% for testing the dependencies of my toolbox

matlabroot=which('addpath');
matlabroot=matlabroot(1:strfind(matlabroot,'toolbox')-1);
p=path;


%keep toolboxs
keep={'local','eml','coder','matlab','simulink','compiler','javabuilder','stm'};
fold=dir(fullfile(matlabroot,'toolbox'));
warning('off','MATLAB:rmpath:DirNotFound');
for i=1:length(fold)
    if(~ismember(fold(i).name,keep) & isempty(strfind(fold(i).name,'.')))
        g=genpath(fullfile(matlabroot,'toolbox',fold(i).name));
        rmpath(g);
    end
end