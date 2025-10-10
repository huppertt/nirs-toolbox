function info=get_toolbox_info()

folder=which('nirs.core.Data');

gitinfo=fullfile(folder(1:strfind(folder,'+nirs')-1),'.git','logs','HEAD');

if(exist(gitinfo,'file')==2)
    % The code is under Git control.  Doing this through gitstatus would be
    
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    tbl=readtable(gitinfo,'FileType','text');
    v=tbl.(tbl.Properties.VariableNames{1}){end};
    [a,v]=strtok(v,' ');
    gitID=strtrim(strtok(v,' '));

    info.location=folder(1:strfind(folder,'+nirs')-1);
    info.GitHub_version=gitID;
else
    info.location=folder(1:strfind(folder,'+nirs')-1);
    info.GitHub_version='unknown';
end

