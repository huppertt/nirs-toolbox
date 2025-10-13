function info=get_toolbox_info()

folder=which('nirs.core.Data');

gitinfo=fullfile(folder(1:strfind(folder,'+nirs')-1),'VERSION.txt');
warning('off','MATLAB:deblank:NonStringInput');

info=struct;
    info.location=folder(1:strfind(folder,'+nirs')-1);

if(exist(gitinfo,'file')==2)
    % The code is under Git control.  Doing this through gitstatus would be
    line=1; 
    fid=fopen(gitinfo,'r');
    while(line~=-1)
        line=fgetl(fid);
        if(line~=1)
            [a,b]=strtok(line,':');
            a=strrep(a,' ','_');
            a=strrep(a,')','');
            a=strrep(a,'(','');
            try;
                 b=strrep(b,': ','');
                 if(~isempty(b) & ~strcmp(b,':'))
                info=setfield(info,a,b);
                end
            end
        end
    end
end

