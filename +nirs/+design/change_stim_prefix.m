function data=change_stim_prefix(data,prefix,aspostfix)
% prefix=@(data)data.demographics('session')
% data=change_stim_prefix(data,prefix)

if(nargin<3)
    aspostfix=false;
end

for idx=1:length(data)
    if(isa(prefix,'function_handle'))
        prefix2=feval(prefix,data(idx));
        if(aspostfix)
            prefix2=['_' prefix2];
        else
            prefix2=[prefix2 '_'];
        end
    else
        prefix2=prefix;
    end

    if(isa(data(idx),'nirs.core.Data'))
        keys=data(idx).stimulus.keys;
    elseif(isa(data(idx),'nirs.core.ChannelStats'))
        keys=data(idx).conditions;
    end
    job=nirs.modules.RenameStims;
    for i=1:length(keys)
        job.listOfChanges{i,1}=keys{i};
        if(aspostfix)
            job.listOfChanges{i,2}=[keys{i} prefix2];
        else
            job.listOfChanges{i,2}=[prefix2 keys{i}];
        end
    end
    data(idx)=job.run(data(idx));
end

        