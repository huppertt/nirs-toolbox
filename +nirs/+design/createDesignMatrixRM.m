function [X,names,offset,metadata] = createDesignMatrixRM( stimulus, t, basis, type )
% creates a table for use in the fitRM (repeated measures) model

if nargin < 4, type = ''; end

if(isa(stimulus,'nirs.core.Data'))
    if nargin < 3, type = ''; else;  type=basis; end
    if nargin < 2,
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    else
        basis=t;
    end
    
    t=stimulus.time;
    stimulus=stimulus.stimulus;
elseif nargin < 3,
    basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
end

if(isstr(basis))
    basis=Dictionary({'default'},{basis});
end


% do the additional metadata
keys=stimulus.keys;
names={};
cats={};
for i=1:length(keys)
    st=stimulus(keys{i});
    names={names{:} st.metadata.Properties.VariableNames{:}};
    
    n=st.metadata.Properties.VariableNames;
    for j=1:length(n)
        cats{end+1}=class(st.metadata.(n{j}));
        if(strcmp(class(st.metadata.(n{j})),'orfinal'))
            st.metadata.(n{j})=addlevel(st.metadata.(n{j}),'undefined');
        end
    end
    
end
[names,i]=unique(names);
cats={cats{i}};

metadata=cell(length(t),length(names));
for i=1:length(t); 
    for j=1:length(names)
        if(strcmp(cats{j},'double'))
            metadata{i,j}=NaN;
        elseif(strcmp(cats{j},'categorical'))
            metadata{i,j}=categorical(NaN);
        elseif(strcmp(cats{j},'ordinal'))
            metadata{i,j}='undefined';
        end
    end
end

for i=1:length(keys)
    for j=1:length(st.onset)
        st2=st;
        st2.onset=st.onset(j);
        st2.dur=st.dur(j);
        st2.amp=st.amp(j);
        D=Dictionary({keys{i}},{st2});
        XX=nirs.design.createDesignMatrix(D, t, basis);
       [lst,~]=find(XX~=0);
       for k=1:length(names)
            if(ismember(names{i},st.metadata.Properties.VariableNames))
                for l=1:length(lst)
                    metadata{lst(l),k}=st.metadata.(names{k})(j);
                end
            end
       end
    end
end
metadata=cell2table(metadata,'VariableNames',names);
    
[X,names,offset]=nirs.design.createDesignMatrix( stimulus, t, basis, type );