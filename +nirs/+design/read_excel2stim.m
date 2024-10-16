function data=read_excel2stim(data,filename)


for id=1:length(data)
    
    if(strcmp(filename,'clipboard'))
        alldata=nirs.util.pasteclip2table;
        alldata=table2cell(alldata(2:end,:));   
        for i=2:size(alldata,1)
            for j=1:size(alldata,2)
                alldata{i,j}=str2num(alldata{i,j});
            end
        end
        
    else
    if(~isempty(data(id).description))
        [~,sheetname]=fileparts(data(id).description);
    else
        sheetname=['data-' num2str(id)];
    end
    try
        [~,~,alldata]=xlsread(filename,sheetname);
        if(any(ismember({alldata{1,:}},{'Array1','Array2','Array3'})))
            alldata=alldata(2:end,:);
        end
    catch
        if(strcmp(MException.last.identifier,'MATLAB:xlsread:WorksheetNotOpened'))
            disp(['No worksheet for: ' sheetname]);
        else
            rethrow(MException.last);
        end
    end
    end
    stimtypes={}; types={}; idx=[];
    for i=1:size(alldata,2)
        if(~isnan(alldata{1,i}))
            stimtypes{end+1}=alldata{1,i}(1:strfind(alldata{1,i},'@')-1);
            types{end+1}=alldata{1,i}(strfind(alldata{1,i},'@')+1:end);
            idx(end+1)=i;
        end
    end
    ustimtypes=unique(stimtypes);
    for i=1:length(ustimtypes)
        lst=find(ismember(stimtypes,ustimtypes{i}));
        
        if(ismember('onset',{types{lst}}))
            st=nirs.design.StimulusEvents;
            a1=find(ismember({types{lst}},'onset'));
            a2=find(ismember({types{lst}},'duration'));
            a3=find(ismember({types{lst}},'amplitude'));
            st.name=ustimtypes{i};
            lst=idx(lst);
            for j=2:size(alldata,1)
                if(~isnan(alldata{j,lst(a1)}) & ~isnan(alldata{j,lst(a2)}) &...
                        ~isnan(alldata{j,lst(a3)}))
                    st.onset(end+1,1)=alldata{j,lst(a1)};
                    st.dur(end+1,1)=alldata{j,lst(a2)};
                    st.amp(end+1,1)=alldata{j,lst(a3)};
                end
            end
         data(id).stimulus(ustimtypes{i})=st;   
            
            
        else
            st=nirs.design.StimulusVector;
        end
        
    end
    
    
    
end


