function save_stim2excel(data,filename)
% This function saves a 


for i=1:length(data)
    if(~isempty(data(i).description))
        [~,sheetname]=fileparts(data(i).description);
    else
        sheetname=['File-' num2str(i)];
    end
    keys=data(i).stimulus.keys;
    Array={}; 
    for j=1:length(keys)
       
        st=data(i).stimulus(keys{j});
        if(isa(st,'nirs.design.StimulusEvents'))
            Array{1,j+3*(j-1)}=[keys{j} '@onset'];
            Array{1,j+1+3*(j-1)}=[keys{j} '@duration'];
            Array{1,j+2+3*(j-1)}=[keys{j} '@amplitude'];
            
            for k=1:length(st.onset)
                Array{1+k,j+3*(j-1)}=st.onset(k);
                Array{1+k,j+1+3*(j-1)}=st.dur(k);
                Array{1+k,j+2+3*(j-1)}=st.amp(k);
            end
        elseif(isa(st,'nirs.design.StimulusVector'))
              Array{1,j+3*(j-1)}=[keys{j} '@vector'];
              Array{1,j+1+3*(j-1)}=[keys{j} '@time'];
            
            for k=1:length(st.vector)
                Array{1+k,j+3*(j-1)}=st.vector(k);
                Array{1+k,j+1+3*(j-1)}=st.time(k);
            end
            
        end
    
        
    end
    for j=1:size(Array,2)
       
        for k=1:size(Array,1)
            if(isempty(Array{k,j}))
                Array{k,j}=NaN;
            end
        end
    end
    
    if(strcmp(filename,'clipboard'))
        nirs.util.copytable2clip(cell2table(Array));
        disp('copied to clipboard');
    else
    
    if(~isempty(sheetname))
        sheetname=['File-' num2str(i)];
    end
        
      nirs.util.write_xls(filename,Array,sheetname);
    end
end