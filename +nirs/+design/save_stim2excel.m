function save_stim2excel(data,filename)
% This function saves a 

root=fileparts(which('nirs.design.save_stim2excel'));

%% Initialisation of POI Libs
% Add Java POI Libs to matlab javapath
if exist('org.apache.poi.ss.usermodel.WorkbookFactory', 'class') ~= 8 ...
        || exist('org.apache.poi.hssf.usermodel.HSSFWorkbook', 'class') ~= 8 ...
        || exist('org.apache.poi.xssf.usermodel.XSSFWorkbook', 'class') ~= 8
    javaaddpath(fullfile(root,'private','poi-3.8-20120326.jar'));
    javaaddpath(fullfile(root,'private','poi-ooxml-3.8-20120326.jar'));
    javaaddpath(fullfile(root,'private','poi-ooxml-schemas-3.8-20120326.jar'));
    javaaddpath(fullfile(root,'private','xmlbeans-2.3.0.jar'));
    javaaddpath(fullfile(root,'private','dom4j-1.6.1.jar'));
    javaaddpath(fullfile(root,'private','stax-api-1.0.1.jar'));
end

for i=1:length(data)
    if(~isempty(data(i).description))
        [~,sheetname]=fileparts(data(i).description);
    else
        sheetname=['data-' num2str(i)];
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
    
    warning('off','xlwrite:AddSheet');
    xlwrite(filename,Array,sheetname);
end