function data = loadDotNirs( filenames,force )

if(nargin<2 || isempty(force))
    force=false;
end

    % if a single filename, put it in a cell
    if ischar( filenames )
        filenames = {filenames};
    end
    
  
    
    data = nirs.core.Data.empty;
    
    % iterate through cell array
    for iFile = 1:length(filenames)
        
        [p f]=fileparts(filenames{iFile});
        if(~isempty(dir(fullfile(p,[f '.wl1']))) & ~force)
            %disp(['Skipping ' filenames{iFile} ': NIRx data found in same folder']);
            continue;
        end
        disp(['Loading ' filenames{iFile}]);
        try
            % load data as a struct
            d = load( filenames{iFile}, '-mat' );
           
            
            d.d=nirs.util.fixnan(d.d);
            if(isfield(d,'s'))
                d.s=nirs.util.fixnan(d.s);
            end
            if(isfield(d,'aux'))
                d.aux=nirs.util.fixnan(d.aux);
            end
            
            % put into data class
            thisFile = nirs.core.Data();
            thisFile.description = filenames{iFile};
                
            nTime = length(d.t);

            % data
            if size(d.d,1) == nTime
                thisFile.data = d.d;
            elseif size(d.d,2) == nTime
                thisFile.data = d.d.';
            else
                error('Data length and time vector don''t match in size.')
            end
            
            if(~isfield(d.SD,'MeasList'))
                d.SD.MeasList=d.ml;
            end

            % time vector
            thisFile.time = d.t(:);

            % probe
            thisFile.probe = nirs.util.sd2probe( d.SD );

            % stimulus vector
            if isfield(d,'StimDesign')
                thisFile.stimulus = nirs.util.convertStimDesignStruct( d.StimDesign,d.t );
            else
                % This will handle the HOMER-2 nirs format
                stims = Dictionary();
                for idx=1:size(d.s,2)
                    
                    if isfield(d,'CondNames')
                        stimname = d.CondNames{idx};
                    else
                        stimname = ['stim_channel' num2str(idx)];
                    end

                    if(islogical(d.s)); d.s=d.s*1; end;
                    d.s(:,idx)=d.s(:,idx)./max(d.s(:,idx));

                    s = nirs.design.vector2event( d.t , d.s(:,idx) , stimname );

                    if(~isempty(s.onset))
                        stims(stimname)=s;
                    end
                        
                end
                
                thisFile.stimulus=stims;
                
            end

            % demographics for group level
            if isfield(d,'Demographics')
                for i = 1:length(d.Demographics)
                    name    = d.Demographics(i).name;
                    value   = d.Demographics(i).value;
                    thisFile.demographics(name) = value;
                end
            end
         
            if(isfield(d,'brainsight'))
                a=d.brainsight.acquisition.auxData;
                t=d.brainsight.acquisition.auxTime;
                 for i=1:size(a,2)
                    name{i,1}=['aux-' num2str(i)];
                end
                aux=nirs.core.GenericData(a,t,table(name,repmat({'aux'},length(name),1),'VariableNames',{'name','type'}));
                thisFile.auxillary('aux')=aux;
                thisFile.auxillary('brainsight')=d.brainsight;
            end
            
            if(isfield(d,'aux') || isfield(d,'aux10')) && (~isempty(d.aux))
                if(isfield(d,'aux'))
                    a=d.aux;
                else
                    a=d.aux10;
                end
                a=reshape(a,size(a,1),[]);
                for i=1:size(a,2)
                    name{i,1}=['aux-' num2str(i)];
                end
                aux=nirs.core.GenericData(a,d.t,table(name,repmat({'aux'},length(name),1),'VariableNames',{'name','type'}));
                aux.description=thisFile.description;
                thisFile.auxillary('aux')=aux;
            end
            
        % append to list of data
        data(end+1) = thisFile.sorted();
            
        catch err
            if(~exist('d') || ~isfield(d,'d') || isempty(d.d))
                 disp('Empty file found (skipping):');
                 disp(filenames{iFile});
            else
                warning(err.message)
            end
            
        end
        
        
    end


end

