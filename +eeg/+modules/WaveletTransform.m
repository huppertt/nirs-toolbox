classdef WaveletTransform < nirs.modules.AbstractModule
%% WaveletTransfrom - Preforms wavelet transfrom on EEG data.
% 
    
    properties
        wname;  % wavelet model name
        scale_smoothing; % add scale smoothing step
        frequencies;  % table of frequency values
        convert2power;
        normalize;
        singleoutput;  % outputs the values all into a single data structure
        nscales;
    end
    
    methods

        function obj =WaveletTransform( prevJob )
           obj.name = 'Wavelet Transform using CWT';
           if nargin > 0
               obj.prevJob = prevJob;
           end
           obj.wname='morl';
          
           lower=[1 6 7.5 9 12.5 25]';
           upper=[4 7 12.5 11 30 50]';
           name={'delta','theta','alpha','beta','mu','gamma'}';
           obj.scale_smoothing=true;
           obj.frequencies = table(lower,upper,name);
           obj.convert2power=true;
           obj.normalize=false;
           obj.singleoutput=false;
           obj.nscales=32;
        end
        
        function data = runThis( obj, data )
            dataOut=data(1);
            dataOut(:)=[];
            for i = 1:length(data)
                FS=data(i).Fs;
                
                if(isa(obj.frequencies,'table'))
                    
                    lstbad=[];
                    for j=1:height(obj.frequencies)
                        obj.frequencies.lower(j)=max(obj.frequencies.lower(j),1);
                        obj.frequencies.upper(j)=min(obj.frequencies.upper(j),FS/2);
                        if(obj.frequencies.lower(j)>FS/2)
                            lstbad=[lstbad; j];
                        end
                    end
                    obj.frequencies(lstbad,:)=[];
                    
                    
                    f=linspace(min(obj.frequencies.lower),max(obj.frequencies.upper),obj.nscales);
                    f=unique(max(f,1));
                    scales= centfrq(obj.wname)./(f* 1/FS);
                    f=scal2frq(scales,obj.wname,1/FS);
                else
                   scales= centfrq(obj.wname)./(obj.frequencies* 1/FS);
                   f=scal2frq(scales,obj.wname,1/FS);
                end
                NTW = fix(FS);
                NSW = fix(length(f)/10);
                
                disp('Progress');
                n=size(data(i).data,2);
                cnt=1;
                str='  0';
                fprintf('%s %',str(end-2:end));
                cfs=zeros(length(f),size(data(i).data,1),size(data(i).data,2));
                 for j=1:size(data(i).data,2)
                    cfs(:,:,j) = nirs.math.wavelet(data(i).data(:,j),FS,scales,obj.wname,obj.normalize,obj.scale_smoothing,NSW,NTW);
                    
                    str=['   ' num2str(round(100*cnt/n))];
                    fprintf('\b\b\b\b%s %',str(end-2:end));
                    cnt=cnt+1;
                 end
                disp('Computing frequency bands');
                
               
               
                tbl=table;
                ele=data(i).probe.link;
                n=size(data(i).data,2);
                
                if(isa(obj.frequencies,'table'))
                    for j=1:height(obj.frequencies)
                        lst=find(f>=obj.frequencies.lower(j) & f<=obj.frequencies.upper(j));
                        
                        if( obj.convert2power)
                            d(j,:,:)=sqrt(mean(cfs(lst,:,:).^2,1));
                        else
                            d(j,:,:)=mean(cfs(lst,:,:),1);
                        end
                        F{j}=[num2str(obj.frequencies.lower(j)) '-' ...
                            num2str(obj.frequencies.upper(j)) 'Hz'];
                    end
                     if(~obj.singleoutput)
                        for id=1:size(d,1)
                            dataOut(end+1)=data(i);
                            dataOut(end).data=squeeze(d(id,:,:));
                            dataOut(end).demographics('freq')=F{id};
                            
                        end
                    else
                        link=data(i).probe.link;
                        link=repmat(link,size(d,1),1);
                        link.type=reshape(repmat(strcat(repmat(cellstr('eeg:'),...
                            length(F),1),F'),1,...
                            height(data(i).probe.link))',[],1);
                        dataOut(end+1)=data(i);
                        dataOut(end).probe.link=link;
                        dataOut(end).data=reshape(permute(d,[2 3 1]),...
                            size(d,2),size(d,1)*size(d,3));
                    end
                        
                    
                else
                    if( obj.convert2power)
                        cfs=abs(cfs);
                   end
                    if(~obj.singleoutput)
                        for id=1:length(f)
                            dataOut(end+1)=data(i);
                            dataOut(end).data=squeeze(cfs(id,:,:));
                            dataOut(end).demographics('freq')=f(id);
                        end
                    else
                        link=data(i).probe.link;
                        link=repmat(link,length(f),1);
                        link.type=reshape(repmat(strcat(repmat(cellstr('eeg:'),length(f),1),...
                            strtrim(cellstr(num2str(f'))),...
                            repmat(cellstr('Hz'),length(f),1)),1,...
                            height(data(i).probe.link))',[],1);
                        dataOut(end+1)=data(i);
                        dataOut(end).probe.link=link;
                        dataOut(end).data=reshape(permute(cfs,[2 3 1]),...
                            size(cfs,2),size(cfs,1)*size(cfs,3));
                    end
                end
                disp(['completed file ' num2str(i) ' of ' num2str(length(data))]);
            end
            data=dataOut;
        end
    end
    
end
