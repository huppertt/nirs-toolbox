classdef WaveletTransform < nirs.modules.AbstractModule
%% WaveletTransfrom - Preforms wavelet transfrom on EEG data.
% 
    
    properties
        wname;  % wavelet model name
        scale_smoothing; % add scale smoothing step
        frequencies;  % table of frequency values
        convert2power;
    end
    
    methods

        function obj =WaveletTransform( prevJob )
           obj.name = 'Wavelet Transform using CWT';
           if nargin > 0
               obj.prevJob = prevJob;
           end
           obj.wname='morl';
          
           lower=[0 6 7.5 9 12.5 25]';
           upper=[4 7 12.5 11 30 50]';
           name={'delta','theta','alpha','beta','mu','gamma'}';
           obj.scale_smoothing=false;
           obj.frequencies = table(lower,upper,name);
           obj.convert2power=true;
        end
        
        function data = runThis( obj, data )
            dataOut=data(1);
            dataOut(:)=[];
            for i = 1:length(data)
                FS=data(i).Fs;
                
                scales=[2:128];
                f=scal2frq(scales,obj.wname,1/FS);
                scales = find(f<.5*FS & f>min(obj.frequencies.lower) & ...
                    f<max(obj.frequencies.upper) & f>3/(FS*size(data(i).data,1)));
                f=scal2frq(scales,obj.wname,1/FS);
                flag_SMOOTH = true;
                NTW = fix(FS);
                NSW = fix(length(f)/10);
                
                disp('Progress');
                n=size(data(i).data,2);
                cnt=1;
                str='  0';
                fprintf('%s %',str(end-2:end));
                cfs=zeros(length(f),size(data(i).data,1),size(data(i).data,2));
                 for j=1:size(data(i).data,2)
                    cfs(:,:,j) = cwt(data(i).data(:,j),scales,obj.wname);
                    if(obj.scale_smoothing)
                        cfs(:,:,j)    = smoothCFS(abs(squeeze(cfs(:,:,j))).^2,flag_SMOOTH,NSW,NTW);
                        cfs(:,:,j)    = sqrt(cfs(:,:,j));
                    end
                    str=['   ' num2str(round(100*cnt/n))];
                    fprintf('\b\b\b\b%s %',str(end-2:end));
                    cnt=cnt+1;
                 end
                disp('Computing frequency bands');
                
               
               
                tbl=table;
                ele=data(i).probe.link;
                n=size(data(i).data,2);
                
                for j=1:height(obj.frequencies)
                    lst=find(f>=obj.frequencies.lower(j) & f<=obj.frequencies.upper(j));
                    
                    if( obj.convert2power)
                        d=sqrt(mean(cfs(lst,:,:).^2,1));
                    else
                        d=mean(cfs(lst,:,:),1);
                    end
                    dataOut(end+1)=data(i);
                    dataOut(end).data=squeeze(d);
                    dataOut(end).demographics('freq')=obj.frequencies.name{j};   
                    
                end
                disp(['completed file ' num2str(i) ' of ' num2str(length(data))]);
            end
            data=dataOut;
        end
    end
    
end




%----------------------------------------------------------------------
function CFS = smoothCFS(CFS,flag_SMOOTH,NSW,NTW)

if ~flag_SMOOTH , return; end
if ~isempty(NTW)
    len = NTW;
    F   = ones(1,len)/len;
    CFS = conv2(CFS,F,'same');
end
if ~isempty(NSW)
    len = NSW;
    F   = ones(1,len)/len;    
    CFS = conv2(CFS,F','same');
end
%----------------------------------------------------------------------
end