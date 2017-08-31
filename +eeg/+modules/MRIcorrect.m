classdef MRIcorrect < nirs.modules.AbstractModule
%% MRI artifact correction for EEG
%

    methods

        function obj = MRIcorrect( prevJob )
           obj.name = 'MRI Artifact correction';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                TRpulse=data(i).stimulus('MRI_Pulse');
                if(isempty(TRpulse))
                    error('data must have MRI times marked in the stimulus field');
                end
                disp(['Processing scan ' num2str(i) ' of ' num2str(length(data))]); 
                
                TR=mean(diff(TRpulse.onset));
                startMR=min(TRpulse.onset);
                endMR=max(TRpulse.onset);
                
                d=data(i).data(find(data(i).time>=startMR & data(i).time<=endMR),:);
                
                [fa,fb]=butter(4,5*2/data(i).Fs,'high');
                
                [U,S,V]=nirs.math.mysvd(filtfilt(fa,fb,d));
                
                k=dsearchn(data(i).time,TRpulse.onset);
                s=zeros(size(data(i).time));
                s(k)=1;
                s=s(find(data(i).time>=startMR & data(i).time<=endMR));
                
                X=convmtx(s,ceil(TR*1.05.*data(i).Fs));
                X=X(1:length(s),:);
                iX=pinv(X);
                
                beta = iX*U;
                beta=beta(50:end-50,:);  % remove a bit of the edges
                yfilt=[]; f={};
                disp('Apply AR filters');
                for j=1:size(beta,2)
                    a=nirs.math.ar_fit(beta(:,j),data(i).Fs*4);
                    f=[1; -a(2:end)];
                    yfilt(:,j) = filter(f, 1, data(i).data*V(:,j));
                    nirs.util.flushstdout(1);
                    fprintf( 'Finished %4i of %4i.\n', j, size(beta,2) )
                end
                
                dd=yfilt*V';
                
                data(i).data=dd;
                
                
            end
        end
    end
    
end

