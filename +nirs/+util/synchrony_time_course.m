function dataOut = synchrony_time_course(data,method,freqband)
% This function computes "instantanious" synchrony of data and returns as a time
% course nirs.core.Dara type
%
% Inputs:
% data - the nirs.core.Data input
% method - One of the following
%     Pearson (default)
%     AR-Pearson
%     Wavelet
%
% freq-band(if wavelet)
%   [lower upper] in Hz
%
% remove_nonsig (string; e.g. "p<0.05")
%   If provided, this is the p< or q< cutoff for keeping the data


if(nargin<2 || isempty(method))
    method='Pearson';
end
if(nargin<3 || isempty(freqband))
    freqband = [0.01 .1];
end
if(nargin<4 || isempty(remove_nonsig))
    remove_nonsig=[];
end


data.data=data.data-ones(size(data.data,1),1)*mean(data.data,1);
data.data=data.data-ones(size(data.data,1),1)*mean(data.data,1);

connections=[];

ynew=[]; %zeros(size(data.data,1),size(data.data,2)*(size(data.data,2)-1)/2);
if(strcmp(lower(method),'pearson'))
    cnt=1;
    for i=1:size(data.data,2)
        for j=i+1:size(data.data,2)
            if(ismember(data.probe.link(i,:).type,data.probe.link(j,:).type))
                connections.start(cnt,1)=i;
                connections.end(cnt,1)=j;
                connections.type(cnt,1)=data.probe.link.type(i);

                ynew(:,cnt)=data.data(:,i).*data.data(:,j);
                cnt=cnt+1;
            end
        end
    end
    d=nirs.math.innovations(data.data,round(data.Fs*4));
    [~,p]=nirs.math.robust_corrcoef(d);

elseif(strcmp(lower(method),'ar-pearson'))
    d=nirs.math.innovations(data.data,round(data.Fs*4));
    cnt=1;
    for i=1:size(data.data,2)
        for j=i+1:size(data.data,2)
            if(ismember(data.probe.link(i,:).type,data.probe.link(j,:).type))
                connections.start(cnt,1)=i;
                connections.end(cnt,1)=j;
                if(iscell(data.probe.type))
                    connections.type{cnt}=data.probe.type{cnt};
                else
                    connections.type(cnt)=data.probe.type(cnt);
                
                end

                ynew(:,cnt)=d(:,i).*d(:,j);
                cnt=cnt+1;
            end
        end
    end
    [~,p]=nirs.math.robust_corrcoef(d);

elseif(strcmp(lower(method),'wavelet'))
    cnt=1;
    for ch=1:size(data.data,2)
        [dw(:,:,ch),f]=nirs.math.robust_cwt(data.data(:,ch));
    end
    lstfq = find(f>=freqband(1) & f<=freqband(2));
    for i=1:size(data.data,2)
        for j=i+1:size(data.data,2)
            if(ismember(data.probe.link(i,:).type,data.probe.link(j,:).type))
                connections.start(cnt,1)=i;
                connections.end(cnt,1)=j;
                connections.type(cnt,1)=data.probe.link.type(i);
                yy=dw(:,:,i).*conj(dw(:,:,j));
                ynew(:,cnt)=median(yy(lstfq,:),1);

                cnt=cnt+1;
            end
        end

    end
    dr=nirs.math.innovations(real(squeeze(median(dw(lstfq,:,:),1))),round(data.Fs*4));
    di=nirs.math.innovations(imag(squeeze(median(dw(lstfq,:,:),1))),round(data.Fs*4));
    [~,p]=nirs.math.robust_corrcoef(dr+di*sqrt(-1));

else
    error('unknown method provide')
end

% 
% 
% [id,jd]=ind2sub([height(data.probe.link),height(data.probe.link)],lstsif);
% connections=[];
% connections.start=id(:);
% connections.end=jd(:);
% connections.type=data.probe.link(id,:).type;
 connections=struct2table(connections);
% lst=find(strcmp(data.probe.link(jd,:).type,data.probe.link(id,:).type));
% connections=connections(lst,:);


dataOut = data;
dataOut.probe=nirs.core.ProbeConnections(data.probe);

if(strcmp(lower(method),'wavelet'))
    dataOut.data = [abs(ynew) pi*unwrap(angle(ynew)/pi,[],1)];
    cm=connections;
    cp=connections;
    for i=1:length(cm.type); cm.type{i}=[cm.type{i} '_magnitude']; end;
    for i=1:length(cp.type); cp.type{i}=[cp.type{i} '_phase']; end;
    dataOut.probe.connections=[cm;cp];
else
    dataOut.data = ynew;
    dataOut.probe.connections=connections;
end




