function [data, truth, truthchan, beta] = simDataImage(fwdModel, noise, stim, beta, basis )
%SIMDATA Simulates NIRS data by adding a task to baseline noise.
% 
% Args:
%     fwdModel - a nirs.forward.() model class
%     noise  -  a raw data file
%     stim   -  dictionary containing stimulus objects using the stim name as key
%     beta   -  vector of magnitude of response for each stim condition (nVox x nCond)
%     basis  -  dictionary containing basis objects using stim condition as key
%     
% Example:
%     noise = nirs.testing.simARNoise();
%     stim  = nirs.testing.randStimDesign(noise.time, 2, 7, 3);
%     beta  = [3 2 1]';
%     
%     sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
%     channels = sd(1:round(end/2),:);
%
%     [data, truth] = simData( noise, stim, beta, channels )
    
    if nargin < 2 || isempty(noise)
        noise = nirs.testing.simARNoise();
    end
    
    if nargin < 3 || isempty(stim)
        stim = nirs.testing.randStimDesign(noise.time, 2, 7, 1);
    end
    
    if nargin < 4 || isempty(beta)
        nVox=size(fwdModel.mesh.nodes,1);
        beta=zeros(nVox,length(stim.keys));
        
        pos=(fwdModel.probe.srcPos(fwdModel.probe.link.source,:)+...
            fwdModel.probe.detPos(fwdModel.probe.link.detector,:))/2;
        
        for idx=1:length(stim.keys)
            p=pos(randi(size(pos,1),1,1),:);
            d=sqrt(sum((fwdModel.mesh.nodes-ones(nVox,1)*p).^2,2));
            lst=find(d<20 & d>5 & abs(fwdModel.mesh.nodes(:,3))>8);
            beta(lst,idx)=7;
        end
    end
    
    
    
    if nargin < 5 || isempty(basis)
        % default to canonical basis
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
 
    
    % loop through and add
    data = noise.sorted();
    
    link = data.probe.link;
    Y    = data.data;
    truth = zeros(size(Y,2), 1);
    
    % optical density
    m = mean(Y);
    Y = bsxfun(@plus, -log(Y), log(m));
    
    J = fwdModel.jacobian('spectral');
    
    
    if length(beta) == size(J.hbo,2)
        % oxy; deoxy
        b = [beta; -beta/4];
    else
        b = beta;
    end
    
    
    Xhbo = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbo' );
    Xhbr = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbr' );
    
    Yact = (J.hbo*b(1:end/2,:)*Xhbo'+J.hbr*b(1+end/2:end,:)*Xhbr');
   
    Y = Y+Yact';
    
    Y = exp( -bsxfun(@minus, Y, log(m)) );
    
    data.data = Y;
    data.stimulus = stim;
    
    truth=(b(:)~=0);
    t=[];
    [un,~,j]=unique(link.type);
    for i=1:length(un)
        t(:,i)=max(abs(Yact(j==i,:)),[],2)/max(max(abs(Yact(j==i,:)),[],2));
    end
    %beta=t;
    t=any(t>.25,2);
    truthchan=zeros(size(Y,2),1);
    for i=1:length(un)
        truthchan(j==i)=t;
    end
    
   