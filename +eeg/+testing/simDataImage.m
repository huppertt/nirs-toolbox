function [data, truth, truthchan, beta] = simDataImage(fwdModel, noise, stim, beta, basis )
%SIMDATA Simulates EEG data by adding a task to baseline noise.
% 
% Args:
%     fwdModel - a eeg.forward.() model class
%     noise  -  a raw data file
%     stim   -  dictionary containing stimulus objects using the stim name as key
%     beta   -  vector of magnitude of response for each stim condition (nVox x nCond)
%     basis  -  dictionary containing basis objects using stim condition as key
%     

    
    if nargin < 2 || isempty(noise)
        noise = eeg.testing.simARNoise();
    end
    
    if nargin < 3 || isempty(stim)
        stim = nirs.testing.randStimDesign(noise.time, .1, 3, 1);
    end
    
    if nargin < 4 || isempty(beta)
        nVox=size(fwdModel.mesh.mesh(end).nodes,1);
        beta=zeros(nVox,length(stim.keys));
        
        pos=[fwdModel.probe.electrodes.X fwdModel.probe.electrodes.Y fwdModel.probe.electrodes.Z];
        
        for idx=1:length(stim.keys)
            p=pos(randi(size(pos,1),1,1),:);
            d=sqrt(sum((fwdModel.mesh.mesh(end).nodes-ones(nVox,1)*p).^2,2));
            lst=find(d<30 & d>5);
            beta(lst,idx)=7;
        end
    end
    
    
    
    if nargin < 5 || isempty(basis)
        % default to canonical basis
        basis = Dictionary({'default'}, {eeg.design.basis.ERP()});
    end
    
 
    
    % loop through and add
    data = noise.sorted();
    
    link = data.probe.link;
    Y    = data.data;
    truth = zeros(size(Y,2), 1);
    
     
    J = fwdModel.jacobian;
        
    
    X = nirs.design.createDesignMatrix( stim, data.time, basis);
    
    Yact = J.eeg*beta*X';
   
    Y = Y+Yact';
    
    data.data = Y;
    data.stimulus = stim;
    
    truth=(beta~=0);
    t=[];
    [un,~,j]=unique(link.type);
    for i=1:length(un)
        t(:,i)=max(abs(Yact(j==i,:)),[],2)/max(max(abs(Yact(j==i,:)),[],2));
    end
    beta=t;
    t=any(t>.25,2);
    truthchan=zeros(size(Y,2),1);
    for i=1:length(un)
        truthchan(j==i)=t;
    end
    
   