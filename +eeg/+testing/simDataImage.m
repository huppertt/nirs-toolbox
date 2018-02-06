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
        beta = {'BA-46_R'};       
    end
    
    if(iscellstr(beta))
        if(length(beta)<stim.count)
            beta=repmat(beta,stim.count,1);
        end
        
        mask=zeros(length(fwdModel.mesh(end).nodes),stim.count);
        for i=1:length(beta)
            for j=1:fwdModel.mesh(end).labels.count
                if(ismember(beta{i},fwdModel.mesh(end).labels.values{j}.Label))
                    
                    ii=find(ismember(fwdModel.mesh(end).labels.values{j}.Label,beta{i}));
                    
                    mask(fwdModel.mesh(end).labels.values{j}.VertexIndex{ii},i)=1;
                   
                    
                end
            end
        end
    else
        mask=beta;  
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
    
    Yact = (J.eeg*mask)*X';
   
    Y = Y+Yact';
    
    data.data = Y;
    data.stimulus = stim;
    
    truth=(mask~=0);
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
    
   