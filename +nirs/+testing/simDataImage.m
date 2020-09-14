function [data, truth, fwdModel, truthchan, nulldata] = simDataImage(fwdModel, noise, stim, beta, basis,SNR )
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
    
    
if(nargin<6)
   SNR=1;
end


    if(nargin<1 || isempty(fwdModel))
        probe = defaultProbe();
    else
        probe=fwdModel.probe;
    end
    

    if nargin < 2 || isempty(noise)
        noise = nirs.testing.simARNoise(probe);
        noise2 = nirs.testing.simARNoise(probe);
    end
    
    if nargin < 3 || isempty(stim)
        stim = nirs.testing.randStimDesign(noise.time, 2, 7, 1);
    end
    
    if(nargin<1 || isempty(fwdModel))
        probe=noise.probe;
        
        Name{1}='FpZ';
        xyz(1,:)=[0 0 0];
        Type{1}='FID-anchor';  % This is an anchor point
        Units{1}='mm';
        
        %Now let's add a few more
        Name{2}='Cz';
        xyz(2,:)=[0 100 0];
        Type{2}='FID-attractor';  % This is an attractor
        Units{2}='mm';
        
        Name{3}='T7';
        xyz(3,:)=[-200 0 0];
        Type{3}='FID-attractor';  % This is an attractor
        Units{3}='mm';
        
        Name{4}='T8';
        xyz(4,:)=[200 0 0];
        Type{4}='FID-attractor';  % This is an attractor
        Units{4}='mm';
        
        % now we need to add these points to the optodes field in the probe.  This
        % is a table, so we need to create a matching format table with the
        % fiducials
        fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        % and concatinate it to the probe
        probe.optodes=[probe.optodes; fid];
        % NOTE- these fiducials are automatically imported from the AtlasViewer
        % format when you use the command nirs.util.sd2probe
        
        probe=nirs.util.registerprobe1020(probe);
        
        lambda=unique(probe.link.type);
        fwdBEM=nirs.registration.Colin27.BEM(lambda);
        
        % Likewise, this will register a mesh onto your probe.  Note- the mesh is
        % the thing that is warped to mathc the head size (not the probe).
        probe=probe.register_mesh2probe(fwdBEM.mesh);
        
        probe.defaultdrawfcn='10-20';
        fwdModel=nirs.forward.ApproxSlab;
        fwdModel.mesh=fwdBEM.mesh(4);
        fwdModel.prop=fwdBEM.prop{4};
        fwdModel.probe=probe;
        
    
    end
    
    
    
    
    if nargin < 4 || isempty(beta)
        beta = {'BA-46_R'};       
    end
    
    if(iscellstr(beta))
        if(length(beta)<stim.count)
            beta=repmat(beta,stim.count,1);
        end
        
        mask=zeros(2*length(fwdModel.mesh(end).nodes),stim.count);
        for i=1:length(beta)
            for j=1:fwdModel.mesh(end).labels.count
                if(ismember(beta{i},fwdModel.mesh(end).labels.values{j}.Label))
                    
                    ii=find(ismember(fwdModel.mesh(end).labels.values{j}.Label,beta{i}));
                    
                    mask(fwdModel.mesh(end).labels.values{j}.VertexIndex{ii},i)=7;
                    mask(end/2+fwdModel.mesh(end).labels.values{j}.VertexIndex{ii},i)=-2;
                    
                end
            end
        end
    else
        mask=beta;  
    end
    
    
    
    if nargin < 5 || isempty(basis)
        % default to canonical basis
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
 
    
    % loop through and add
    data = noise;
    Y    = data.data;
    
    % optical density
    m = mean(Y);
    Y = bsxfun(@plus, -log(Y), log(m));
    
    J = fwdModel.jacobian('spectral');
    
    
    if length(mask) == size(J.hbo,2)
        % oxy; deoxy
        b = [mask; -mask/4];
    else
        b = mask;
    end
    
    
    Xhbo = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbo' );
    Xhbr = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbr' );
    
    Yact = (J.hbo*b(1:end/2,:)*Xhbo'+J.hbr*b(1+end/2:end,:)*Xhbr');
   
    noise.probe=fwdModel.probe;
    noise.data = exp( -bsxfun(@minus, Y, log(m)) );
    noise.stimulus = stim;
    nulldata=noise;
   
  
    Yact=Yact./max(Yact(:))*std(Y(:));
    
    data.probe=fwdModel.probe;
    data.data = exp( -bsxfun(@minus, Y+SNR*Yact', log(m)) );
    data.stimulus = stim;
%     
%     data2 = noise2;
%     Y2    = data2.data;
%     
%     % optical density
%     m2 = mean(Y2);
%     Y2 = bsxfun(@plus, -log(Y2), log(m2));
%     noise2.probe=fwdModel.probe;
%     noise2.data = exp( -bsxfun(@minus, Y2, log(m2)) );
%     noise2.stimulus = stim;
%     nulldata=noise2;
    
    
    truth=(b(:)~=0);
    t=[];
    [un,~,j]=unique(data.probe.link.type);
    for i=1:length(un)
        t(:,i)=max(abs(Yact(j==i,:)),[],2)/max(max(abs(Yact(j==i,:)),[],2));
    end
    %beta=t;
    t=any(t>.25,2);
    truthchan=zeros(size(Y,2),1);
    for i=1:length(un)
        truthchan(j==i)=t;
    end
    
end
    function probe = defaultProbe()

    
    
    srcPos(:,1) = (-80:20:80)';
    srcPos(:,2:3) = 0;
    
    detPos(:,1) = (-70:20:70)';
    detPos(:,2) = 25;
    detPos(:,3) = 0;
    
    probe = nirs.core.Probe(srcPos,detPos);
    
    link = [1	1	690
        2	1	690
        2	2	690
        3	2	690
        3	3	690
        4	3	690
        4	4	690
        5	4	690
        5	5	690
        6	5	690
        6	6	690
        7	6	690
        7	7	690
        8	7	690
        8	8	690
        9	8	690
        1	1	830
        2	1	830
        2	2	830
        3	2	830
        3	3	830
        4	3	830
        4	4	830
        5	4	830
        5	5	830
        6	5	830
        6	6	830
        7	6	830
        7	7	830
        8	7	830
        8	8	830
        9	8	830];
    
    link = sortrows(link);
    
    probe.link = table(link(:,1), link(:,2), link(:,3), ...
        'VariableNames', {'source', 'detector', 'type'});
    
end
    
   