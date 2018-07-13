function [data, truth] = simData_registered( noise, stim, beta, channels, basis )
%SIMDATA Simulates NIRS data by adding a task to baseline noise.
% 
% Args:
%     noise  -  a raw data file
%     stim   -  dictionary containing stimulus objects using the stim name as key
%     beta   -  vector of magnitude of response for each stim conditon
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
    
    if nargin < 1 || isempty(noise)
        noise = nirs.testing.simARNoise();
    end
    
    if nargin < 2 || isempty(stim)
        stim = nirs.testing.randStimDesign(noise.time, 2, 7, 1);
    end
    
    if nargin < 3 || isempty(beta)
        beta = 7*ones( length(stim.keys), 1 );
    end
    
    if length(beta) == length(stim.keys)
        % oxy; deoxy
        b = [beta; -beta/2];
    else
        b = beta;
    end
    
    if nargin < 5 || isempty(basis)
        % default to canonical basis
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
    if nargin < 4 || isempty(channels)
        % default to first half of channels
        sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
        channels = sd(1:round(end/2),:);
    end
    
    [data, truth] = nirs.testing.simData( noise, stim, beta, channels, basis );
    
    
    probe=data.probe;
    
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
   
    probe1020=nirs.util.registerprobe1020(probe);
    
    lambda=unique(probe1020.link.type);
    fwdBEM=nirs.registration.Colin27.BEM(lambda);
    
     
    % Likewise, this will register a mesh onto your probe.  Note- the mesh is
    % the thing that is warped to mathc the head size (not the probe).
    probe1020=probe1020.regsister_mesh2probe(fwdBEM.mesh);
    
    probe1020.defaultdrawfcn='10-20';
    data.probe=probe1020;
    
end

