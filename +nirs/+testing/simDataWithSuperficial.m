function [data, truth] = simDataWithSuperficial( noise, stim, beta, channels, basis )
% Simulate data with superficial systemic noise

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
    
    
    %%

noise = noise.sorted();
chan=[noise.probe.link.source noise.probe.link.detector];

[data, truth] = nirs.testing.simData( noise, stim,beta,chan);

% Make the forward model
probe=noise.probe;

braindepth=10;
n=6;

% Compute the optical forward model based on the slab model
minX = min(probe.optodes.X);
maxX = max(probe.optodes.X);
dX = (maxX-minX)/3;
minY = min(probe.optodes.Y);
maxY = max(probe.optodes.Y);
dY = (maxY-minY)/3;

[X,Y,Z]=meshgrid([minX-dX:dX/10:maxX+dX],[minY-dY:dY/10:maxY+dY],[-1:-2:-braindepth]);

mesh=nirs.core.Mesh;
mesh.nodes=[X(:) Y(:) Z(:)];

lambda=unique(probe.link.type);
if(iscell(lambda));
    disp('using 808nm as appromation to hemoglobin fwd-model')
    lambda=808;
    probe.link.type=repmat(lambda,height(probe.link),1);
end;

FwdModel=nirs.forward.ApproxSlab;
FwdModel.mesh=mesh;
FwdModel.prop=nirs.media.tissues.brain(.7,50,lambda);
FwdModel.Fm=0;

FwdModel.probe=probe;
Jacob=FwdModel.jacobian('spectral');
L=[Jacob.hbo Jacob.hbr];
[U,s,V]=nirs.math.mysvd(L);
Lsm = U(:,1:n)*s(1:n,1:n)*V(:,1:n)';
PSF = L*pinv(Lsm);

j=nirs.modules.OpticalDensity;
dOD=j.run(data);
dOD.data=(PSF*dOD.data')';

j=nirs.modules.OpticalDensity2Intensity;
data = j.run(dOD);
data = data.sorted();
 
[data, truth] = nirs.testing.simData(data, stim, beta, channels, basis );


