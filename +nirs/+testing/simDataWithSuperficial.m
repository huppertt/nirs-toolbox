function [data, truth] = simDataWithSuperficial( noise, stim, beta, channels, basis )
% Simulate data with superficial systemic noise

sigma=120;
SNR=2;

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
        channels = sd; %(randi(length(sd),1,round(end/2)),:);
    end
    
    
    %%

noise = noise.sorted();
chan=[noise.probe.link.source noise.probe.link.detector];

[data, truth] = nirs.testing.simData( noise, stim,beta,chan);

% Make the forward model
probe=noise.probe;

braindepth=10;

% Compute the optical forward model based on the slab model
minX = min(probe.optodes.X);
maxX = max(probe.optodes.X);
dX = (maxX-minX)/10;
minY = min(probe.optodes.Y);
maxY = max(probe.optodes.Y);
dY = (maxY-minY);

[X,Y,Z]=meshgrid([minX-dX*2:dX:maxX+dX*2],[minY-dY*2:dY:maxY+dY*2],[-1:-2:-braindepth]);

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

m=size(mesh.nodes,1);
X = repmat(mesh.nodes(:,1),1,m);
Y = repmat(mesh.nodes(:,2),1,m);
Z = repmat(mesh.nodes(:,3),1,m);

X=X-X';
Y=Y-Y';
Z=Z-Z';
Dist = X.^2 + Y.^2 + Z.^2;
Smoother=exp(-Dist/sigma^2);
Smoother=Smoother.*(Smoother>1E-3);
Smoother=sparse(blkdiag(Smoother,Smoother));

FwdModel.probe=probe;
Jacob=FwdModel.jacobian('spectral');
L=[Jacob.hbo Jacob.hbr];
PSF = L*Smoother*L';
PSF=PSF./normest(PSF);

j=nirs.modules.OpticalDensity;
dOD=j.run(data);
dOD.data=(PSF*dOD.data')';

j=nirs.modules.OpticalDensity2Intensity;
data = j.run(dOD);
data = data.sorted();

[X,Y,Z]=meshgrid([minX-dX:dX/3:maxX+dX],[minY-dY:dY/3:maxY+dY],[-15]);
mesh=nirs.core.Mesh;
mesh.nodes=[X(:) Y(:) Z(:)];
FwdModel.mesh=mesh;

[data2, ~, truth] = nirs.testing.simDataImage(FwdModel, data, stim, [], basis );
%[data2, truth] = nirs.testing.simData([], stim, beta, channels, basis );

d=data2.data-data.data;
d=d*SNR/(max(abs(d(:)))/sqrt(mean(var(data.data,[],1))));
data.data=data.data+d;



