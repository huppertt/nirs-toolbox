function data = simARNoise_shortsep(probe, t, P, sigma, N_files )
    
    if (nargin < 5 || isempty(N_files)), N_files=1; end
    if (nargin < 4 || isempty(sigma)), sigma=.33; end;
    if (nargin < 3  || isempty(P)), P = 10; end
    if (nargin < 2 || isempty(t)), t = (0:1/10:300)'; end
    if (nargin < 1 || isempty(probe)), probe = defaultProbe(); end
    
    if(strcmp(class(probe),'double'))
        probe = defaultProbe(probe);
    end
    
    if (N_files>1)
        data(1:N_files,1) = nirs.core.Data;
        for i = 1:N_files
            data(i) = nirs.testing.simARNoise( probe, t, P, sigma, 1 );
        end
        return;
    end
    
    lambda=unique(probe.link.type)';  %List of wavelengths from the probe
 
    if(0)
        Slab = nirs.core.Image;
        Slab.dim = [5 5 5];
        Slab.origin = [-100 -100 0];
        Slab.description = 'Slab model for FEM/BEM models';
        Slab.vol = ones(41,41,11);  % SLab from -100:5:100, -100:5:100, 0:5:50
        
        % This command will create a nirs.core.Mesh data type
        % This requires iso2mesh to be installed.
        % Cite:
        % Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary
        % and gray-scale images," Proceedings of IEEE International Symposium on Biomedical
        % Imaging 2009, pp. 1142-1145, 2009
        %
        % http://iso2mesh.sourceforge.net
        
        % If the code is not installed, you will get an error and instructions on
        % downloading from the next line.  Iso2Mesh is not distributed with the
        % nirs-toolbox and needs to be seperately downloaded
        
        mesh = Slab.convertFEMmesh();   % This will convert the (binary) Image/Vol to a mesh
        mesh = mesh.reducemesh(.2);  % This will resample the mesh using Iso2Mesh
        
        fwdSlab = nirs.forward.ApproxSlab;
        fwdSlab.Fm=0;
        fwdSlab.probe=probe;
        fwdSlab.prop=nirs.media.tissues.brain(lambda,0.7, 50);
        fwdSlab.mesh=mesh;
    else
        load('SlabModel.mat');
        fwdSlab.probe=probe;
        
    end
    
    
    voxellist=find(fwdSlab.mesh.nodes(:,3)<5);
    nchan=length(voxellist);
    
    % noise mean and spatial covariance
    mu = zeros(nchan,1);
    S = toeplitz( [1 sigma*ones(1,nchan-1)] );
    
    e = mvnrnd( mu, S, length(t) );
    %e = mvnrnd( mu, eye(nchan), length(t) );

    % add temporal covariance
    for i = 1:size(e,2)
        a = randAR( P );
        e(:,i) = filter(1, [1; -a], e(:,i));
    end
    J=fwdSlab.jacobian('spectral');
    e=(J.hbo(:,voxellist)*e')';
    e=e./(ones(length(t),1)*std(e,[],1));
    
    % output
    data = nirs.core.Data();
    data.data   = 100 * exp( - e * 5e-3 );
    data.probe  = probe;
    data.time   = t;
  
end

function a = randAR( P )
    % random Pth order AR coef    
    a = flipud( cumsum( rand(P, 1) ) );
    a = a / sum(a) * 0.99;
end

function probe = defaultProbe(lambda)

if(nargin==0)
    lambda=[690 830];
end
    
    
    srcPos(:,1) = (-80:20:80)';
    srcPos(:,2:3) = 0;
    
    detPos(:,1) = (-70:20:70)';
    detPos(:,2) = 25;
    detPos(:,3) = 0;
    detPos=[detPos; srcPos-ones(size(srcPos,1),1)*[-5 5 0]];
    
    probe = nirs.core.Probe(srcPos,detPos);
    
    link = [1	1	0
            2	1   0	
            2	2   0	
            3	2   0	
            3	3   0	
            4	3   0	
            4	4   0	
            5	4	0
            5	5	0
            6	5	0
            6	6	0
            7	6	0
            7	7	0
            8	7	0
            8	8	0
            9	8	0
            1   9   1
            2   10  1
            3   11  1
            4   12  1
            5   13  1
            6   14  1
            7   15  1
            8   16  1
            9   17  1];
        
    link=[repmat(link,length(lambda),1) reshape(repmat(lambda(:)',size(link,1),1),1,[])'];
    
    link = sortrows(link);
    short = false(size(link,1),1);
    short(end-size(srcPos,1)+1:end)=true;
    probe.link = table(link(:,1), link(:,2), link(:,4),(link(:,3)==1), ...
        'VariableNames', {'source', 'detector', 'type','ShortSeperation'});
        
    
end

function [r,T]=mvnrnd(mu,sigma,cases)


[T,err] = chol(sigma);
r = randn(cases,size(T,1)) * T + ones(cases,1)*mu';

end

