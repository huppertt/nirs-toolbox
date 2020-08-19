function [data,truth]=simData_connectivity_shortsep
sigma=150; %units-(mm) % spatial smoothing kernel for the skin layer 
pmax=10;  % model order to use for the AR model
t = (0:1/10:300)';

SNR=1;  % ratio of skin to brain signals

probe = defaultProbe;

   lambda=unique(probe.link.type)';  %List of wavelengths from the probe
 
   % This will rerun the forward model, but there is a saved version that
   % will be used
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

    % find all the superficial (<5mm) voxels
    voxellist=find(fwdSlab.mesh.nodes(:,3)<5);
    nchan=length(voxellist);
    
    % noise mean and spatial covariance
    mu = zeros(nchan,1);
    % Noise model based on spatial distance
    S=exp(-squareform(pdist(fwdSlab.mesh.nodes(voxellist,:))).^2/sigma^2);
    
    
    e = mvnrnd( mu, S, length(t) );
    %e = mvnrnd( mu, eye(nchan), length(t) );

    % add temporal covariance
    a = randAR( pmax );
    for i = 1:size(e,2) 
        e(:,i) = filter(1, [1; -a], e(:,i));
    end
    
      
    J=fwdSlab.jacobian('spectral');
    e=(J.hbo(:,voxellist)*e')';
    e=6*e./(ones(length(t),1)*std(e,[],1));
    
    % now the connectivity part
    types=probe.types;
    et=zeros(length(t),height(probe.link));
    T=zeros(height(probe.link));
    for id=1:length(types)
        lst=find(~probe.link.ShortSeperation & ...
            probe.link.type==types(id));
        if(id==1)
            fract=.1;
            truth=(rand(length(lst))>(1-fract/sqrt(length(lst))));
            truth=truth+eye(length(lst));
            
            truth=truth*truth';
            truth=min(truth,1);
            truth=max(truth,-1);
        end
        et(:,lst) = mvnrnd( zeros(length(lst),1),truth, length(t) );
        T(lst,lst)=truth;
        a = randAR( pmax );
        for i = 1:length(lst)
            et(:,lst(i)) = filter(1, [1; -a], et(:,lst(i)));
        end
    end
lst=find(~probe.link.ShortSeperation);
et(:,lst)=SNR*6*et(:,lst)./(ones(length(t),1)*std(et(:,lst),[],1));

y=(e+et)+200;
y=y./(ones(size(y,1),1)*mean(y,1));

data = nirs.core.Data();
data.data   = 100 * exp( -y);
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


[T,err] = cholcov(sigma);
r = randn(cases,size(T,1)) * T + ones(cases,1)*mu';

end
