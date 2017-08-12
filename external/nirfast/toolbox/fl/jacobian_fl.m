function [J,data]=jacobian_fl(mesh,frequency,datax)

% [J,data]=jacobian_fl_new_new(fn,frequency,mesh,data)
%
% Calculats jacobian for fluorescence yield.
% Ensure data input is excitation field data!!
%
% mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% datax is the excitation field data (variable)


% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

% modulation frequency
omega = 2*pi*frequency*1e6;

% Create Emission FEM Matrix
if mesh.dimension == 2
    [i,j,s] = gen_matrices_2d(mesh.nodes(:,1:2),...
        sort(mesh.elements')', ...
        mesh.bndvtx,...
        mesh.muam,...
        mesh.kappam,...
        mesh.ksi,...
        mesh.c,...
        omega);
    
elseif mesh.dimension ==3
    [i,j,s] = gen_matrices_3d(mesh.nodes,...
        sort(mesh.elements')', ...
        mesh.bndvtx,...
        mesh.muam,...
        mesh.kappam,...
        mesh.ksi,...
        mesh.c,...
        omega);
end
junk = length(find(i==0));
MASS_m = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
clear junk i j s;
  
% If the fn.ident exists, then we must modify the FEM matrices to
% account for refractive index mismatch within internal boundaries
if isfield(mesh,'ident') == 1
    disp('Modifying for refractive index')
    M = bound_int(MASS,mesh);
    MASS = M;
    clear M
end

% Calculate the RHS (the source vectors) for the Emission.
source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
nsource = length(source);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qvec = zeros(nnodes,nsource);
% Simplify the RHS of emission equation
beta = mesh.gamma.*(1-(sqrt(-1).*omega.*mesh.tau));
% get rid of any zeros!
if frequency == 0
    beta(beta==0) = 1e-20;
else
    beta(beta==0) = complex(1e-20,1e-20);
end

if mesh.dimension == 2
  for i = 1 : nsource
    val = beta.*datax.phi(:,i);
    qvec(:,i) = gen_source_fl(mesh.nodes(:,1:2),...
			      sort(mesh.elements')',...
			      mesh.dimension,...
			      val);
  end
elseif mesh.dimension == 3
  for i = 1 : nsource
    val = beta.*datax.phi(:,i);
    qvec(:,i) = gen_source_fl(mesh.nodes,...
			      sort(mesh.elements')',...
			      mesh.dimension,...
			      val);
  end
end

clear junk i nnodes nsource val beta;

% Calculate EMISSION field for all sources
[data.phim,mesh.R]=get_field(MASS_m,mesh,qvec);
clear qvec;

% Calculate Adjoint source vector
[qvec] = gen_source_adjoint(mesh);

% Calculate adjoint field for all detectors
[data.aphim]=get_field(conj(MASS_m),mesh,conj(qvec));
clear qvec R MASS_m;

% Calculate boundary data
[data.complexm]=get_boundary_data(mesh,data.phim);
data.link = mesh.link;

% Map complex data to amplitude and phase
data.amplitudem = abs(data.complexm);

data.phasem = atan2(imag(data.complexm),...
		   real(data.complexm));
data.phasem(data.phasem<0) = data.phasem(data.phasem<0) + (2*pi);
data.phasem = data.phasem*180/pi;

data.paam = [data.amplitudem data.phasem];
data.phix = datax.phi;

% Build the Emission jacobian
data2 = data;
ind = data.link(:,3) == 0;
data2.complexm(ind,:)=[];
    
if omega == 0
    [J] = build_jacobian_cw_fl(mesh,data2,omega);
else
    [J] = build_jacobian_fl(mesh,data2,omega);
end
