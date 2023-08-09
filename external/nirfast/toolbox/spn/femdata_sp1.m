function [data,mesh]=femdata_sp1(mesh,frequency)

% [data,mesh]=femdata_sp1(mesh,frequency)
% Calculates data (phase and amplitude) for a given
% standard mesh at a given frequency (MHz) based on SP1 approximation.
% outputs phase and amplitude in structure data
% and mesh information in mesh


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
nvtx=length(mesh.nodes);



% Calculate boundary condition coefficients
[f1] = ksi_calc_sp1(mesh);

% Create FEM matrices
if mesh.dimension == 2
  [i,j,k1,c,B,F1,ib,jb]=gen_matrices_2d_sp1(mesh.nodes(:,1:2),...
    sort(mesh.elements')',...
    mesh.bndvtx,...
    mesh.mua,...
    mesh.mus,...
    mesh.g,...
    f1,...
    mesh.c,...
    omega);

nz_i = nonzeros(i);
K1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),k1(1:length(nz_i)));
C = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c(1:length(nz_i)));
B = sparse(i(1:length(nz_i)),j(1:length(nz_i)),B(1:length(nz_i)));
nz_i = nonzeros(ib);
F1 = sparse(ib(1:length(nz_i)),jb(1:length(nz_i)),F1(1:length(nz_i)));

elseif mesh.dimension == 3
     [i,j,k1,c,B,F1,ib,jb]=gen_matrices_3d_sp1(mesh.nodes,...
    sort(mesh.elements')',...
   mesh.bndvtx,...
    mesh.mua,...
    mesh.mus,...
    mesh.g,...
    f1,...
    mesh.c,...
    omega);
nz_i = nonzeros(i);
K1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),k1(1:length(nz_i)));
C = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c(1:length(nz_i)));
B = sparse(i(1:length(nz_i)),j(1:length(nz_i)),B(1:length(nz_i)));
F1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),F1(1:length(nz_i)));

end

clear i j k1 c f1 g1 ib jb nz_i omega 

MASS = K1 + C + B +F1;

clear C* B* F* G*

% Now calculate source vector
% NOTE last term in mex file 'gen_source' is the source FWHM
%
source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);
[nsource,junk]=size(source);
qvec = spalloc(nnodes,nsource,nsource*100);
if mesh.dimension == 2
  for i = 1 : nsource
    s_ind = mesh.source.num == source(i);
    if mesh.source.fwhm(s_ind) == 0
        qvec(:,i) = gen_source_point(mesh,mesh.source.coord(s_ind,1:2));
    else
      qvec(:,i) = gen_source(mesh.nodes(:,1:2),...
			   sort(mesh.elements')',...
			   mesh.dimension,...
			   mesh.source.coord(s_ind,1:2),...
			   mesh.source.fwhm(s_ind));
    end
  end
elseif mesh.dimension == 3
  for i = 1 : nsource
    s_ind = mesh.source.num == source(i);
    if mesh.source.fwhm(s_ind) == 0
        qvec(:,i) = gen_source_point(mesh,mesh.source.coord(s_ind,1:3));
    else
    qvec(:,i) = gen_source(mesh.nodes,...
			   sort(mesh.elements')',...
			   mesh.dimension,...
			   mesh.source.coord(s_ind,:),...
			   mesh.source.fwhm(s_ind));
    end
  end
end

clear junk i nnodes nsource;

% Catch zero frequency (CW) here
if frequency == 0
  MASS = real(MASS);
  qvec = real(qvec);
end

% ======================================
% Optimise MASS matrix
 
% [MASS_opt,Q_opt,invsort]=optimise(MASS,qvec);
MASS_opt = MASS;
Q_opt = qvec;


% =======================================
% 
 phi_all=get_field(MASS_opt,mesh,Q_opt);

% Re-order elements
% phi_all=phi_all(invsort,:);

% Extract composite moments of phi
data.phi=phi_all(1:nvtx,:);


clear qvec* M* C;


% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);
data.link = mesh.link;

% Map complex data to amplitude and phase
data.amplitude = abs(data.complex);

data.phase = atan2(imag(data.complex),...
		   real(data.complex));
data.phase(find(data.phase<0)) = data.phase(find(data.phase<0)) + (2*pi);
data.phase = data.phase*180/pi;

data.paa = [data.amplitude data.phase];
