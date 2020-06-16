function [J,data,mesh]=jacobian_sp1(mesh,frequency,mesh2)

% [J,data,mesh]=jacobian_sp1(mesh,frequency,mesh2)
%
% Calculates the jacobian for the total momentof phase and amplitude for
% scatter and absoption
% See any of Dartmouth Publications regarding the
% structure. Also calculates data (phase and amplitude) for a given
% problem (fn) at a given frequency (MHz).
%
% mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% mesh2 is optional mesh basis

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

% Calculate boundary coefficients
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
% NOTE last term in mex file 'qvec' is the source FWHM

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
 
 [phi,R]=get_field(MASS,mesh,qvec);

% Extract composite moments of phi
  data.phi=phi(1:nvtx,:);

%=======================================
% Now calculate Adjoint source vector
[qvec] = gen_source_adjoint(mesh);

% Catch zero frequency (CW) here
if frequency == 0
  qvec = real(qvec);
end

 phi=get_field(conj(MASS),mesh,conj(qvec));

% Extract composite moments of phi
data.aphi=phi(1:nvtx,:);

clear qvec* Q* R MASS* phi_all;

% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);

% Map complex data to amplitude and phase
data.amplitude = abs(data.complex);

data.phase = atan2(imag(data.complex),...
		   real(data.complex));
data.phase(find(data.phase<0)) = data.phase(find(data.phase<0)) + (2*pi);
data.phase = data.phase*180/pi;

data.paa = [data.amplitude data.phase];

data2 = data;
ind = data.link(:,3)==0;
data2.complex(ind,:)=[];

% Calculate Jacobian
% Catch zero frequency (CW) here
if nargin == 3 % use second mesh basis for jacobian
    data3 = interpolatef2r(mesh,mesh2,data2);
    data3.complex = data2.complex;
    % Calculate Jacobian
    % Catch zero frequency (CW) here
    if frequency == 0
        [J] = build_jacobian_cw(mesh2,data3);
    else
        [J] = build_jacobian(mesh2,data3);
    end
elseif nargin == 2
    % Calculate Jacobian
    % Catch zero frequency (CW) here
    if frequency == 0
        [J] = build_jacobian_cw(mesh,data2);
    else
        [J] = build_jacobian(mesh,data2);
    end
end

% We do not want jacobian for kappa, but for mus
% which is d_phi/d_mus = d_phi/d_kappa * d_kappa/d_musp * d_musp/d_mus
if nargin == 3 % use second mesh basis for jacobian
    musp = mesh2.mus.*(1-mesh2.g);
    kappa = 1./(3.*(mesh2.mua+musp));
    d_musp_d_mus = (1-mesh2.g);
elseif nargin == 2
    musp = mesh.mus.*(1-mesh.g);
    kappa = 1./(3.*(mesh.mua+musp));
    d_musp_d_mus = (1-mesh.g);
end
for i = 1 : size(J.complex,1)
    if frequency == 0
        %J.complex(i,:) = J.complex(i,:).*(-3.*kappa.*kappa)'.*d_musp_d_mus';
    else
        J.complex(i,1:end/2) = J.complex(i,1:end/2).*(-3.*kappa.*kappa)'.*d_musp_d_mus';
    end
end
for i = 1 : size(J.complete,1)
    if frequency == 0
        %J.complete(i,:) = J.complete(i,:).*(-3.*kappa.*kappa)'.*d_musp_d_mus';
    else
        J.complete(i,1:end/2) = J.complete(i,1:end/2).*(-3.*kappa.*kappa)'.*d_musp_d_mus';
    end
end

clear musp kappa d_mus*



function data2 = interpolatef2r(fwd_mesh,recon_mesh,data)
% This function interpolates fwd_mesh data into recon_mesh
% Used to calculate the Jacobian on second mesh

for i = 1 : length(recon_mesh.nodes)
    if fwd_mesh.fine2coarse(i,1) ~= 0
        data2.phi(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
            data.phi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
        data2.aphi(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
            data.aphi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
    elseif fwd_mesh.fine2coarse(i,1) == 0
        dist = distance(mesh.nodes,...
            mesh.bndvtx,...
            pixel.nodes(i,:));
        mindist = find(dist==min(dist));
        mindist = mindist(1);
        data2.phi(i,:) = data.phi(mindist,:);
        data2.aphi(i,:) = data.phi(mindist,:);
    end
end

