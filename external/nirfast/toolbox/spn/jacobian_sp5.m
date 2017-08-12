function [J,data,mesh]=jacobian_sp5(mesh,frequency,mesh2)

% [J,data,mesh]=jacobian_sp5(mesh,frequency,mesh2)
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
[f1,f2,f3,g1,g2,g3,h1,h2,h3] = ksi_calc_sp5(mesh);

% Create FEM matricex
if mesh.dimension==2
[i,j,k1,k3,k5,c,c23,c49,c64,c8,c16,c2_59,c2_16,c2_49,c4_9,B,F1,F2,F3,G1,G2,G3,H1,H2,H3,ib,jb]=gen_matrices_2d_sp5(mesh.nodes(:,1:2),...
    sort(mesh.elements')',...
    mesh.bndvtx,...
    mesh.mua,...
    mesh.mus,...
    mesh.g,...
    f1,...
    f2,...
    f3,...
    g1,...
    g2,...
    g3,...
    h1,...
    h2,...
    h3,...
    mesh.c,...
    omega);

nz_i=length(nonzeros(i));
K1 = sparse(i(1:nz_i),j(1:nz_i),k1(1:nz_i));
K3 = sparse(i(1:nz_i),j(1:nz_i),k3(1:nz_i));
K5 = sparse(i(1:nz_i),j(1:nz_i),k5(1:nz_i));
C  = sparse(i(1:nz_i),j(1:nz_i),c(1:nz_i));
C23= sparse(i(1:nz_i),j(1:nz_i),c23(1:nz_i));
C49= sparse(i(1:nz_i),j(1:nz_i),c49(1:nz_i));
C64= sparse(i(1:nz_i),j(1:nz_i),c64(1:nz_i));
C8 = sparse(i(1:nz_i),j(1:nz_i),c8(1:nz_i));
C16= sparse(i(1:nz_i),j(1:nz_i),c16(1:nz_i));
C2_59 = sparse(i(1:nz_i),j(1:nz_i),c2_59(1:nz_i));
C2_16 = sparse(i(1:nz_i),j(1:nz_i),c2_16(1:nz_i));
C2_49 = sparse(i(1:nz_i),j(1:nz_i),c2_49(1:nz_i));
C4_9 = sparse(i(1:nz_i),j(1:nz_i),c4_9(1:nz_i));
B = sparse(i(1:nz_i),j(1:nz_i),B(1:nz_i));
nz_i = length(nonzeros(ib));
F1 = sparse(ib(1:nz_i),jb(1:nz_i),F1(1:nz_i));
F2 = sparse(ib(1:nz_i),jb(1:nz_i),F2(1:nz_i));
F3 = sparse(ib(1:nz_i),jb(1:nz_i),F3(1:nz_i));
G1 = sparse(ib(1:nz_i),jb(1:nz_i),G1(1:nz_i));
G2 = sparse(ib(1:nz_i),jb(1:nz_i),G2(1:nz_i));
G3 = sparse(ib(1:nz_i),jb(1:nz_i),G3(1:nz_i));
H1 = sparse(ib(1:nz_i),jb(1:nz_i),H1(1:nz_i));
H2 = sparse(ib(1:nz_i),jb(1:nz_i),H2(1:nz_i));
H3 = sparse(ib(1:nz_i),jb(1:nz_i),H3(1:nz_i));

elseif mesh.dimension==3
[i,j,k1,k3,k5,c,c23,c49,c64,c8,c16,c2_59,c2_16,c2_49,c4_9,B,F1,F2,F3,G1,G2,G3,H1,H2,H3,ib,jb]=gen_matrices_3d_sp5(mesh.nodes,...
    sort(mesh.elements')',...
    mesh.bndvtx,...
    mesh.mua,...
    mesh.mus,...
    mesh.g,...
    f1,...
    f2,...
    f3,...
    g1,...
    g2,...
    g3,...
    h1,...
    h2,...
    h3,...
    mesh.c,...
    omega);

nz_i=length(nonzeros(i));
K1 = sparse(i(1:nz_i),j(1:nz_i),k1(1:nz_i));
K3 = sparse(i(1:nz_i),j(1:nz_i),k3(1:nz_i));
K5 = sparse(i(1:nz_i),j(1:nz_i),k5(1:nz_i));
C  = sparse(i(1:nz_i),j(1:nz_i),c(1:nz_i));
C23= sparse(i(1:nz_i),j(1:nz_i),c23(1:nz_i));
C49= sparse(i(1:nz_i),j(1:nz_i),c49(1:nz_i));
C64= sparse(i(1:nz_i),j(1:nz_i),c64(1:nz_i));
C8 = sparse(i(1:nz_i),j(1:nz_i),c8(1:nz_i));
C16= sparse(i(1:nz_i),j(1:nz_i),c16(1:nz_i));
C2_59 = sparse(i(1:nz_i),j(1:nz_i),c2_59(1:nz_i));
C2_16 = sparse(i(1:nz_i),j(1:nz_i),c2_16(1:nz_i));
C2_49 = sparse(i(1:nz_i),j(1:nz_i),c2_49(1:nz_i));
C4_9 = sparse(i(1:nz_i),j(1:nz_i),c4_9(1:nz_i));
B = sparse(i(1:nz_i),j(1:nz_i),B(1:nz_i));
F1 = sparse(i(1:nz_i),j(1:nz_i),F1(1:nz_i));
F2 = sparse(i(1:nz_i),j(1:nz_i),F2(1:nz_i));
F3 = sparse(i(1:nz_i),j(1:nz_i),F3(1:nz_i));
G1 = sparse(i(1:nz_i),j(1:nz_i),G1(1:nz_i));
G2 = sparse(i(1:nz_i),j(1:nz_i),G2(1:nz_i));
G3 = sparse(i(1:nz_i),j(1:nz_i),G3(1:nz_i));
H1 = sparse(i(1:nz_i),j(1:nz_i),H1(1:nz_i));
H2 = sparse(i(1:nz_i),j(1:nz_i),H2(1:nz_i));
H3 = sparse(i(1:nz_i),j(1:nz_i),H3(1:nz_i));

end

clear junk i j k1 k3 k5 c23 c49 c64 c2_59 c2_16 c2_8 c4_9 omega ib jb nz_i

% Add complex component to absorption moments due to 
% frequency dependence
C=C+B;
C23=C23+B;
C49=C49+B;
C64=C64+B;
C8=C8+B;
C16=C16+B;
C2_59=C2_59+B;
C2_16=C2_16+B;
C2_49=C2_49+B;
C4_9=C4_9+B;
Cs=C16+C2_49;


M1=K1+C+F1;
M2=K3+C49+C2_59+F2;
M3=K5+C64+C2_16+C4_9+F3;

% Calculate the RHS (the source vectors. For simplicity, we are
% just going to use a Gaussian Source, The width of the Gaussian is
% changeable (last argument). The source is assumed to have a
% complex amplitude of complex(cos(0.15),sin(0.15));

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

MASS=[M1 -(C23-G1) (C8-H1);-(C23-G2) M2 -(C16+C2_49-H2);
     (C8-G3) -(C16+C2_49-H3) M3];
qvec = [qvec; -(2/3).*qvec; (8/15).*qvec];

% Catch zero frequency (CW) here
if frequency == 0
    MASS = real(MASS);
    qvec = real(qvec);
end

[phi,R]=get_field(MASS,mesh,qvec);

% Re-order elements
data.phi1=phi(1:nvtx,:);
data.phi2=phi((nvtx+1):(2*nvtx),:);
data.phi3=phi(((2*nvtx)+1):3*nvtx,:);

% Calculate overall scalar flux due to contributions from Phi 1,2,3
data.phi = data.phi1 - (2/3).*data.phi2 + (8/15).*data.phi3;

% ==========================================

% Generate adjoint source vector
[qvec] = gen_source_adjoint(mesh);
qvec = [qvec;-(2/3)*qvec;(8/15)*qvec];

% Catch zero frequency (CW) here
if frequency == 0
  qvec = real(qvec);
end

% Calculate adjoint field
[phi]=get_field(conj(MASS),mesh,conj(qvec));
data.aphi1=phi(1:nvtx,:);
data.aphi2=phi((nvtx+1):(2*nvtx),:);
data.aphi3=phi(((2*nvtx)+1):3*nvtx,:);
data.aphi = data.aphi1 - (2/3).*data.aphi2 + (8/15).*data.aphi3;

clear phi* qvec* R* MASS*;
% =======================================

% Calculate boundary data for each composite moment and total
[data.complex1]=get_boundary_data(mesh,data.phi1);
[data.complex2]=get_boundary_data(mesh,data.phi2);
[data.complex3]=get_boundary_data(mesh,data.phi3);
[data.complex]=get_boundary_data(mesh,data.phi);

% Map complex data to amplitude and phase
% Phi 1
data.amplitude1 = abs(data.complex1);
data.phase1 = atan2(imag(data.complex1),...
		   real(data.complex1));
data.phase1(find(data.phase1<0)) = data.phase1(find(data.phase1<0)) + (2*pi);
data.phase1 = data.phase1*180/pi;

data.paa1 = [data.amplitude1 data.phase1];

% Phi 2
data.amplitude2 = abs(data.complex2);
data.phase2 = atan2(imag(data.complex2),...
		   real(data.complex2));
data.phase2(find(data.phase2<0)) = data.phase2(find(data.phase2<0)) + (2*pi);
data.phase2 = data.phase2*180/pi;

data.paa2 = [data.amplitude2 data.phase2];

% Phi 3
data.amplitude3 = abs(data.complex3);
data.phase3 = atan2(imag(data.complex3),...
		   real(data.complex3));
data.phase3(find(data.phase3<0)) = data.phase3(find(data.phase3<0)) + (2*pi);
data.phase3 = data.phase3*180/pi;

data.paa3 = [data.amplitude3 data.phase3];

% Total
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
