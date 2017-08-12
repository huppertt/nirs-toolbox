function [J,data,mesh]=jacobian_sp3(mesh,frequency,mesh2)

% [J,data,mesh]=jacobian_sp3(mesh,frequency,mesh2)
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
[f1,f2,g1,g2] = ksi_calc_sp3(mesh);

% Create FEM matrices
if mesh.dimension==2
    [i,j,k1,k3,c,c23,c49,c2_59,B,F1,F2,G1,G2,ib,jb]=gen_matrices_2d_sp3(mesh.nodes(:,1:2),...
        sort(mesh.elements')',...
        mesh.bndvtx,...
        mesh.mua,...
        mesh.mus,...
        mesh.g,...
        f1,...
        f2,...
        g1,...
        g2,...
        mesh.c,...
        omega);
    
    nz_i = nonzeros(i);
    K1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),k1(1:length(nz_i)));
    K3 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),k3(1:length(nz_i)));
    C = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c(1:length(nz_i)));
    C23 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c23(1:length(nz_i)));
    C49 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c49(1:length(nz_i)));
    C2_59 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c2_59(1:length(nz_i)));
    B = sparse(i(1:length(nz_i)),j(1:length(nz_i)),B(1:length(nz_i)));
    nz_i = nonzeros(ib);
    F1 = sparse(ib(1:length(nz_i)),jb(1:length(nz_i)),F1(1:length(nz_i)));
    F2 = sparse(ib(1:length(nz_i)),jb(1:length(nz_i)),F2(1:length(nz_i)));
    G1 = sparse(ib(1:length(nz_i)),jb(1:length(nz_i)),G1(1:length(nz_i)));
    G2 = sparse(ib(1:length(nz_i)),jb(1:length(nz_i)),G2(1:length(nz_i)));
    
elseif mesh.dimension==3
    [i,j,k1,k3,c,c23,c49,c2_59,B,F1,F2,G1,G2,ib,jb]=gen_matrices_3d_sp3(mesh.nodes,...
        sort(mesh.elements')',...
        mesh.bndvtx,...
        mesh.mua,...
        mesh.mus,...
        mesh.g,...
        f1,...
        f2,...
        g1,...
        g2,...
        mesh.c,...
        omega);
    
    nz_i = nonzeros(i);
    K1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),k1(1:length(nz_i)));
    K3 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),k3(1:length(nz_i)));
    C = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c(1:length(nz_i)));
    C23 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c23(1:length(nz_i)));
    C49 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c49(1:length(nz_i)));
    C2_59 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),c2_59(1:length(nz_i)));
    B = sparse(i(1:length(nz_i)),j(1:length(nz_i)),B(1:length(nz_i)));
    F1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),F1(1:length(nz_i)));
    F2 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),F2(1:length(nz_i)));
    G1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),G1(1:length(nz_i)));
    G2 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),G2(1:length(nz_i)));
    
end
clear i* j* k* c* f1 f2 g* nz_i omega

% Add complex component to absorption moments due to frequency dependence
C=C+B;
C23=C23+B;
C49=C49+B;
C2_59=C2_59+B;

% Create MASS matrices
M1 = K1 + C + F1 ;
M2 = K3 + C49 + C2_59 + F2;
MASS=[M1 (G1-C23);(C23-G2) -M2];
clear F* K* M1 M2 C*

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
qvec=[qvec;-(2/3)*qvec];
clear junk i nnodes nsource;

% Catch zero frequency (CW) here
if frequency == 0
    MASS = real(MASS);
    qvec = real(qvec);
end

[phi,R]=get_field(MASS,mesh,qvec);

% Extract composite moments of phi
data.phi1=phi(1:nvtx,:);
data.phi2=phi((nvtx+1):(2*nvtx),:);
clear qvec* phi;

% Calculate scalar flux from contributions of Phi1 Phi2
data.phi = data.phi1-(2/3)*data.phi2;

% Now calculate Adjoint source vector
[qvec] = gen_source_adjoint(mesh);
qvec = [qvec; -(2/3).*qvec];

% Catch zero frequency (CW) here
if frequency == 0
    qvec = real(qvec);
end

% Calculate adjoint field
[phi]=get_field(conj(MASS),mesh,conj(qvec));

data.aphi1=phi(1:nvtx,:);
data.aphi2=phi(nvtx+1:end,:);
data.aphi = data.aphi1-((2/3).*data.aphi2);
clear phi* qvec* R* MASS* *sort;

% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);
[data.complex1]=get_boundary_data(mesh,data.phi1);
[data.complex2]=get_boundary_data(mesh,data.phi2);

% Map complex data to amplitude and phase
% Phi 1
data.amplitude1 = abs(data.complex1);
data.phase1 = atan2(imag(data.complex1),real(data.complex1));
data.phase1(find(data.phase1<0)) = data.phase1(find(data.phase1<0)) + (2*pi);
data.phase1 = data.phase1*180/pi;
data.paa1 = [data.amplitude1 data.phase1];

% Phi 2
data.amplitude2 = abs(data.complex2);
data.phase2 = atan2(imag(data.complex2),real(data.complex2));
data.phase2(find(data.phase2<0)) = data.phase2(find(data.phase2<0)) + (2*pi);
data.phase2 = data.phase2*180/pi;
data.paa2 = [data.amplitude2 data.phase2];

% Total
data.amplitude = abs(data.complex);
data.phase = atan2(imag(data.complex),real(data.complex));
data.phase(find(data.phase<0)) = data.phase(find(data.phase<0)) + (2*pi);
data.phase = data.phase*180/pi;
data.paa = [data.amplitude data.phase];

data2 = data;
ind = data.link(:,3)==0;
data2.complex(ind,:)=[];

% Calculate Jacobian
% Note this is for total fluence wrt kappa and mua

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