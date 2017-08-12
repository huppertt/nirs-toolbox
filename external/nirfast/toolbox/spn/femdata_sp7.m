function [data,mesh]=femdata_sp7(mesh,frequency)

% [data,mesh]=femdata_sp7(mesh,frequency)
% Calculates data (phase and amplitude) for a given
% standard mesh at a given frequency (MHz).
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

% Calculate boundary coefficients
[f1,f2,f3,f4,g1,g2,g3,g4,h1,h2,h3,h4,i1,i2,i3,i4] = ksi_calc_sp7(mesh);

% Create FEM matrices
if mesh.dimension==2

  [i,j,k1,k3,k5,k7,c,c23,c49,c64,c8,c16,c32,c128,c256,c2_59,c2_16,c2_49,...
    c2_8,c2_32,c2_64,c4_9,c4_54,c4_324,c6_13,B,F1,F2,F3,F4,G1,G2,G3,G4,...
    H1,H2,H3,H4,I1,I2,I3,I4,ib,jb]=gen_matrices_2d_sp7(mesh.nodes(:,1:2),...
    sort(mesh.elements')',...
    mesh.bndvtx,...
    mesh.mua,...
    mesh.mus,...
    mesh.g,...
    f1,...
    f2,...
    f3,...
    f4,...
    g1,...
    g2,...
    g3,...
    g4,...
    h1,...
    h2,...
    h3,...
    h4,...
    i1,...
    i2,...
    i3,...
    i4,...
    mesh.c,...
    omega);

nz_i=length(nonzeros(i));
K1 = sparse(i(1:nz_i),j(1:nz_i),k1(1:nz_i));
K3 = sparse(i(1:nz_i),j(1:nz_i),k3(1:nz_i));
K5 = sparse(i(1:nz_i),j(1:nz_i),k5(1:nz_i));
K7= sparse(i(1:nz_i),j(1:nz_i),k7(1:nz_i));
C = sparse(i(1:nz_i),j(1:nz_i),c(1:nz_i));
C23 = sparse(i(1:nz_i),j(1:nz_i),c23(1:nz_i));
C49 = sparse(i(1:nz_i),j(1:nz_i),c49(1:nz_i));
C64 = sparse(i(1:nz_i),j(1:nz_i),c64(1:nz_i));
C8 = sparse(i(1:nz_i),j(1:nz_i),c8(1:nz_i));
C16 = sparse(i(1:nz_i),j(1:nz_i),c16(1:nz_i));
C32 = sparse(i(1:nz_i),j(1:nz_i),c32(1:nz_i));
C128 = sparse(i(1:nz_i),j(1:nz_i),c128(1:nz_i));
C256 = sparse(i(1:nz_i),j(1:nz_i),c256(1:nz_i));
C2_59 = sparse(i(1:nz_i),j(1:nz_i),c2_59(1:nz_i));
C2_16 = sparse(i(1:nz_i),j(1:nz_i),c2_16(1:nz_i));
C2_49 = sparse(i(1:nz_i),j(1:nz_i),c2_49(1:nz_i));
C2_8 = sparse(i(1:nz_i),j(1:nz_i),c2_8(1:nz_i));
C2_32 = sparse(i(1:nz_i),j(1:nz_i),c2_32(1:nz_i));
C2_64 = sparse(i(1:nz_i),j(1:nz_i),c2_64(1:nz_i));
C4_9 = sparse(i(1:nz_i),j(1:nz_i),c4_9(1:nz_i));
C4_54 = sparse(i(1:nz_i),j(1:nz_i),c4_54(1:nz_i));
C4_324 = sparse(i(1:nz_i),j(1:nz_i),c4_324(1:nz_i));
C6_13 = sparse(i(1:nz_i),j(1:nz_i),c6_13(1:nz_i));
B = sparse(i(1:nz_i),j(1:nz_i),B(1:nz_i));

[i,j]=size(C);
nz_i = length(nonzeros(ib));
F1 = sparse(ib(1:nz_i),jb(1:nz_i),F1(1:nz_i));
F2 = sparse(ib(1:nz_i),jb(1:nz_i),F2(1:nz_i));
F3 = sparse(ib(1:nz_i),jb(1:nz_i),F3(1:nz_i));
F4 = sparse(ib(1:nz_i),jb(1:nz_i),F4(1:nz_i));
G1 = sparse(ib(1:nz_i),jb(1:nz_i),G1(1:nz_i));
G2 = sparse(ib(1:nz_i),jb(1:nz_i),G2(1:nz_i));
G3 = sparse(ib(1:nz_i),jb(1:nz_i),G3(1:nz_i));
G4 = sparse(ib(1:nz_i),jb(1:nz_i),G4(1:nz_i));
H1 = sparse(ib(1:nz_i),jb(1:nz_i),H1(1:nz_i));
H2 = sparse(ib(1:nz_i),jb(1:nz_i),H2(1:nz_i));
H3 = sparse(ib(1:nz_i),jb(1:nz_i),H3(1:nz_i));
H4 = sparse(ib(1:nz_i),jb(1:nz_i),H4(1:nz_i));
I1 = sparse(ib(1:nz_i),jb(1:nz_i),I1(1:nz_i));
I2 = sparse(ib(1:nz_i),jb(1:nz_i),I2(1:nz_i));
I3 = sparse(ib(1:nz_i),jb(1:nz_i),I3(1:nz_i));
I4 = sparse(ib(1:nz_i),jb(1:nz_i),I4(1:nz_i));
    
elseif mesh.dimension==3
  [i,j,k1,k3,k5,k7,c,c23,c49,c64,c8,c16,c32,c128,c256,c2_59,c2_16,c2_49,...
    c2_8,c2_32,c2_64,c4_9,c4_54,c4_324,c6_13,B,F1,F2,F3,F4,G1,G2,G3,G4,...
    H1,H2,H3,H4,I1,I2,I3,I4,ib,jb]=gen_matrices_3d_sp7(mesh.nodes,...
    sort(mesh.elements')',...
    mesh.bndvtx,...
    mesh.mua,...
    mesh.mus,...
    mesh.g,...
    f1,...
    f2,...
    f3,...
    f4,...
    g1,...
    g2,...
    g3,...
    g4,...
    h1,...
    h2,...
    h3,...
    h4,...
    i1,...
    i2,...
    i3,...
    i4,...
    mesh.c,...
   omega);

nz_i=length(nonzeros(i));
K1 = sparse(i(1:nz_i),j(1:nz_i),k1(1:nz_i));
K3 = sparse(i(1:nz_i),j(1:nz_i),k3(1:nz_i));
K5 = sparse(i(1:nz_i),j(1:nz_i),k5(1:nz_i));
K7= sparse(i(1:nz_i),j(1:nz_i),k7(1:nz_i));
C = sparse(i(1:nz_i),j(1:nz_i),c(1:nz_i));
C23 = sparse(i(1:nz_i),j(1:nz_i),c23(1:nz_i));
C49 = sparse(i(1:nz_i),j(1:nz_i),c49(1:nz_i));
C64 = sparse(i(1:nz_i),j(1:nz_i),c64(1:nz_i));
C8 = sparse(i(1:nz_i),j(1:nz_i),c8(1:nz_i));
C16 = sparse(i(1:nz_i),j(1:nz_i),c16(1:nz_i));
C32 = sparse(i(1:nz_i),j(1:nz_i),c32(1:nz_i));
C128 = sparse(i(1:nz_i),j(1:nz_i),c128(1:nz_i));
C256 = sparse(i(1:nz_i),j(1:nz_i),c256(1:nz_i));
C2_59 = sparse(i(1:nz_i),j(1:nz_i),c2_59(1:nz_i));
C2_16 = sparse(i(1:nz_i),j(1:nz_i),c2_16(1:nz_i));
C2_49 = sparse(i(1:nz_i),j(1:nz_i),c2_49(1:nz_i));
C2_8 = sparse(i(1:nz_i),j(1:nz_i),c2_8(1:nz_i));
C2_32 = sparse(i(1:nz_i),j(1:nz_i),c2_32(1:nz_i));
C2_64 = sparse(i(1:nz_i),j(1:nz_i),c2_64(1:nz_i));
C4_9 = sparse(i(1:nz_i),j(1:nz_i),c4_9(1:nz_i));
C4_54 = sparse(i(1:nz_i),j(1:nz_i),c4_54(1:nz_i));
C4_324 = sparse(i(1:nz_i),j(1:nz_i),c4_324(1:nz_i));
C6_13 = sparse(i(1:nz_i),j(1:nz_i),c6_13(1:nz_i));
B = sparse(i(1:nz_i),j(1:nz_i),B(1:nz_i));
F1 = sparse(i(1:nz_i),j(1:nz_i),F1(1:nz_i));
F2 = sparse(i(1:nz_i),j(1:nz_i),F2(1:nz_i));
F3 = sparse(i(1:nz_i),j(1:nz_i),F3(1:nz_i));
F4 = sparse(i(1:nz_i),j(1:nz_i),F4(1:nz_i));
G1 = sparse(i(1:nz_i),j(1:nz_i),G1(1:nz_i));
G2 = sparse(i(1:nz_i),j(1:nz_i),G2(1:nz_i));
G3 = sparse(i(1:nz_i),j(1:nz_i),G3(1:nz_i));
G4 = sparse(i(1:nz_i),j(1:nz_i),G4(1:nz_i));
H1 = sparse(i(1:nz_i),j(1:nz_i),H1(1:nz_i));
H2 = sparse(i(1:nz_i),j(1:nz_i),H2(1:nz_i));
H3 = sparse(i(1:nz_i),j(1:nz_i),H3(1:nz_i));
H4 = sparse(i(1:nz_i),j(1:nz_i),H4(1:nz_i));
I1 = sparse(i(1:nz_i),j(1:nz_i),I1(1:nz_i));
I2 = sparse(i(1:nz_i),j(1:nz_i),I2(1:nz_i));
I3 = sparse(i(1:nz_i),j(1:nz_i),I3(1:nz_i));
I4 = sparse(i(1:nz_i),j(1:nz_i),I4(1:nz_i));

end

clear i* j* k* c* omega nz_i 

% Add complex component to absorption moments due to 
% frequency dependence
C=C+B;
C23=C23+B;
C49=C49+B;
C64=C64+B;
C8=C8+B;
C16=C16+B;
C32=C32 +B;
C128 = C128+B;
C256=C256+B;
C2_59=C2_59+B;
C2_16=C2_16+B;
C2_49=C2_49+B;
C2_8=C2_8+B;
C2_32=C2_32+B;
C2_64=C2_64+B;
C4_9=C4_9+B;
C4_54 = C4_54+B;
C4_324 = C4_324+B;
C6_13 = C6_13+B;

M1=K1+C+F1;
M2=K3+C49+C2_59+F2;
M3=K5+C64+C2_16+C4_9+F3;
M4=K7+C256+C2_64+C4_324+C6_13+F4;


%==========================================

if isfield(mesh,'ident') == 1
  disp('Modifying for refractive index')
  M = bound_int(MASS,mesh);
  MASS = M;
  clear M
end

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

%i=nnodes;
qvec=[qvec;(-2/3)*qvec;(8/15)*qvec;(-16/35)*qvec];

MASS=[M1 -(C23-G1) +(C8-H1) -(C16-I1);...
     -(C23-G2) M2 -(C16-C2_49-H2) (C32+C2_8-I2);...
     (C8-G3) -(C16+C2_49-H3) M3 -(C128+C2_32+C4_54-I3);...
     (-C16-G4) (C32+C2_8-H4) -(C128+C2_32+C4_54-I4) M4];
 
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

% %=======================================

phi_all=get_field(MASS_opt,mesh,Q_opt);

% Re-order elements
%phi_all=phi_all(invsort,:);
% Extract various orders of phi
data.phi1=phi_all(1:nvtx,:);
data.phi2=phi_all((nvtx+1):(2*nvtx),:);
data.phi3=phi_all(((2*nvtx)+1):3*nvtx,:);
data.phi4=phi_all(((3*nvtx)+1):4*nvtx,:);

data.phi = data.phi1 - (2/3).*data.phi2+...
        (8/15).*data.phi3-(16/35).*data.phi4;
clear phi*

% Calculate boundary data
[data.complex1]=get_boundary_data(mesh,data.phi1);
[data.complex2]=get_boundary_data(mesh,data.phi2);
[data.complex3]=get_boundary_data(mesh,data.phi3);
[data.complex4]=get_boundary_data(mesh,data.phi4);
[data.complex]=get_boundary_data(mesh,data.phi);
data.link = mesh.link;

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

% Phi 4
data.amplitude4 = abs(data.complex4);
data.phase4 = atan2(imag(data.complex4),...
		   real(data.complex4));
data.phase4(find(data.phase4<0)) = data.phase4(find(data.phase4<0)) + (2*pi);
data.phase4 = data.phase4*180/pi;
data.paa4 = [data.amplitude4 data.phase4];

% Total
data.amplitude = abs(data.complex);
data.phase = atan2(imag(data.complex),...
		   real(data.complex));
data.phase(find(data.phase<0)) = data.phase(find(data.phase<0)) + (2*pi);
data.phase = data.phase*180/pi;
data.paa = [data.amplitude data.phase];
