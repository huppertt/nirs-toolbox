function [data,mesh]=femdata_fl_tr(mesh,t,dt,flag)

% [data,mesh]=femdata_fl_tr(mesh,t,dt);
% Calculates Time Resolved (TR) data, i.e. TPSF at detectors, 
% for Fluorescence meshes (i.e. excitation and FL emission)
% for total time t at time interval dt in seconds

% e.g. data = femdata_fl_tr('circle2000_86_fl',1e-9,10e-12);

% in order to save memory, only one source will be run at each time, 
% i.e. multiple sources are not vectorised and each source will be taken 
% in turn to run within the code. 

% if last argument is given as 1, it will also output internal fluence

% Part of NIRFAST package
% Excitation by H Dehghani 2008, Emission by Q ZHU April 2010

% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

% flag to see if internal fluence is needed as output
if nargin == 4
    flag = 1;
else
    flag = 0;
end

% Create FEM matricex
% note dt for this routine needed in pico-seconds
%Excitation
if mesh.dimension == 2
    [i,j,s1,s2] = gen_matrices_2d_TR(mesh.nodes(:,1:2),...
        sort(mesh.elements')', ...
        mesh.bndvtx,...
        mesh.muax+mesh.muaf,...
        mesh.kappax,...
        mesh.ksi,...
        dt*mesh.c);
elseif mesh.dimension ==3
    [i,j,s1,s2] = gen_matrices_3d_TR(mesh.nodes,...
        sort(mesh.elements')', ...
        mesh.bndvtx,...
        mesh.muax+mesh.muaf,...
        mesh.kappax,...
        mesh.ksi,...
        dt*mesh.c);
end
nz_i = nonzeros(i);
MASSx1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),s1(1:length(nz_i)));
MASSx2 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),s2(1:length(nz_i)));
clear junk i j s1 s2 nz_i

%Emission
if mesh.dimension == 2
    [i,j,s1,s2] = gen_matrices_2d_TR(mesh.nodes(:,1:2),...
        sort(mesh.elements')', ...
        mesh.bndvtx,...
        mesh.muam,...
        mesh.kappam,...
        mesh.ksi,...
        dt*mesh.c);
elseif mesh.dimension ==3
    [i,j,s1,s2] = gen_matrices_3d_TR(mesh.nodes,...
        sort(mesh.elements')', ...
        mesh.bndvtx,...
        mesh.muam,...
        mesh.kappam,...
        mesh.ksi,...
        dt*mesh.c);
end
nz_i = nonzeros(i);
MASSm1 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),s1(1:length(nz_i)));
MASSm2 = sparse(i(1:length(nz_i)),j(1:length(nz_i)),s2(1:length(nz_i)));
clear junk i j s1 s2 nz_i

% If the fn.ident exists, then we must modify the FEM matrices to
% account for refractive index mismatch within internal boundaries
if isfield(mesh,'ident') == 1
    disp('Modifying for refractive index')
    M = bound_int(MASS,mesh);
    MASS = M;
    clear M
end

% Calculate the RHS (the source vectors. For simplicity, we are
% just going to use a Gaussian Source, The width of the Gaussian is
% changeable (last argument). The source is assumed to have a
% complex amplitude of complex(cos(0.15),sin(0.15));

% Now calculate source vector
% NOTE last term in mex file 'qvec' is the source FWHM
[nnodes]=size(mesh.nodes,1);
[nsource]=size(mesh.source.coord,1);
if max(mesh.source.num)~=nsource
    disp('Error in maximun numbering source')
    errordlg('The maximun source numbering is wrong!','NIRFAST Error');
    error('The maximun source numbering is wrong!');
end

qvec = spalloc(nnodes,nsource,nsource*100);
if mesh.dimension == 2
    for j = 1:nsource 
        i = mesh.source.num(j);
        if mesh.source.fwhm(j) == 0
            qvec(:,i) = gen_source_point(mesh,mesh.source.coord(j,1:2));
        else
            qvec(:,i) = gen_source(mesh.nodes(:,1:2),...
                sort(mesh.elements')',...
                mesh.dimension,...
                mesh.source.coord(j,1:2),...
                mesh.source.fwhm(j));
        end
    end
elseif mesh.dimension == 3
    for j = 1:nsource 
        i = mesh.source.num(j);
        if mesh.source.fwhm(j) == 0
            qvec(:,i) = gen_source_point(mesh,mesh.source.coord(j,1:3));
        else
            qvec(:,i) = gen_source(mesh.nodes,...
                sort(mesh.elements')',...
                mesh.dimension,...
                mesh.source.coord(j,:),...
                mesh.source.fwhm(j));
        end
    end
end
if (j)~=nsource
    disp('Error in sequencing qvec')
    errordlg('Check the sequence of the qvec!','NIRFAST Error');
    error('Check the sequence of the qvec!');
end

% only need real part
qvec = real(qvec);

clear junk i;

% catch error in source vector
junk = sum(qvec);
junk = find(junk==0);
if ~isempty(junk)
    display(['WARNING...Check the FWHM of Sources ' num2str(junk)]);
end
clear junk

% calculte size of data output
ndt = length(0:dt:t);
nm = length(find(mesh.link(:,3)==1));
data.tpsfx = zeros(nm,ndt);
data.tpsfm = zeros(nm,2*ndt-1);
data.phix = zeros(nnodes,nsource,ndt);
data.phim = zeros(nnodes,nsource,2*ndt-1);


for iii = 1:nsource
    %%%%% Calculate time resolved field for one source for Excitation %%%%%
    [phi_ex,mesh.R]=get_field_tpsf_ex(MASSx1,MASSx2,mesh,qvec(:,iii),dt,t);
    
    %memory allocation for source of fluorescent
    qvec_m = zeros(nnodes,2*ndt-1);
    % predefine variable for fluorescent life-time convolution
    cf1 = zeros(1,ndt);
    
    %generation of fluorescent sources
    mf2 = (mesh.eta.*mesh.muaf).*(dt./mesh.tau);
    % need to ensure will work with zero lifetime
    ind = find(mesh.tau==0);
    mf2(ind) = (mesh.eta(ind).*mesh.muaf(ind));
    
    % create rhs for emission part
    for i = 1:nnodes
        if mf2(i) ~= 0 && ~isnan(mf2(i))
            cf1(1,:) = exp((-(0:ndt-1).*dt)/mesh.tau(i));
            % need to ensure will work with zero lifetime
            if mesh.tau(i) == 0
                cf1(1,1) = 1;
                cf1(1,2:end) = 0;
            end
            %convolution, for fluorescent sources generation
            qvec_m(i,:) = conv(phi_ex(i,:),cf1(1,:)).*mf2(i);
        end
    end
    
    data.phix(:,iii,:) = phi_ex;

    clear phi_ex cf1
    
    %%%%% Calculate time resolved field for one source for Emission %%%%%%%
    [phi_em,mesh.R]=get_field_tpsf_em(MASSm1,MASSm2,mesh,qvec_m,dt,t);
    clear qvec_m;
    

    data.phim(:,iii,:) = phi_em;
    clear phi_em 
end


% Calculate boundary data for Excitation %%%%%%%%%%%%%
tpsf_tem = zeros(size(mesh.link,1),ndt);
for i = 1 : ndt
        tpsf_tem(:,i)= ...
            get_boundary_data(mesh,data.phix(:,:,i));         
end
rw = find(~isnan(tpsf_tem(:,1)));
data.tpsfx = tpsf_tem(rw,:);
clear tpsf_tem

    
% Calculate boundary data for Emission %%%%%%%%%%%%%
tpsf_tem = zeros(size(mesh.link,1),2*ndt-1);
for i = 1 : 2*ndt-1
        tpsf_tem(:,i)= ...
            get_boundary_data(mesh,data.phim(:,:,i));         
end
rw = find(~isnan(tpsf_tem(:,1)));
data.tpsfm = tpsf_tem(rw,:);
data.paa = data.tpsfx;
data.wv = 0:dt:t;
data.link = mesh.link;
clear tpsf_tem
     
            
if flag == 0
    data = rmfield(data,'phix');
    data = rmfield(data,'phim');
end

clear qvec dt flag i ii iii k ndt nm nnodes nsource t MASS*


%%%%%%%%%%%%%%time resolved field for one source for Emission
function [phi,R]=get_field_tpsf_em(MASS1,MASS2,mesh,qvec_m,dt,t)

[nnodes,n_conv_dt] = size(qvec_m);
phi=zeros(nnodes,n_conv_dt);
msg=[];
flag = 0;

% calculate the preconditioner
if isfield(mesh,'R') == 0 % for spec, this may have been calculated
    if length(mesh.nodes) >= 3800
        R = cholinc(MASS1,1e-3);
    elseif length(mesh.nodes) < 3800
        R = [];
    end
else
    R = mesh.R;
end

% calculating 'phi' for the first time step
if nnodes >= 3800
%     tic
    [x,flag] = bicgstab(MASS1,qvec_m(:,1),1e-12,100,R',R);
    msg = [msg flag];
    phi(:,1) = x;
    k = 2;
    for j = dt : dt : t
        [x,flag]=bicgstab(MASS1,(qvec_m(:,k)-MASS2*phi(:,k-1)),1e-12,100,R',R);
        msg = [msg flag];
        phi(:,k) = x;
        k = k + 1;
    end

else
    phi(:,1) = MASS1\qvec_m(:,1);
    k = 2;
    for j = dt : dt : t
        phi(:,k) = MASS1\(qvec_m(:,k)-MASS2*phi(:,k-1));
        k = k + 1;
    end
end

if isempty(msg)
%   disp('Used backslash!')
elseif any(msg==1)
    disp('some solutions did not converge')
    errordlg('Some solutions did not converge; this could be caused by noisy/bad data','NIRFAST Error');
    error('Some solutions did not converge; this could be caused by noisy/bad data');
elseif any(msg==2)
    disp('some solutions are unusable')
    errordlg('Some solutions are unusable; this could be caused by noisy/bad data','NIRFAST Error');
    error('Some solutions are unusable; this could be caused by noisy/bad data');
elseif any(msg==3)
    disp('some solutions from stagnated iterations')
    errordlg('Some solutions are unusable; this could be caused by noisy/bad data','NIRFAST Error');
    error('Some solutions are unusable; this could be caused by noisy/bad data');
end
clear x msg flag


%%%%%%%%%%%%%%time resolved field for one source for Excitation
function [phi,R]=get_field_tpsf_ex(MASS1,MASS2,mesh,qvec,dt,t)

[nnodes,nsource]=size(qvec);
ndt = length(0:dt:t);
phi=zeros(nnodes,ndt);
msg=[];
flag = 0;

% calculate the preconditioner

if isfield(mesh,'R') == 0 % for spec, this may have been calculated
    if nnodes >= 3800
        R = cholinc(MASS1,1e-3);
    elseif nnodes < 3800
        R = [];
    end
else
    R = mesh.R;
end

if nnodes >= 3800

    [x,flag] = bicgstab(MASS1,qvec,1e-12,100,R',R);
    msg = [msg flag];
    phi(:,1) = x;
    k = 2;
    for j = dt : dt : t
        [x,flag]=bicgstab(MASS1,(-MASS2*phi(:,k-1)),1e-12,100,R',R);
        msg = [msg flag];
        phi(:,k) = x;
        k = k + 1;
    end

else
    phi(:,1) = MASS1\qvec;
    k = 2;
    for j = dt : dt : t
        phi(:,k) = MASS1\(-MASS2*phi(:,k-1));
        k = k + 1;
    end
end

if isempty(msg)
    %   disp('Used backslash!')
elseif any(msg==1)
    disp('some solutions did not converge')
    errordlg('Some solutions did not converge; this could be caused by noisy/bad data','NIRFAST Error');
    error('Some solutions did not converge; this could be caused by noisy/bad data');
elseif any(msg==2)
    disp('some solutions are unusable')
    errordlg('Some solutions are unusable; this could be caused by noisy/bad data','NIRFAST Error');
    error('Some solutions are unusable; this could be caused by noisy/bad data');
elseif any(msg==3)
    disp('some solutions from stagnated iterations')
    errordlg('Some solutions are unusable; this could be caused by noisy/bad data','NIRFAST Error');
    error('Some solutions are unusable; this could be caused by noisy/bad data');
end
clear x msg flag