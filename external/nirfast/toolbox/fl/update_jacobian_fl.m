function [J,data]=update_jacobian_fl(Jin,mesh,frequency,datax,MASS_m)

% Calculates the fluroescence forward field with the current value of muaf
% and calculates the Jacobian by dividing the "pre-Jacobian", Jin, by the
% fluorescence field.
% This is required since we reconstruct on the log of the data, and
%  therefore a 1/phifl at each iteration must be used at each iteration.
%  The main part of the Jacobian, calculated with the direct excitation and
%  adjoint emission fields, does not change with each iteration and is
%  therefore calculated before the iteration loop and input as Jin here.

% Inputs
% mesh - the fwd_mesh you are using to calculate the Jacobian
% ferquency - alwasy set = 0 for fluor.
% datax - excitation field calculated previously in the form datax.phi
% MASS_m (optinal) - may use this if previously calculated.

% Output
% J is the jacobian for the current iteration
% data.amplitudem is the relevant forward fluorescence emission data

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

% Create Emission FEM Matrix if not an input
if ~exist('MASS_m','var')
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
        M = bound_int(MASS_m,mesh);
        MASS_m = M;
        clear M
    end
end

% Calculate the RHS (the source vectors) for the Emission.
source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);
%[nsource,junk]=size(source);

ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
nsource = length(source);

%qvec = zeros(nnodes,nsource);
qvec = spalloc(nnodes,nsource,nsource*100);

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

% Update the Emission jacobian
data2 = data;
ind = data.link(:,3) == 0;
data2.complexm(ind,:)=[];
[nsd, msd] = size(mesh.link);

k = 1;
for i = 1 : nsd
    if mesh.link(i,3) == 1
        J.completem(k,:) = ...
            real(Jin.complexm(k,1:end/2)./data2.complexm(k));
        %this used to be outside the loop but then was getting too high on
        %the k and the values of Jin didn't exist
         k = k + 1;
    %   ??? Index exceeds matrix dimensions.
    % Error in ==> update_jacobian_fl at 134
    %         J.completem(k,:) = ... 
    end
%     k = k + 1;
end