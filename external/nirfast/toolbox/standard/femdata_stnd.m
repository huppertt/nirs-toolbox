function [data,mesh]=femdata_stnd(mesh,frequency)

% [data,mesh]=femdata_stnd(mesh,frequency)
% Calculates data (phase and amplitude) for a given
% standard mesh at a given frequency (MHz).
% outputs phase and amplitude in structure data
% and mesh information in mesh

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

% Create FEM matricex
if mesh.dimension == 2
  [i,j,s] = gen_matrices_2d(mesh.nodes(:,1:2),...
			    sort(mesh.elements')', ...
			    mesh.bndvtx,...
			    mesh.mua,...
			    mesh.kappa,...
			    mesh.ksi,...
			    mesh.c,...
			    omega);
elseif mesh.dimension ==3
  [i,j,s] = gen_matrices_3d(mesh.nodes,...
			    sort(mesh.elements')', ...
			    mesh.bndvtx,...
			    mesh.mua,...
			    mesh.kappa,...
			    mesh.ksi,...
			    mesh.c,...
			    omega);
end
%junk = length(find(i==0));
%MASS = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
%%next two lines modified by subha, 6/11/07 to make it more memory efficient
nz_i = nonzeros(i);
MASS = sparse(i(1:length(nz_i)),j(1:length(nz_i)),s(1:length(nz_i)));
clear junk i j s omega nz_i

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

% Check for distributed source
if mesh.source.distributed == 1
    qvec = sum(qvec,2);
end

% catch error in source vector
junk = sum(qvec);
junk = find(junk==0);
if ~isempty(junk)
    display(['WARNING...Check the FWHM of Sources ' num2str(junk)]);
end
clear junk

% Catch zero frequency (CW) here
if frequency == 0
  MASS = real(MASS);
  qvec = real(qvec);
end


% Calculate field for all sources
[data.phi,mesh.R]=get_field(MASS,mesh,qvec);
clear qvec;
% clear preconditionare if not a spec mesh!
if all(mesh.type=='spec') == 0
  mesh=rmfield(mesh,'R');
end


% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);
data.link = mesh.link;

% Map complex data to amplitude and phase
data.amplitude = abs(data.complex);

data.phase = atan2(imag(data.complex),...
		   real(data.complex));
data.phase(data.phase<0) = data.phase(data.phase<0) + (2*pi);
data.phase = data.phase*180/pi;

data.paa = [data.amplitude data.phase];
