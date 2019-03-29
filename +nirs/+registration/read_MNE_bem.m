function varargout = read_MNE_bem(bemfile,srcspace,lambda)
% Conversion code for reading MNE BEM/SrcSpace fif files

if(nargin<3)
    lambda=[690 830];
end

bem=mne_read_bem_surfaces(bemfile);
for i=1:length(bem)
    % Make the Mesh models
    mesh(i,1)=nirs.core.Mesh(bem(i).rr*1000,bem(i).tris);
    mesh(i).transparency=.1;
end

if(nargin>1)
    src=mne_read_source_spaces(srcspace);
    [~,f]=ismember(src(1).use_tris,src(1).vertno);
    v=src(1).rr(src(1).vertno,:);
    [~,f2]=ismember(src(2).use_tris,src(2).vertno);
    v2=src(2).rr(src(2).vertno,:);
    mesh(end+1)=nirs.core.Mesh([v; v2]*1000,[f; f2+length(v)]);
    mesh(end).transparency=1;
end

if(nargin<3)
    varargout{1}=mesh;
else
    % Now make the forward model
    fwdModel=nirs.forward.NirfastBEM;
    fwdModel.mesh=mesh;
    fwdModel.prop={nirs.media.tissues.skin(lambda)...
        nirs.media.tissues.bone(lambda)...
        nirs.media.tissues.water(lambda)...
        nirs.media.tissues.brain(lambda,.7,60)};
    varargout{1}=fwdModel;
end