function [vertexnew,facenew,fnew] = semireg_sphere_resample(vertex,face,f,Jnew)
% this assumes that the orginal vertex was a sphere

ico = [12 42 162 642 2562 10242 40962 163842];
oct = [6 18 66 258 1026 4098 16386 65538];

options.keep_subdivision=true;
if(ismember(length(vertex),ico))
    options.base_mesh='ico';
    J=find(ico==length(vertex));
elseif(ismember(length(vertex),oct))
    options.base_mesh='oct';
    J=find(oct==length(vertex));
else
    error('unknown mesh size');
end

[V,F] = compute_semiregular_sphere(J,options);
[Vnew,Fnew] = compute_semiregular_sphere(Jnew,options);

[th,phi,R]=cart2sph(V{end}(1,:),V{end}(2,:),V{end}(3,:));
[~, idx]=sortrows([th; phi]',[1 2]);

[th,ph,R]=cart2sph(vertex(1,:),vertex(2,:),vertex(3,:));
[~, idx2]=sortrows([th; phi]',[1 2]);

f(idx)=f(idx2);
fw = perform_wavelet_mesh_transform(V,F, f, 1);

fwnew=zeros(length(Vnew{end}),1);
if(J>Jnew)
    fwnew=fw(1:length(fwnew));
else
    fwnew(1:length(fw))=fw;
end
fnew = perform_wavelet_mesh_transform(Vnew,Fnew, fwnew, -1);

facenew = Fnew{end};
vertexnew = Vnew{end};