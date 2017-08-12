function [fnew,facenew,vertexnew] = semireg_sphere_resample(vertex,face,f,Jnew)
% this assumes that the orginal vertex was a sphere

if(size(vertex,1)~=3)
    vertex=vertex';
    face=face';
end

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
    warning('unknown mesh size')
    options.base_mesh='ico';
    J=min(find(ico>length(vertex)));
end

%[V,F] = compute_semiregular_sphere(J,options);
[Vnew,Fnew] = compute_semiregular_sphere(Jnew,options);

%[th,phi,R]=cart2sph(V{end}(1,:),V{end}(2,:),V{end}(3,:));
[th2,phi2,R]=cart2sph(vertex(1,:),vertex(2,:),vertex(3,:));
[th3,phi3,R]=cart2sph(Vnew{end}(1,:),Vnew{end}(2,:),Vnew{end}(3,:));

f=double(f);
face=double(face);
vertex=double(vertex);   
if(J>Jnew)
    f = perform_mesh_smoothing(face,vertex,f,struct('niter_averaging',J-Jnew));
end
for i=1:size(f,2)
    disp(i)
    ff=scatteredInterpolant(double(th2'),double(phi2'),f(:,i));
    fnew(:,i) =ff(double(th3),double(phi3))';
end

% 
% f=ff(double(th),double(phi))';
% 
% fw = perform_wavelet_mesh_transform(V,F, f, 1);
% 
% fwnew=zeros(length(Vnew{end}),1);
% if(J>Jnew)
%     fwnew=fw(1:length(fwnew));
% else
%     fwnew(1:length(fw))=fw;
% end
% fnew = perform_wavelet_mesh_transform(Vnew,Fnew, fwnew, -1);
% 
facenew = Fnew{end}';
vertexnew = Vnew{end}';