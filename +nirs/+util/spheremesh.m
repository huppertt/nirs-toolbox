function mesh = spheremesh(J)
% This makes a spherical mesh using the wavelet tools

if(nargin==0)
    J=5;
end

[vertex,face] = compute_semiregular_sphere(J);

mesh=nirs.core.Mesh;
mesh.nodes=vertex';
mesh.faces=face';

end
