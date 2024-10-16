function varargout=random_rotation(max_angle,max_trans,mesh)
% max_angle ; max angle (in degrees) of rotation about X,Y,Z 
% max_trans

angles=rand(3,1).*max_angle(:);
trans=rand(3,1).*max_trans(:);

T=eul2tform(angles'*pi/180,'XYZ');
T(1:3,4)=trans;

varargout{1}=T;

if(nargin>2)
    mesh=nirs.registration.applytform(mesh,T);
end
if(nargout==2)
    varargout{2}=mesh;
end