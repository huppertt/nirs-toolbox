function Ylm = sphericalharmonics(nodes,l,m)

[theta,phi,R] = cart2sph(nodes(:,1),...
    nodes(:,2),nodes(:,3));

Pnm = legendre(l,cos(theta));

if(m<0)    
    Ylm = sqrt(2)*sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m))))*...
        Pnm(abs(m)+1,:)'.*sin(phi*abs(m));
elseif(m==0)
    Ylm = sqrt((2*l+1)/(4*pi))*...
        Pnm(1,:)';
else
   Ylm = sqrt(2)*sqrt((2*l+1)*factorial(l-m)/(4*pi*factorial(l+m)))*...
        Pnm(m+1,:)'.*cos(phi*m);
end