function H = inth(t,a,b,nkr,nkz,lk,rk,zk,rk2,zk2)
r=((rk2+rk)/2)-(lk*t*nkz/2);
z=((zk2+zk)/2)+(lk*t*nkr/2);
k=sqrt(4*a*r/((r+a).^2+(z-b).^2));
[K,E]=ellipke(k);
dpsisr = (r./(2*pi*sqrt((r+a).^2+(z-b).^2))).*(K-E.*((r.^2-a^2+(z-b).^2)./((r-a).^2+(z-b).^2)));
dpsisz = ((z-b)./(2*pi*sqrt((r+a).^2+(z-b).^2))).*(K-E.*((r.^2+a^2+(z-b).^2)./((r-a).^2+(z-b).^2)));
H=(dpsisr*nkr+dpsisz*nkz)./r;
