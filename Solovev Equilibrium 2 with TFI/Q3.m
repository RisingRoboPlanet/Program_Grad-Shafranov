%fungsi mencari parameter bagian 3 dari Q,untuk batas,k!=m
function Y = Q3(t,a,b,nkx,nky,lk,xk,yk,xk2,yk2,l,m,yo,Ly,Lx)
x=((xk2+xk)/2)-(lk*t*nky/2);
y=((yk2+yk)/2)+(lk*t*nkx/2);
k=sqrt((4*a*x)./((x+a).^2+(y-b).^2));
[K,E]=ellipke(k);
% K=ellipticK(k.^2);
% E=ellipticE(k.^2);
psis = (sqrt(a*x)./(pi.*k)).*((1-((k.^2)/2)).*K-E);

%% menentukan phi, dphir, dan dphiz
dphirpart=0;
phipart=0;
dphizpart=0;
for j=1:5
    C=1;
    for s=1:j
        C=C.*(-((l-2*s+2)*(l-2*s)/((m+2*s+1)*(m+2*s+2))).*(((y-yo)./x).^2));
    end
    dphirpart=dphirpart+((l-2*j)*C);
    phipart=phipart+C;
    dphizpart=dphizpart+((m+2*j+2)*C); 
end
dphir=-((x/Lx).^(l-1)).*((((y-yo)/Ly).^(m+2))/((m+1)*(m+2)))*((Ly.^2)/Lx).*(l+dphirpart);
dphiz=-((x/Lx).^l).*((((y-yo)/Ly).^(m+1))/((m+1)*(m+2)))*((Ly.^2)/Ly).*(m+2+dphizpart);
phi = -(Ly^2)*((x./Lx).^l).*((((y-yo)/Ly).^(m+2))/((m+1)*(m+2))).*(1+phipart);
%% menentukan dphin
dphin = dphir*nkx+dphiz*nky;

%% menentukan dpsisn
dpsisr = (x./(2*pi*sqrt(((x+a).^2)+((y-b).^2)))).*(K-E.*(((x.^2)-(a^2)+((y-b).^2))./(((x-a).^2)+((y-b).^2))));
dpsisz = ((y-b)./(2*pi*sqrt(((x+a).^2)+((y-b).^2)))).*(K-E.*(((x.^2)+(a^2)+((y-b).^2))./(((x-a).^2)+((y-b).^2))));
dpsisn=(dpsisr*nkx+dpsisz*nky);

%% menentukan Q3
Y= (psis./x).*dphin-(phi./x).*dpsisn;