% fungsi untuk mencari parameter Q k!=m
function Y = Q2(a,b,l,m,yo,Ly,Lx,c)
phipart=0;
for k=1:5
    C=1;
    for s = 1:k
        C=C*(-(l-2*s+2)*(l-2*s)/((m+2*s+1)*(m+2*s+2))*((b-yo)/a)^2);
    end
    phipart=phipart+C;
end
phi = -((Ly^2)*((a/Lx)^l)*((((b-yo)/Ly)^(m+2)))/((m+1)*(m+2)))*(1+phipart);
Y=c*phi;
