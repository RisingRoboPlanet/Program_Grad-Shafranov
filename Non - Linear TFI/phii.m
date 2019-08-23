% fungsi untuk mencari parameter phi
function y = phii(a,b,l,m,zo,Lz,Lr,c)
phipart=0;
for k=1:5
    C=1;
    for s = 1:k
        C=C*(-(l-2*s+2)*(l-2*s)/((m+2*s+1)*(m+2*s+2))*((b-zo)/a)^2);
    end
    phipart=phipart+C;
end
phi = -((Lz^2)*((a/Lr)^l)*((((b-zo)/Lz)^(m+2)))/((m+1)*(m+2)))*(1+phipart);
y=phi;