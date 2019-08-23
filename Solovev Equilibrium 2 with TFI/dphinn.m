%fungsi mencari parameter bagian 3 dari Q,untuk batas,k!=m
function y = dphinn(a,b,nkr,nkz,l,m,zo,Lz,Lr)
%% menentukan phi, dphir, dan dphiz
dphirpart=0;
dphizpart=0;
for j=1:5
    C=1;
    for s=1:j
        C=C.*(-((l-2*s+2)*(l-2*s)/((m+2*s+1)*(m+2*s+2))).*(((b-zo)./a).^2));
    end
    dphirpart=dphirpart+((l-2*j)*C);
    dphizpart=dphizpart+((m+2*j+2)*C);
end
dphir=-(((b-zo)/Lz).^(l-1)).*(((a/Lr).^(m+2))/((m+1)*(m+2)))*((Lz.^2)/Lr).*(l+dphirpart);
dphiz=-(((b-zo)/Lz).^l).*(((a/Lr).^(m+1))/((m+1)*(m+2)))*((Lz.^2)/Lr).*(m+2+dphizpart);
%% menentukan dphin
dphin = dphir*nkr+dphiz*nkz;
%% menentukan Q3
y = dphin;