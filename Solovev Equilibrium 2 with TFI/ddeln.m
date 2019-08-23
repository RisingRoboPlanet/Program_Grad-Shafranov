function Y = ddeln(t,xk,yk,xk2,yk2,lk,nkx,nky,Lx,Ly,yo,l,m)
x=((xk2+xk)/2)-(lk*t*nky/2);
y=((yk2+yk)/2)+(lk*t*nkx/2);
% menentukan ddelr
p=0;
pp=0;
cs(1)=Lx^2/(2*(l+2)*(l+1));
ds(1)=Ly^2/(2*(m+2)*(m+1));

K=(m+2)/2;
for s=1:K
    if(s==1)
        ptot=cs(s)*(l+2*s)*(x/Lx).^(l+2*s-1).*((y-yo)/Ly).^(m-2*s+2)/(Lx);
        p=p+ptot;
        ptott=cs(s)*(m-2*s+2)*(x/Lx).^(l+2*s).*((y-yo)/Ly).^(m-2*s+1)/(Ly);
        pp=pp+ptott;
    else
        cs(s)=-(m-2*s+4)*(m-2*s+3)*(Lx^2)/((l+2*s)*(l+2*s-1)*Ly^2)*cs(s-1);
        ptot=cs(s)*(l+2*s)*(x/Lx).^(l+2*s-1).*((y-yo)/Ly).^(m-2*s+2)/(Lx);
        p=p+ptot;
        ptott=cs(s)*(m-2*s+2)*(x/Lx).^(l+2*s).*((y-yo)/Ly).^(m-2*s+1)/(Ly);
        pp=pp+ptott;
    end
end

K=(l+2)/2;
for s=1:K
    if(s==1)
        ptot=ds(s)*(l-2*s+2)*(x/Lx).^(l-2*s+1)/(Lx).*((y-yo)/Ly).^(m+2*s);
        p=p+ptot;
        ptott=ds(s)*(m+2*s)*(x/Lx).^(l-2*s+2)/(Ly).*((y-yo)/Ly).^(m+2*s-1);
        pp=pp+ptott;
    else
        ds(s)=-(l-2*s+4)*(l-2*s+3)*(Ly^2)/((m+2*s)*(m+2*s-1)*Lx^2)*ds(s-1);
        ptot=ds(s)*(l-2*s+2)*((x/Lx).^(l-2*s+1))/(Lx).*((y-yo)/Ly).^(m+2*s);
        p=p+ptot;
        ptott=ds(s)*(m+2*s)*(x/Lx).^(l-2*s+2)/(Ly).*((y-yo)/Ly).^(m+2*s-1);
        pp=pp+ptott;
    end
end
ddelx=p;
ddely=pp;
Y=(ddelx.*nkx+ddely.*nky);