function Y = findQ3(a,b,nkx,nky,lk,xk,yk,xk2,yk2,l,m,zo,Ly,Lx)
Y1=(lk/2)*(quadgk(@(t) Q3(t,a,b,nkx,nky,lk,xk,yk,xk2,yk2,l,m,zo,Ly,Lx),-1,-0.00001,'AbsTol',1e-10));
Y2=(lk/2)*(quadgk(@(t) Q3(t,a,b,nkx,nky,lk,xk,yk,xk2,yk2,l,m,zo,Ly,Lx),0.00001,1,'AbsTol',1e-10));
Y=Y1+Y2;