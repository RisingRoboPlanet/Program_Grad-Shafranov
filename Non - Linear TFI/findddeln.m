function Y = findddeln(xk,yk,xk2,yk2,lk,nkx,nky,Lx,Ly,yo,l,m)
    Y=(lk/2)*quadgk(@(t) ddeln(t,xk,yk,xk2,yk2,lk,nkx,nky,Lx,Ly,yo,l,m),-1,1,'AbsTol',1e-10);
