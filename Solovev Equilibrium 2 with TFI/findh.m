function F = findh(a,b,nkr,nkz,lk,rk,zk,rk2,zk2)
F1=(lk/2)*quadgk(@(t) inth(t,a,b,nkr,nkz,lk,rk,zk,rk2,zk2),-1,-0.00001,'AbsTol',1e-10);
F2=(lk/2)*quadgk(@(t) inth(t,a,b,nkr,nkz,lk,rk,zk,rk2,zk2),0.00001,1,'AbsTol',1e-10);
F=F1+F2;
