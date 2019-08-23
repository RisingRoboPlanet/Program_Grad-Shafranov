function [K1,K2,E1,E2] = findKE(a,b,nkr,nkz,lk,rk,zk,rk2,zk2)
k1=(quadgk(@(t) intKE(t,a,b,nkr,nkz,lk,rk,zk,rk2,zk2),-1,0,'AbsTol',1e-10));
k2=(quadgk(@(t) intKE(t,a,b,nkr,nkz,lk,rk,zk,rk2,zk2),0,1,'AbsTol',1e-10));
[K1,E1]=ellipke(k1);
[K2,E2]=ellipke(k2);
