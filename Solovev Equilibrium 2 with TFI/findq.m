function Y = findq(nkz,lk,rk,rk2)
Y=(1/(1.25663706e-6))*(lk/2)*2*pi;
%quadgk(@(t) intq(t,nkz,lk,rk,rk2),-1,1,'AbsTol',1e-10);