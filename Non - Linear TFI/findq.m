function Y = findq(nkz,lk,rk,rk2)
Y=(lk/2)*quadgk(@(t) intq(t,nkz,lk,rk,rk2),-1,1,'AbsTol',1e-10);