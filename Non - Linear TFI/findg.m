 function Y = findg(a,b,nkx,nky,lk,xk,yk,xk2,yk2)
 Y1=(lk/2)*(quadgk(@(t) intg(t,a,b,nkx,nky,lk,xk,yk,xk2,yk2),-1,-0.0001,'AbsTol',1e-10));
 Y2=(lk/2)*(quadgk(@(t) intg(t,a,b,nkx,nky,lk,xk,yk,xk2,yk2),0.0001,1,'AbsTol',1e-10));
 Y=Y1+Y2;

