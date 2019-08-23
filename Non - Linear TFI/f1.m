function Y = f1(x,psi,Ro,bp,alpha,gamma)
X=((psi));
c=nthroot(X^8,10);
Y=0.5*(0.2*(x)+(1-0.2)*(1/x))*c;