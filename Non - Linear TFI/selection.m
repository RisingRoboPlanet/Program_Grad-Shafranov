function y = selection(x1,y1,x2,y2,xp,yp)
Ax=x1-xp;
Bx=x2-xp;
Ay=y1-yp;
By=y2-yp;

Cp=Ax*By-Ay*Bx;
sgn = sign(Cp);

if(abs(Cp)>1e-6)
    Dot=Ax*Bx+Ay*By;
    A2=Ax*Ax+Ay*Ay;
    B2=Bx*Bx+By*By;
    prod = A2*B2;
    Ang=0;
    if (prod ~=0)
        Ang=acos(Dot/sqrt(prod));
        Ang=Ang*sgn;
    end
    y=Ang;
else
    y=0;
end