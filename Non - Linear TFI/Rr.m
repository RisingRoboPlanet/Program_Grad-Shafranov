function xyz = Rr(s)
inv = 0.32;
elon = 1.7;
tri = 0.33;
r = 1+inv*cos((pi/2*s)+asin(tri)*sin((pi/2*s)));
z = elon*inv*sin((pi/2*s));

xyz = [r;z];