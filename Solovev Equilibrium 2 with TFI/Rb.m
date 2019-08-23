function xyz = Rb(s)
inv = 0.32;
elon = 1.7;
tri = 0.33;
r = 1+inv*cos((pi/2*s)+(pi/2*3)+asin(tri)*sin((pi/2*s)+(pi/2*3)));
z = elon*inv*sin((pi/2*s)+(pi/2*3));

xyz = [r;z];

