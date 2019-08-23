function xyz = Rt(s)
inv = 0.32;
elon = 1.7;
tri = 0.33;
r = 1+inv*cos((pi/2*(1-s))+(pi/2)+asin(tri)*sin((pi/2*(1-s))+(pi/2)));
z = elon*inv*sin((pi/2*(1-s))+(pi/2));

xyz = [r;z];