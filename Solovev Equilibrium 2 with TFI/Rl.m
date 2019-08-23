function xyz = Rl(s)
inv = 0.32;
elon = 1.7;
tri = 0.33;
r = 1+inv*cos((pi/2*(1-s))+(pi)+asin(tri)*sin((pi/2*(1-s))+(pi)));
z = elon*inv*sin((pi/2*(1-s))+(pi));

xyz = [r;z];