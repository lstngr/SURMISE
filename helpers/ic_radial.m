function [x,v,m] = ic_radial(sz,xl,xh,vrad,vstd,ml,mh)

dx = (xh - xl)/2;
r  = unifrnd(0,dx,sz,1);
vr = normrnd(vrad,vstd,sz,1);
a  = unifrnd(0,2*pi,sz,1);
x  = [50,50] + cat(2,r .*cos(a), r .*sin(a));
v  = cat(2,vr.*sin(a),-vr.*cos(a));
m  = unifrnd(ml,mh,sz,1); 

end