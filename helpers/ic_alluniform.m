function [x,v,m] = ic_alluniform(sz,xl,xh,vl,vh,ml,mh)
% IC_ALLUNIFORM Generates uniformly distributed variables for
% positions/velocities/masses. The input arguments have l=low and h=high
% for setting boundaries. sz sets the sample size.

x = unifrnd(xl,xh,sz,2);
v = unifrnd(vl,vh,sz,2);
m = unifrnd(ml,mh,sz,1);

end