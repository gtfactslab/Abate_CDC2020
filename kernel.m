function kern = kernel(a,b,param)
[A,B]=meshgrid(a,b);
sqdist = (A-B).^2;
kern=exp(-1/2*(1/param)*sqdist);