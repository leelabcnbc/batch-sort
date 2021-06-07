function max_M=sac_distribution_inverse(nu,len,P);

x=[0.1:.2:len*5];
y=cumsum(betad(1./(1+x/nu),nu/2+2,len/2));
y=y/y(end);
max_M=spline(y,x,P);