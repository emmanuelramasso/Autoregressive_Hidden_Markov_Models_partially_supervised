function [x,pl]=signalstep(nu)
x=zeros(500,1);
x(126:250)=1;
x(251:375)=2;
x(376:500)=3;
x=x+0.1*rand(size(x));

pl=zeros(500,4)+nu;
pl(1:125,1)=1;
pl(126:250,2)=1;
pl(251:375,3)=1;
pl(376:500,4)=1;