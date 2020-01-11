%%%% DAUBECHIES 5/3 TAP FILTER FOR INVERSE TRANSFORM WITHOUT ELIMINATION %%%%

function [x]=idbmw(a,h,v,d)
[f,g]=size(a);

% MERGING a&h and v&d bands
for j=1:g
for i=1:f
p(2*i,j)=a(i,j);
p(2*i-1,j)=h(i,j);
q(2*i,j)=v(i,j);
q(2*i-1,j)=d(i,j);
i=i+1;
end
j=j+1;
end

% INVERSE COLUMN TRANSFORM
for j=1:g
for i=1:f
if i==1
m(2*i-1,j)=(6*p(2*i-1,j)+2*p(2*i,j)-p(2*i+1,j)+4)/8;
m(2*i,j)=(-p(2*i-1,j)+2*p(2*i,j)-p(2*i+1,j))/4;
n(2*i-1,j)=(6*q(2*i-1,j)+2*q(2*i,j)-q(2*i+1,j)+4)/8;
n(2*i,j)=(-q(2*i-1,j)+2*q(2*i,j)-q(2*i+1,j))/4;
elseif i==f
m(2*i-1,j)=(-p(2*i-3,j)+2*p(2*i-2,j)+6*p(2*i-1,j)+2*p(2*i,j)+4)/8;
m(2*i,j)=(-p(2*i-1,j)+2*p(2*i,j))/4;
n(2*i-1,j)=(-q(2*i-3,j)+2*q(2*i-2,j)+6*q(2*i-1,j)+2*q(2*i,j)+4)/8;
n(2*i,j)=(-q(2*i-1,j)+2*q(2*i,j))/4;
else 
m(2*i-1,j)=(-p(2*i-3,j)+2*p(2*i-2,j)+6*p(2*i-1,j)+2*p(2*i,j)-p(2*i+1,j)+4)/8;
m(2*i,j)=(-p(2*i-1,j)+2*p(2*i,j)-p(2*i+1,j))/4;
n(2*i-1,j)=(-q(2*i-3,j)+2*q(2*i-2,j)+6*q(2*i-1,j)+2*q(2*i,j)-q(2*i+1,j)+4)/8;
n(2*i,j)=(-q(2*i-1,j)+2*q(2*i,j)-q(2*i+1,j))/4;
end
i=i+1;
end
j=j+1;
end
m=floor(m);
n=floor(n);

% MERGING m&n bands
for i=1:2*f
for j=1:g
w(i,2*j)=m(i,j);
w(i,2*j-1)=n(i,j);
j=j+1;
end
i=i+1;
end

% INVERSE ROW TRANSFORM
for i=1:2*f
for j=1:g
if j==1
x(i,2*j-1)=(6*w(i,2*j-1)+2*w(i,2*j)-w(i,2*j+1)+4)/8;
x(i,2*j)=(-w(i,2*j-1)+2*w(i,2*j)-w(i,2*j+1))/4;
elseif j==g
x(i,2*j-1)=(-w(i,2*j-3)+2*w(i,2*j-2)+6*w(i,2*j-1)+2*w(i,2*j)+4)/8;
x(i,2*j)=(-w(i,2*j-1)+2*w(i,2*j))/4;
else 
x(i,2*j-1)=(-w(i,2*j-3)+2*w(i,2*j-2)+6*w(i,2*j-1)+2*w(i,2*j)-w(i,2*j+1)+4)/8;
x(i,2*j)=(-w(i,2*j-1)+2*w(i,2*j)-w(i,2*j+1))/4;
end
j=j+1;
end
i=i+1;
end
x=floor(x);