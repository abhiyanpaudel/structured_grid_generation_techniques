clc
clear all
close all

B=load('NASASC20714AirfoilCoordinatesSelig.txt');

x1=B(1:end,1);
y1=B(1:end,2);

%center airfoil at the origin
for i=1:length(x1)
x1(i) = x1(i) - 0.5;
end
R=5;
L=length(y1);

theta=linspace(0,2*pi,L);
x2=R*cos(theta);
y2=R*sin(theta);

 figure(1)
 plot(x1,y1,x2,y2)

x=zeros(L,L);
y=x;
z=y;

for i=1:L
x(i,1)=x1(i);
x(i,end)=x2(i);
y(i,1)=y1(i);
y(i,end)=y2(i);
end

%Create algebraic grid by linearly spaced elements between known x
%coordinates on airfoil and circle then linearly interpolate to find values
%of y
for i=1:L
deltaX=linspace(x(i,1),x(i,end),L);
for j=2:L-1
x(i,j)=deltaX(j);

y(i,j)=y(i,1)+(y(i,end)-y(i,1))*((x(i,j)-x(i,1))/(x(i,end)-x(i,1)));
end
end

figure(2)
surf(x,y,z)
view(2)

errX=1;
errY=1;
err = 0.001;

xold=x;
yold=y;

xi=linspace(0,1,L);
dxi=xi(2)-xi(1);

eta=linspace(1,abs(R/max(y1)),L);
deta=eta(2)-eta(1);

iter=0;

while errX > err || errY > err

for i=2:L-1
for j=2:L-1
alpha1= (1/(2*deta))*(xold(i,j+1)-x(i,j-1));
alpha2= (1/(2*deta))*(yold(i,j+1)-y(i,j-1));

gamma1 = (1/(2*dxi))*(xold(i+1,j)-x(i-1,j));
gamma2 = (1/(2*dxi))*(yold(i+1,j)-y(i-1,j));

Alpha = alpha1^2 + alpha2^2;
beta = alpha1*gamma1 + alpha2*gamma2;
gamma = gamma1^2 + gamma2^2;

factor = 1/(4*(Alpha*deta^2 + gamma*dxi^2));

x(i,j) = factor*(2*Alpha*deta^2*(xold(i+1,j)+x(i-1,j))-beta*dxi*deta*(xold(i+1,j+1)-xold(i+1,j-1)-xold(i-1,j+1)+x(i-1,j-1)) + 2*gamma*dxi^2*(xold(i,j+1)+x(i,j-1)));

y(i,j) = factor*(2*Alpha*deta^2*(yold(i+1,j)+y(i-1,j))-beta*dxi*deta*(yold(i+1,j+1)-yold(i+1,j-1)-yold(i-1,j+1) +y(i-1,j-1)) + 2*gamma*dxi^2*(yold(i,j+1)+y(i,j-1)));
end
end
iter=iter+1;
if mod(iter,10)==0
iter
errX = norm(xold-x,1)/norm(xold,1)
errY = norm(yold-y,1)/norm(yold,1)
end

xold=x;
yold=y;
end

figure(3)
surf(x,y,z)
view(2)

figure(4)
for i = 1:L
    for j = 1:L
        xxi(i,j) = xi(i);
        eeta(i,j) = eta(j);
    end
end
    
plot(xxi,eeta,'-r',xxi',eeta','-r')