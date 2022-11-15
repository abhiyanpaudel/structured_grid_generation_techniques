clc
clear all
close all


R = 1.0;
H = 6.0;
L = 6.0;
IM = 71;
JM = 71;


theta = linspace(pi/4,(2*pi+pi/4),IM);
aa = 0.0;
bb = 0.0;
for i = 1:IM
    if theta(i)>=pi/4 && theta(i)<=3*pi/4
    aa = aa+1;
    elseif theta(i)>3*pi/4 && theta(i)<5*pi/4
    bb = bb+1;
    end
    x(i,1) = 0.5*L+R*cos(theta(i));
    y(i,1) = 0.5*H+R*sin(theta(i));
end
dxx = L/(aa-1);
dyy = H/(bb+1);

% Discretize boundary of rectangular domain 

for i = 1:IM
    if i<=aa
        x(i,JM) = (aa-i)*dxx;
        y(i,JM) = H;
    elseif i>aa && i<=aa+bb
        x(i,JM) = 0.0;
        y(i,JM) = (aa+bb+1-i)*dyy;
    elseif i>aa+bb && i<=aa+bb+aa
         x(i,JM) = dxx*(i-(aa+bb+1));
        y(i,JM) =  0.0;
    else
        x(i,JM) = L;
        y(i,JM) = dyy*(i-(aa+bb+aa));
    end
end


dxi = 1.0;
deta = 1.0;
for i = 1:IM
    for j = 1:JM
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
    end
end

beta = 1.2;
beta1 = (beta+1)/(beta-1);

for i = 1:IM
    dx = x(i,JM)-x(i,1);
    dy = y(i,JM)-y(i,1);
    delta = sqrt(dx^2+dy^2);
    phi = atan2(dy,dx); 
    for j = 1:JM
         gamma = (j-1)/(JM-1);
        num = 1-(beta1)^(1-gamma);
        den = 1+(beta1)^(1-gamma);
        beta2 = num/den;
        c = delta*(1+beta*beta2);
        x(i,j) = x(i,1)+c * cos(phi);
        y(i,j) = y(i,1)+c * sin(phi);
    end
end

figure(1)
plot(x(:,1),y(:,1),'-r',x(:,JM),y(:,JM),'-r')
axis([0 L 0 H])

figure(2)
plot(x,y,'-r',x',y','-r')
axis([0 L 0 H])
