tic
clc
clear all
close all

c = 1.0;
IT1 = 30;
IL = 48;
IM = 2*IL-1;
JM = 50;
R = c;
L = 1.5;
n = IM-2*(IT1-1)-1;                    % Number of nodes on airfoil
%Airfoil coordinates
nacaseries = input('Enter the 4-digit naca series = ','s');


 % creating points on airfoil

 s1 = str2double(nacaseries(1));
 s2 = str2double(nacaseries(2));
 s3 = str2double(nacaseries(3));
 s4 = str2double(nacaseries(4));
  m = s1*0.01; p = s2*0.1 ; t = (10*s3+s4)*0.01;


for i= 1:n
    
    theta = (i-1)*2*pi/n;
    xc = 0.5*c*(1+cos(theta));
if(xc/c)<p
    yc(i) = m*c/p^2*(2*p*(xc/c)-(xc/c)^2);
    dydx(i) = (2*m/p^2)* (p-xc/c);
    beta(i) = atan(dydx(i));
else
    yc(i) = m*c/(1-p)^2 * ((1-2*p)+2*p*(xc/c)-(xc/c)^2);
    dydx(i) = (2*m/(1-p)^2)* (p-xc/c);
    beta(i) = atan(dydx(i));
end
yt=5*t*c*(0.2969*sqrt(xc/c)-0.1260*(xc/c)...
    -0.3516*(xc/c)^2+0.2843*(xc/c)^3-0.1036*(xc/c)^4);

if(i<(0.5*n+1))
    xa(i)=xc - yt*sin(beta(i));
    ya(i)=yc(i)+yt*cos(beta(i));
else
    xa(i)=xc + yt*sin(beta(i));
    ya(i)=yc(i)-yt*cos(beta(i));
end

end
xa(n+1)= xa(1) ; 
ya(n+1) = ya(1); 
yc(n+1) = yc(1);  % trailing edge

figure(1)
plot(xa,ya,'*-r')
axis equal
% compute x(i,1) starting from rear of physical domain
dx = L/(IT1-1);
aa = 0;
bb = 0;
cc = 0;
for i = 1:IM
    if i<IT1
        aa = aa+1;
        x(i,1) = xa(n+1)+dx*(IT1-i);
        y(i,1) = ya(n+1);
    elseif i>=IT1 && i<=IT1+n
        bb = bb+1;
        x(i,1) = xa(i+1-IT1);
        y(i,1) = ya(i+1-IT1);
    else
        cc = cc+1;
        x(i,1) = xa(n+1)+dx*(i-IT1-n);
        y(i,1) = ya(n+1);
    end
    
end

% Discretize outer boundary

I2 = 9;
c_m = IM-2*(I2+IT1-1);
phi = linspace(pi/2,3*pi/2,c_m);

dx1 = c/I2; 
sc_n = IT1+I2-1+c_m;
for i = 1:IM
    if i<IT1
        x(i,JM) = x(i,1);
        y(i,JM) = y(i,1)+R;
    elseif i>=IT1 && i<=(IT1+I2-1)
        x(i,JM) =  dx1*(IT1+I2-i);
        y(i,JM) = ya(n+1)+R;
    elseif i>(IT1+I2-1) && i<=sc_n
        x(i,JM) = R*cos(phi(i-(IT1+I2-1)));
        y(i,JM) = R*sin(phi(i-(IT1+I2-1)));
    elseif i>sc_n && i<=(sc_n+I2)
        x(i,JM) = dx1*(i-sc_n);
        y(i,JM) = ya(n+1)-R;
    else
        x(i,JM) = c + dx*(i-(sc_n+I2));
        y(i,JM) = ya(n+1)-R;
    end
end
        
    figure(2)
    plot(x(:,1),y(:,1),x(:,JM),y(:,JM),'-r')
    
        dxi = 1.0;
        deta = 1.0;
        for i = 1:IM
            for j = 1:JM
                xi(i,j) = (i-1)*dxi;
                eta(i,j) = (j-1)*deta;
            end
        end
        
        beta = 1.1;
        beta1 = (beta+1)/(beta-1);
        
        for i = 1:IM
            dx = x(i,JM)-x(i,1);
            dy = y(i,JM)-y(i,1);
            delta = sqrt(dx^2+dy^2);
            alpha = atan2(dy,dx);
            for j = 1:JM
                gamma = (j-1)/(JM-1);
                num = 1-(beta1)^(1-gamma);
                den = 1+(beta1)^(1-gamma);
                beta2 = num/den;
                cij = delta*(1+beta*beta2);
                x(i,j) = x(i,1)+cij * cos(alpha);
                y(i,j) = y(i,1)+cij * sin(alpha);
            end
        end
        
        figure(3)
        plot(x,y,'-r',x',y','-r')


        
        
        
