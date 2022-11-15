clc
clear all
close all

A1 = 3.0;
A2 = 2.0;
D = 2.0;
R = D/2;
H = 4.0;
I1 = 10;
IMID = 14;
I2 = 18;
IM = 28;
JM = 30;

% Place the half cylinder at the distance of A1 from inlet 

theta = linspace(pi,0,IM);

for i = 1:IM
    x(i,1) = A1+R+R*cos(theta(i));
    y(i,1) = R*sin(theta(i));
end

dy1 = H/(I1-1);
dy2 = H/(IM-I2);
dx = (A1+D+A2)/(I2-I1);

for i = 1:IM
    if i<=I1 
        x(i,JM) = 0.0;
        y(i,JM) = (i-1)*dy1;
    elseif i>I1 && i<=I2
        x(i,JM) = (i-I1)*dx;
        y(i,JM) = H;
    else
        x(i,JM) = A1+D+A2;
        y(i,JM) = (IM-i)*dy2;
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

beta = 1.07;
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
axis([0 A1+A2+D 0 H])

figure(2)
plot(x,y,'-r',x',y','-r')
axis([0 A1+A2+D 0 H])

errormax = 0.00001;
count = 1;
iteration = 0;

% compute grid by elliptic PDE scheme
while count>0
    count = 0;
    prev_x = x;
    prev_y = y;

%    Equation (9.68) & (9.69)    
%     computation of coefficients a,b & c lag by one iterative level
    for i = 2:IM-1
        for j = 2:JM-1
           
     a = ((x(i,j+1)-x(i,j-1))/(2*deta))^2+((y(i,j+1)-y(i,j-1))/(2*deta))^2;
     b = ((x(i+1,j)-x(i-1,j))/(2*dxi))*((x(i,j+1)-x(i,j-1))/(2*deta))+((y(i+1,j)-y(i-1,j))/(2*dxi))*((y(i,j+1)-y(i,j-1))/(2*deta));
     c = ((x(i+1,j)-x(i-1,j))/(2*dxi))^2+((y(i+1,j)-y(i-1,j))/(2*dxi))^2;
     a_1(i,j) = a/(dxi)^2;
     b_1(i,j) = b/(2*dxi*deta);
     c_1(i,j) = c/(deta)^2;
        end
    end
    
    

    for i = 2:IM-1
        for j = 2:JM-1
            ax_2 = a_1(i,j)*(x(i+1,j)+x(i-1,j));
            bx_2 = b_1(i,j)*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1));
            cx_2 = c_1(i,j)*(x(i,j+1)+x(i,j-1));
            ay_2 = a_1(i,j)*(y(i+1,j)+y(i-1,j));
            by_2 = b_1(i,j)*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1));
            cy_2 = c_1(i,j)*(y(i,j+1)+y(i,j-1));   
            x(i,j) = (ax_2+cx_2-bx_2)/(2*(a_1(i,j)+c_1(i,j)));
            y(i,j) = (ay_2+cy_2-by_2)/(2*(a_1(i,j)+c_1(i,j)));
        end
    end
   
%     Error monitoring 
    errorX = 0 ;
    errorY = 0 ;
    

    for i = 1:IM
        for j = 1:JM
            errorX = errorX + abs(x(i,j)-prev_x(i,j));
            errorY = errorY + abs(y(i,j)-prev_y(i,j));
            
        end
    end
    errorT = errorX+errorY;
    count = count+1;
    iteration = iteration+1,  errorT
    
    if (errorT<errormax)
        break;
    end
end

figure(3)
plot(x,y,'-r',x',y','-r')
axis([0 A1+A2+D 0 H])

for i = 1:IM
    for j = 1:JM
        if i == 1 
        x_xi = (-3*x(i,j)+4*x(i+1,j)-x(i+2,j))/(2*dxi);
        y_xi = (-3*y(i,j)+4*y(i+1,j)-y(i+2,j))/(2*dxi);
        
        elseif i == IM 
        x_xi = (3*x(i,j)-4*x(i-1,j)+x(i-2,j))/(2*dxi);
        y_xi = (3*y(i,j)-4*y(i-1,j)+x(i-2,j))/(2*dxi);
        
        else 
        x_xi = (x(i+1,j)-x(i-1,j))/(2*dxi);
        y_xi = (y(i+1,j)-y(i-1,j))/(2*dxi);
        end
            
        
        if j == 1
        x_eta = (-3*x(i,j)+4*x(i,j+1)-x(i,j+2))/(2*deta);
        y_eta = (-3*y(i,j)+4*y(i,j+1)-y(i,j+2))/(2*deta);
        
        elseif j == JM
        x_eta = (3*x(i,j)-4*x(i,j-1)+x(i,j-2))/(2*deta);
        y_eta = (3*y(i,j)-4*y(i,j-1)+y(i,j-2))/(2*deta);

        else
        x_eta = (x(i,j+1)-x(i,j-1))/(2*deta);
        y_eta = (y(i,j+1)-y(i,j-1))/(2*deta);
        end 
        inv = (x_xi*y_eta)-(y_xi*x_eta);
        J = 1/inv;
        xi_x(i,j) = J*y_eta;
        xi_y(i,j) = -J*x_eta;
        eta_x(i,j) = -J*y_xi;
        eta_y(i,j) = J*x_xi;
        
    end
end


figure(4)
plot(xi,eta,'-r',xi',eta','-r')
xlabel('\xi')
ylabel('\eta')

figure(5)
plot3(xi,eta,xi_x,'-r',xi',eta',xi_x','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\xi_{x}')
view(135,48)
figure(6)
plot3(xi,eta,xi_y,'-r',xi',eta',xi_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\xi_{y}')
view(135,48)
figure(7)
plot3(xi,eta,eta_x,'-r',xi',eta',eta_x','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{x}')
view(135,48)
figure(8)
plot3(xi,eta,eta_y,'-r',xi',eta',eta_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{y}')
view(135,48)
toc

