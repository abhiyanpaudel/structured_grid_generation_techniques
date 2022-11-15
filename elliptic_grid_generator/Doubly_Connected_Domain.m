tic
clc
clear all
close all

c = 1.0;
IM = 59;
JM = 41;
R = 2;



%Airfoil coordinates
nacaseries = input('Enter the 4-digit naca series = ','s');


 % creating points on airfoil

 s1 = str2double(nacaseries(1));
 s2 = str2double(nacaseries(2));
 s3 = str2double(nacaseries(3));
 s4 = str2double(nacaseries(4));
  m = s1*0.01; p = s2*0.1 ; t = (10*s3+s4)*0.01;
 n = IM-1;

for i= n:-1:1
    
    theta = (i-n)*2*pi/n;                 
     x(n+1-i,JM) = R*cos(theta);
     y(n+1-i,JM) = R*sin(theta);
    xc = 0.5*c*(1+cos(theta));
if(xc/c)<p
    yc(n+1-i) = m*c/p^2*(2*p*(xc/c)-(xc/c)^2);
    dydx(n+1-i) = (2*m/p^2)* (p-xc/c);
    beta(n+1-i) = atan(dydx(n+1-i));
else
    yc(n+1-i) = m*c/(1-p)^2 * ((1-2*p)+2*p*(xc/c)-(xc/c)^2);
    dydx(n+1-i) = (2*m/(1-p)^2)* (p-xc/c);
    beta(n+1-i) = atan(dydx(n+1-i));
end
yt=5*t*c*(0.2969*sqrt(xc/c)-0.1260*(xc/c)...
    -0.3516*(xc/c)^2+0.2843*(xc/c)^3-0.1036*(xc/c)^4);

if(i<(0.5*n+1))
    xa(n+1-i)=xc - yt*sin(beta(n+1-i));
    ya(n+1-i)=yc(n+1-i)+yt*cos(beta(n+1-i));
else
    xa(n+1-i)=xc + yt*sin(beta(n+1-i));
    ya(n+1-i)=yc(n+1-i)-yt*cos(beta(n+1-i));
end

end
x(n+1,JM) = R;
y(n+1,JM) = 0;
xa(n+1)= c ; 
ya(n+1) = 0; 
yc(n+1) = 0;  % trailing edge

%Place midchord of the airfoil at the origin
for i = 1:IM
    x(i,1) = xa(i)-0.5;
    y(i,1) = ya(i);
end

dxi = 1.0;
deta = 1.0;

% computing ccoordinates of computational domain
for i = 1:IM
    for j = 1:JM
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
    end
end

beta = 1.02;
beta1 = (beta+1)/(beta-1);


% Algebraic grid system generator 

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
plot(x,y,'-r',x',y','-r')
axis equal 
errormax = 0.01;
count = 1;
iteration = 0;

% Elliptic grid generator 
while count>0
    count = 0;
    prev_x = x;
    prev_y = y;
    
%     Update x(1,j),y(1,j) and x(IM,j),y(IM,j)
%     computation of coefficients a,b & c lag by one iterative level

 for i = 1
        for j = 2:JM-1
           
    a = ((x(i,j+1)-x(i,j-1))/(2*deta))^2+((y(i,j+1)-y(i,j-1))/(2*deta))^2;
    b = ((x(i+1,j)-x(IM-1,j))/(2*dxi))*((x(i,j+1)-x(i,j-1))/(2*deta))+((y(i+1,j)-y(IM-1,j))/(2*dxi))*((y(i,j+1)-y(i,j-1))/(2*deta));
    c = ((x(i+1,j)-x(IM-1,j))/(2*deta))^2+((y(i+1,j)-y(IM-1,j))/(2*deta))^2;
     a_1(i,j) = a/(dxi)^2;
     b_1(i,j) = b/(2*dxi*deta);
     c_1(i,j) = c/(deta)^2;
        end
    end
for i = 1
    for j = 2:JM-1
        ax_2 = a_1(i,j)*(x(i+1,j)+x(IM-1,j));
        bx_2 = b_1(i,j)*(x(i+1,j+1)-x(i+1,j-1)+x(IM-1,j-1)-x(IM-1,j+1));
        cx_2 = c_1(i,j)*(x(i,j+1)+x(i,j-1));
        ay_2 = a_1(i,j)*(y(i+1,j)+y(IM-1,j));
        by_2 = b_1(i,j)*(y(i+1,j+1)-y(i+1,j-1)+y(IM-1,j-1)-y(IM-1,j+1));
        cy_2 = c_1(i,j)*(y(i,j+1)+y(i,j-1));
        x(i,j) = (ax_2+cx_2-bx_2)/(2*(a_1(i,j)+c_1(i,j)));
        y(i,j) = (ay_2+cy_2-by_2)/(2*(a_1(i,j)+c_1(i,j)));
    end
end
x(IM,:) = x(1,:);
y(IM,:) = y(1,:);
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


% computation of metrics 
for i = 1:IM
    for j = 1:JM
        if i == 1 
        x_xi(i,j) = (-3*x(i,j)+4*x(i+1,j)-x(i+2,j))/(2*dxi);
        y_xi(i,j) = (-3*y(i,j)+4*y(i+1,j)-y(i+2,j))/(2*dxi);
        
        elseif i == IM 
        x_xi(i,j) = (3*x(i,j)-4*x(i-1,j)+x(i-2,j))/(2*dxi);
        y_xi(i,j) = (3*y(i,j)-4*y(i-1,j)+x(i-2,j))/(2*dxi);
        
        else 
        x_xi(i,j) = (x(i+1,j)-x(i-1,j))/(2*dxi);
        y_xi(i,j) = (y(i+1,j)-y(i-1,j))/(2*dxi);
        end
            
        
        if j == 1
        x_eta(i,j) = (-3*x(i,j)+4*x(i,j+1)-x(i,j+2))/(2*deta);
        y_eta(i,j) = (-3*y(i,j)+4*y(i,j+1)-y(i,j+2))/(2*deta);
        
        elseif j == JM
        x_eta(i,j) = (3*x(i,j)-4*x(i,j-1)+x(i,j-2))/(2*deta);
        y_eta(i,j) = (3*y(i,j)-4*y(i,j-1)+y(i,j-2))/(2*deta);

        else
        x_eta(i,j) = (x(i,j+1)-x(i,j-1))/(2*deta);
        y_eta(i,j) = (y(i,j+1)-y(i,j-1))/(2*deta);
        end 

    end
end


        
figure(2)
plot(x,y,'-r',x',y','-r')
xlabel('x')
ylabel('y')


figure(3)
plot(xi,eta,'-r',xi',eta','-r')
xlabel('\xi')
ylabel('\eta')

figure(4)
plot3(xi,eta,x_xi,'-r',xi',eta',x_xi','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('x_{\xi}')
view(139,26)

figure(5)
plot3(xi,eta,y_xi,'-r',xi',eta',y_xi','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('y_{\xi}')
view(139,26)

figure(6)
plot3(xi,eta,x_eta,'-r',xi',eta',x_eta','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('x_{\eta}')
view(139,26)

figure(7)
plot3(xi,eta,y_eta,'-r',xi',eta',y_eta','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('y_{\eta}')
view(139,26)
toc

% Author comment : y_xi has to be fixed...Others are perfect..



