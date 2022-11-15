tic
clc
clear all
close all

c = 1.0;
IM = 39;
JM = 21;
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


dxi = 1;
deta = dxi;

for i = 1:IM
deltax = linspace(x(i,1),x(i,JM),JM);
    for j = 2:JM-1
        x(i,j) = deltax(j);
        y(i,j) = y(i,1)+((y(i,JM)-y(i,1))/(x(i,JM)-x(i,1)))*(x(i,j)-x(i,1));
    end
end

figure(1)
plot(x,y,'-r',x',y','-r')
errormax = 0.01;
count = 1;
iteration = 0;


while count>0
    count = 0;
    prev_x = x;
    prev_y = y;
    
%     Update x(1,j),y(1,j) and x(IM,j),y(IM,j)
for i = 1
    for j = 2:JM-1
        a = ((x(i,j+1)-x(i,j-1))/(2*deta))^2+((y(i,j+1)-y(i,j-1))/(2*deta))^2;
        b = ((x(i+1,j)-x(IM-1,j))/(2*dxi))*((x(i,j+1)-x(i,j-1))/(2*deta))+((y(i+1,j)-y(IM-1,j))/(2*dxi))*((y(i,j+1)-y(i,j-1))/(2*deta));
        c = ((x(i+1,j)-x(IM-1,j))/(2*deta))^2+((y(i+1,j)-y(IM-1,j))/(2*deta))^2;
        a_1 = a/(dxi)^2;
        b_1 = b/(2*dxi*deta);
        c_1 = c/(deta)^2;
        ax_2 = a_1*(x(i+1,j)+x(IM-1,j));
        bx_2 = b_1*(x(i+1,j+1)-x(i+1,j-1)+x(IM-1,j-1)-x(IM-1,j+1));
        cx_2 = c_1*(x(i,j+1)+x(i,j-1));
        ay_2 = a_1*(y(i+1,j)+y(IM-1,j));
        by_2 = b_1*(y(i+1,j+1)-y(i+1,j-1)+y(IM-1,j-1)-y(IM-1,j+1));
        cy_2 = c_1*(y(i,j+1)+y(i,j-1));
        x(i,j) = (ax_2+cx_2-bx_2)/(2*(a_1+c_1));
        y(i,j) = (ay_2+cy_2-by_2)/(2*(a_1+c_1));
    end
end
x(IM,:) = x(1,:);
y(IM,:) = y(1,:);
    for i = 2:IM-1
        for j = 2:JM-1
            a = ((x(i,j+1)-x(i,j-1))/(2*deta))^2+((y(i,j+1)-y(i,j-1))/(2*deta))^2;
            b = ((x(i+1,j)-x(i-1,j))/(2*dxi))*((x(i,j+1)-x(i,j-1))/(2*deta))+((y(i+1,j)-y(i-1,j))/(2*dxi))*((y(i,j+1)-y(i,j-1))/(2*deta));
            c = ((x(i+1,j)-x(i-1,j))/(2*deta))^2+((y(i+1,j)-y(i-1,j))/(2*deta))^2;
            a_1 = a/(dxi)^2;
            b_1 = b/(2*dxi*deta);
            c_1 = c/(deta)^2;
            ax_2 = a_1*(x(i+1,j)+x(i-1,j));
            bx_2 = b_1*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1));
            cx_2 = c_1*(x(i,j+1)+x(i,j-1));
            ay_2 = a_1*(y(i+1,j)+y(i-1,j));
            by_2 = b_1*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1));
            cy_2 = c_1*(y(i,j+1)+y(i,j-1));
            x(i,j) = (ax_2+cx_2-bx_2)/(2*(a_1+c_1));
            y(i,j) = (ay_2+cy_2-by_2)/(2*(a_1+c_1));
        end
    end
   
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

for i = 1:IM
    for j = 1:JM
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
    end
end
        
figure(2)
plot(x,y,'-r',x',y','-r')
xlabel('x')
ylabel('y')

figure(3)
plot(xi,eta,'-r',xi',eta','-r')
toc




