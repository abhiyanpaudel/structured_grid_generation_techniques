clc
clear all
close all

IM = 16;
JM = 12;
dxi = 1.0;
deta = 1.0;

for i = 1:IM
    for j = 1:JM
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
    end
end

R = 1.0;
a1 = 2.0;
b1 = 1.25;
a2 = 3.0;
b2 = 3.5;
beta = 5;
beta1 = (beta+1)/(beta-1);
theta = linspace(0,pi/2,IM);
for i = 1:IM
    RDIS(i,1) = 1/sqrt((cos(theta(i))^2/a1^2)+(sin(theta(i))^2/b1^2));
    rdis(i,JM) = 1/sqrt((cos(theta(i))^2/a2^2)+(sin(theta(i))^2/b2^2));
    delta(i) = rdis(i,JM)-RDIS(i,1);
end
for j = 1:JM
    for i = 1:IM
        gamma = (j-1)/(JM-1);
        nr = 1-(beta1)^(1-gamma);
        dr = 1+(beta1)^(1-gamma);
        beta2 = nr/dr;
        c(i,j) = delta(i)*(1+beta*beta2);
        DL(i,j) = RDIS(i,1)+c(i,j);
        x(i,j) = -DL(i,j) * cos(theta(i));
        y(i,j) = DL(i,j) * sin(theta(i));
    end
end


errormax = 0.01;
count = 1;
iteration = 0;
while count>0
    count = 0;
    prev_x = x;
    prev_y = y;
    
    % Discretized equation (5.16)
    % Form tridiagonal matrix
    
    for j = 2:JM-1
        p = 0;
        for i = 2:IM-1
            p = p+1;
%             x_eta = (prev_x(i,j+1)-prev_x(i,j-1))/(2*deta);
%             y_eta = (prev_y(i,j+1)-prev_y(i,j-1))/(2*deta);
%             x_xi = (prev_x(i+1,j)-prev_x(i-1,j))/(2*dxi);
%             y_xi = (prev_y(i+1,j)-prev_y(i-1,j))/(2*dxi);
%             a = x_eta^2+y_eta^2;
%             b = x_xi*x_eta+y_xi*y_eta;
%             c = x_xi^2+y_xi^2;
           a = ((prev_x(i,j+1)-prev_x(i,j-1))/(2*deta))^2+((prev_y(i,j+1)...
                -prev_y(i,j-1))/(2*deta))^2;
            b = ((prev_x(i+1,j)-prev_x(i-1,j))/(2*dxi))*((prev_x(i,j+1)...
                -prev_x(i,j-1))/(2*deta))+((prev_y(i+1,j)-prev_y(i-1,j))...
                /(2*dxi))*((prev_y(i,j+1)-prev_y(i,j-1))/(2*deta));
            c = ((prev_x(i+1,j)-prev_x(i-1,j))/(2*dxi))^2+((prev_y(i+1,j)-prev_y(i-1,j))...
                /(2*dxi))^2;
            a_1 = a/(dxi)^2;
            b_1 = b/(2*dxi*deta);
            c_1 = c/(deta)^2;
            d1(p) = a_1;
            d2(p) = -2*(a_1+c_1);
            d3(p) = a_1;
            
            % RHS of equation (5.16)
            if (i==2)
                Q1(p) = b_1*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1))...
                    -c_1*(x(i,j+1)+x(i,j-1))-a_1*x(1,j);
                Q2(p) = b_1*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1))...
                    -c_1*(y(i,j+1)+y(i,j-1))-a_1*y(1,j);
            elseif (i==IM-1)
                Q1(p) = b_1*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1))...
                    -c_1*(x(i,j+1)+x(i,j-1))-a_1*x(IM,j);
                Q2(p) = b_1*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1))...
                    -c_1*(y(i,j+1)+y(i,j-1))-a_1*y(IM,j);
            else
                Q1(p) = b_1*(x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1))...
                    -c_1*(x(i,j+1)+x(i,j-1));
                Q2(p) = b_1*(y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1))...
                    -c_1*(y(i,j+1)+y(i,j-1));
            end
            
        end
        P = diag(d1(1:p-1),-1)+diag(d2(1:p))+diag(d3(1:p-1),1);
        R1 = P\Q1';
        R2 = P\Q2';
        p = 0;
        for i = 2:IM-1
            p = p+1;
            x(i,j) = R1(p);
            y(i,j) = R2(p);
        end
        
    end
    errorX = 0 ;
    errorY = 0 ;
    for i = 2:IM-1
        for j = 2:JM-1
            errorX = errorX + abs(x(i,j)-prev_x(i,j));
            errorY = errorY + abs(y(i,j)-prev_y(i,j));
        end
    end
    errorT = errorX+errorY;
    count = count+1;
    iteration = iteration+1 , errorT
    
    if (errorT<errormax)
        break;
    end
end

% for i = 1:IM
%     for j = 1:JM
%         if i == 1
%         x_xi = (-3*x(i,j)+4*x(i+1,j)-x(i+2,j))/(2*dxi);
%         y_xi = (-3*y(i,j)+4*y(i+1,j)-y(i+2,j))/(2*dxi);
%
%         elseif i == IM
%         x_xi = (3*x(i,j)-4*x(i-1,j)+x(i-2,j))/(2*dxi);
%         y_xi = (3*y(i,j)-4*y(i-1,j)+x(i-2,j))/(2*dxi);
%
%         else
%         x_xi = (x(i+1,j)-x(i-1,j))/(2*dxi);
%         y_xi = (y(i+1,j)-y(i-1,j))/(2*dxi);
%         end
%
%
%         if j == 1
%         x_eta = (-3*x(i,j)+4*x(i,j+1)-x(i,j+2))/(2*deta);
%         y_eta = (-3*y(i,j)+4*y(i,j+1)-y(i,j+2))/(2*deta);
%
%         elseif j == JM
%         x_eta = (3*x(i,j)-4*x(i,j-1)+x(i,j-2))/(2*deta);
%         y_eta = (3*y(i,j)-4*y(i,j-1)+y(i,j-2))/(2*deta);
%
%         else
%         x_eta = (x(i,j+1)-x(i,j-1))/(2*deta);
%         y_eta = (y(i,j+1)-y(i,j-1))/(2*deta);
%         end
%         inv = (x_xi*y_eta)-(y_xi*x_eta);
%         J = 1/inv;
%         xi_x(i,j) = J*y_eta;
%         xi_y(i,j) = -J*x_eta;
%         eta_x(i,j) = -J*y_xi;
%         eta_y(i,j) = J*x_xi;
%
%     end
% end
% figure(1)
% plot(x,y,'-r',x',y','-r')
% xlabel('x')
% ylabel('y')

%
% figure(2)
% plot(xi,eta,'-r',xi',eta','-r')
% xlabel('\xi')
% ylabel('\eta')
%
% figure(3)
% plot3(xi,eta,xi_x,'-r',xi',eta',xi_x','-r')
% xlabel('\xi')
% ylabel('\eta')
% zlabel('\xi_{x}')
%
% figure(4)
% plot3(xi,eta,xi_y,'-r',xi',eta',xi_y','-r')
% xlabel('\xi')
% ylabel('\eta')
% zlabel('\xi_{y}')
%
% figure(5)
% plot3(xi,eta,eta_x,'-r',xi',eta',eta_x','-r')
% xlabel('\xi')
% ylabel('\eta')
% zlabel('\eta_{x}')
%
% figure(6)
% plot3(xi,eta,eta_y,'-r',xi',eta',eta_y','-r')
% xlabel('\xi')
% ylabel('\eta')
% zlabel('\eta_{y}')


toc






