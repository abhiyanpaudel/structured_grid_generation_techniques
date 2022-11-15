clc
clear all
close all

R = 1.0;
a1 = 2.0;
b1 = 1.8;
b2 = 3.0;
beta = 1.4;
beta1 = (beta+1)/(beta-1);

l_xi = 1.0;
l_eta = 1.0;
KM = 37;
JM = 26;
dxi = l_xi/(KM-1);
deta = l_eta/(JM-1);
delthet = pi/(KM-1);

for k = 1:KM
    for j = 1:JM
%         compute coordinates of computational domain 
        xi(k,j) = (k-1)*dxi;
        eta(k,j) = (JM-j)*deta;
   
 %         compute coordinates of physical domain using equations (9.55),(9.56),(9.57) ,(9.58),(9.59) and (9.60)   
        theta = (k-1)*delthet;
        beta2 = beta1^(eta(k,j));
        if j == 1
            y(k,1) = -R*cos(theta);
            z(k,1) = R*sin(theta);
        end
        if k<=KM/2
            r = 1/sqrt((sin(theta)^2/a1^2)+(cos(theta)^2/b1^2));
        else
            r = 1/sqrt((sin(theta)^2/a1^2)+(cos(theta)^2/b2^2));
        end
        delta = r-R;
        c(k,j) = delta*(1-((beta*(beta2-1))/(beta2+1)));
        y(k,j) = y(k,1)-c(k,j)*cos(theta);
        z(k,j) = z(k,1)+c(k,j)*sin(theta);
    end
end

% compute metric distribution numerically 
for k = 1:KM
    for j = 1:JM
        if k == 1 
        z_xi = (-3*z(k,j)+4*z(k+1,j)-z(k+2,j))/(2*dxi);
        y_xi = (-3*y(k,j)+4*y(k+1,j)-y(k+2,j))/(2*dxi);
        
        elseif k == KM 
        z_xi = (3*z(k,j)-4*z(k-1,j)+z(k-2,j))/(2*dxi);
        y_xi = (3*y(k,j)-4*y(k-1,j)+z(k-2,j))/(2*dxi);
        
        else 
        z_xi = (z(k+1,j)-z(k-1,j))/(2*dxi);
        y_xi = (y(k+1,j)-y(k-1,j))/(2*dxi);
        end
            
        
        if j == 1
        z_eta = (-3*z(k,j)+4*z(k,j+1)-z(k,j+2))/(2*deta);
        y_eta = (-3*y(k,j)+4*y(k,j+1)-y(k,j+2))/(2*deta);
        
        elseif j == JM
        z_eta = (3*z(k,j)-4*z(k,j-1)+z(k,j-2))/(2*deta);
        y_eta = (3*y(k,j)-4*y(k,j-1)+y(k,j-2))/(2*deta);

        else
        z_eta = (z(k,j+1)-z(k,j-1))/(2*deta);
        y_eta = (y(k,j+1)-y(k,j-1))/(2*deta);
        end 
        inv = (z_xi*y_eta)-(y_xi*z_eta);
        J = 1/inv;
        xi_z(k,j) = J*y_eta;
        xi_y(k,j) = -J*z_eta;
        eta_z(k,j) = -J*y_xi;
        eta_y(k,j) = J*z_xi;
        
    end
end

% plot physical domain
figure(1)
plot(z,y,'-r',z',y','-r')
xlabel('z')
ylabel('y')

% plot computational domain 
figure(2)
plot(xi,eta,'-r',xi',eta','-r')
xlabel('\xi')
ylabel('\eta')

% plot distribution of metrics 
figure(3)
plot3(xi,eta,xi_z,'-r',xi',eta',xi_z','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\xi_{z}')
view(-42,40)

figure(4)
plot3(xi,eta,xi_y,'-r',xi',eta',xi_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\xi_{y}')
view(-42,40)

figure(5)
plot3(xi,eta,eta_z,'-r',xi',eta',eta_z','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{z}')
view(-42,40)

figure(6)
plot3(xi,eta,eta_y,'-r',xi',eta',eta_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{y}')
view(-42,40)
