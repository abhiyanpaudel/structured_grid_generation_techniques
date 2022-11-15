clc
clear all
close all

l_xi = 1.0;
l_eta = 1.0;
H = 0.8; 
D = H/4;                       % y-coordinate where clustering is desired           
IM = 21;
JM = 24;
dxi = l_xi/(IM-1);
deta = l_eta/(JM-1);
beta =5;                                   % clustering parameter 
D_H = D/H;
a = (1+(exp(beta)-1)*D_H)/(1+(exp(-beta)-1)*D_H)
A = 0.5/beta*log(a)

for i = 1:IM 
    for j = 1:JM 
  %  compute coordinates of computational domain 
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
  %compute the coordinates of physical domain acoording to the
  %relation given by (9.49) and (9.50)    
        x(i,j) = xi(i,j);
        y(i,j) = D*(1+sinh(beta*(eta(i,j)-A))/sinh(beta*A));
   % computing metrics using equations(9.54) 
        y_term  = 1+(((y(i,j)/D)-1))^2*(sinh(beta*A))^2;
        eta_y(i,j) = sinh(beta*A)/(beta*D*y_term^0.5);
        
    end
end


% Plot computational domain 
figure(1)
  plot(xi,eta,'-r',xi',eta','-r')  
  axis image
  
 
% Plot physical domain 
  figure(2)
  plot(x,y,'-r',x',y','-r')
axis([0 1 0 0.9])
  
% Plot distribution of metrics
  figure(3)
  plot3(xi,eta,eta_y,'-r',xi',eta',eta_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{y}')
view(-70,40)

 
  