clc
clear all
close all

l_xi = 1.0;
l_eta = 1.0;
H = 0.8;
IM = 21;
JM = 24;
dxi = l_xi/(IM-1);
deta = l_eta/(JM-1);
beta = 0.4;                            % clustering parameter                                     
beta1 = (beta+1)/(beta-1);
for i = 1:IM 
    for j = 1:JM 
%         compute coordinates of computational domain 
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
        beta2 = beta1^(1-eta(i,j));
%         computing the coordinates of physical domain acoording to the
%         relation given by (9.33) and (9.34)     
        x(i,j) = xi(i,j);
        y(i,j) = H*((beta+1)-((beta-1)*beta2))/(beta2+1);
        y_term  = (1-(y(i,j)/H))^2;
        
 %        computing metrics using equations(9.38) 
        eta_y(i,j) = 2*beta/(H*(beta*beta-y_term)*log(beta1));
        
    end
end

% Plotting computational domain 
figure(1)
  plot(xi,eta,'-r',xi',eta','-r')  
  axis image
  xlabel('\xi')
ylabel('\eta')
  
% Plotting physical domain 
  figure(2)
  plot(x,y,'-r',x',y','-r')
  axis image
  xlabel('x')
ylabel('y')
  
% Plotting distribution of metrics 
  figure(3)
  plot3(xi,eta,eta_y,'-r',xi',eta',eta_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{y}')
view(111,60)
 
  