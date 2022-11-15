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
beta = 1.05;                                    %clustering parameter 
alpha = 0.5;                                    
beta1 = (beta+1)/(beta-1);
bet_alph = 2*alpha+beta;
for i = 1:IM 
    for j = 1:JM 
   %compute coordinates of computational domain 
        xi(i,j) = (i-1)*dxi;
        eta(i,j) = (j-1)*deta;
  %computing the coordinates of physical domain acoording to the
  %relation given by (9.41) and (9.42)       
        expo = (eta(i,j)-alpha)/(1-alpha);
        beta2 = beta1^expo;
        x(i,j) = xi(i,j);
        y(i,j) = H*(bet_alph*beta2+2*alpha-beta)/((2*alpha+1)*(1+beta2));
   % computing metrics using equations(9.46)      
        y_term  = ((2*alpha+1)*(y(i,j)/H)-2*alpha)^2;
        eta_y(i,j) = (2*beta*(2*alpha+1)*(1-alpha))/(H*(beta*beta-y_term)*log(beta1));
        
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
view(-66,54)

 
  