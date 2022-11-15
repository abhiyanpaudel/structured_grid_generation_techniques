clc
clear all
close all

L = 4;
H1 = 2;
H2 = 4;
IM = 17;
JM = 13;
dgsi = L/(IM-1);
deta = 1/(JM-1);
H = (H2-H1)/L;

for i = 1:IM
    for j = 1:JM 
%         computing the coordinates of computational domain 
        gsi(i,j) = (i-1)*dgsi;
        eta(i,j) = (j-1)*deta;
        
%         computing the coordinates of physical domain acoording to the
%         relation given by (9.23) and (9.24) 
        x(i,j) = gsi(i,j);
        y(i,j) = (H1+H*gsi(i,j))*eta(i,j);
        
%        computing metrics using equations(9.27),(9.28),(9.29) and (9.30) 
        gsi_x(i,j) = 1;
        gsi_y(i,j) = 0;
        eta_x(i,j) = -H*eta(i,j)/(H1+H*gsi(i,j));
        eta_y(i,j) = 1/(H1+H*gsi(i,j));
     
    end
end 

% plotting coordinates of computational domain
figure(1)
 plot(gsi,eta,'-r',gsi',eta','-r')
 xlabel('\xi')
 ylabel('\eta')

 
% plotting coordinates of physical domain
figure(2)
plot(x,y,'-r',x',y','-r')
xlabel('x')
ylabel('y')

% plotting distribution of metrics
figure(3)
plot3(gsi,eta,gsi_x,'-r',gsi',eta',gsi_x','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\xi_{x}')
view(36,54)

figure(4)
plot3(gsi,eta,gsi_y,'-r',gsi',eta',gsi_y','-r')
axis([0 4 0 1 0 0.5])
xlabel('\xi')
ylabel('\eta')
zlabel('\xi_{y}')
view(36,54)

figure(5)
plot3(gsi,eta,eta_x,'-r',gsi',eta',eta_x','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{x}')
set(gca,'ZDir','Reverse')
view(36,54)

figure(6)
plot3(gsi,eta,eta_y,'-r',gsi',eta',eta_y','-r')
xlabel('\xi')
ylabel('\eta')
zlabel('\eta_{y}')
view(36,54)


% figure(7)
% plot(x,y,'-r',x',y','-r')
% xlabel('x')
% ylabel('y')
% 
% figure(8)
%  plot(gsi,eta,'-ro',gsi',eta','-ro')%,'LineWidth',2,...
% %                        'MarkerEdgeColor','k',...
%                        'MarkerFaceColor','r',...
%                        'MarkerSize',7)
% set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off')
% axis([-0.25 4.25 -0.0833 1.0833])
% axis off
% % xO = 0.1;  
% % yO = 0.05;
% % axis([-0.3 4.7 -0.2 1.7])
% % text(4.5+xO,    0.0, '\xi', 'rotation', 0, 'fontsize', 17)
% % text(0,    1.6, '\eta', 'rotation', 90, 'fontsize', 14)
% % set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off')
% % 
% % % the arrows
% % 
% % patch(...
% %     [4.5-xO -xO; 4.5-xO +xO; 4.5 0.0], ...
% %     [-yO 1.5; yO 1.5; 0 1.5+yO], 'k', 'clipping', 'off')




