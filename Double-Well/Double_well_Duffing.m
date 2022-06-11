% Programma per la simulazione numerica di un oscillatore nonlineare di Double-well Duffing su
% sfera unitaria di I. Cervigni, M. Ippoliti, C. Menotta (UnivPM), 
% Ultima revisione: 24 Aprile 2019 (Fiori)
clear all; close all; clc; 

% Mappa esponenziale sulla sfera
exph = @(x,v,h) x*cos(h*norm(v)) + h*v*sinc(h*norm(v)/pi); 
% Trasporto parallelo sulla sfera
trasp = @(x,y,w) w - ( y*(y'*w) - x*(x'*y)*(y'*w) ) / (1+x'*y) - x*(y'*w);
% Distanza Riemanniana
dist = @(x,y) abs(acos(x'*y));
% Logarithm
logh = @(x,y) (y - x*(x'*y)) / sinc(dist(x,y)/pi);

%opengl software % Per le trasparenze

n = 500000; % (n=100000 nondamped)(n=500000 damped) Numero di punti della simulazione
h = 0.0001; % Passo del metodo numerico
e = 1.3; % Esponenziale per il nonlinear damping
kappa = 1.4; %(kappa=0.8 nondamped)(kappa=1.4 damped) Coefficiente nel potenziale di Duffing
mu = 0; %mu=0.2 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
x = zeros(3,n); v = zeros(3,n); % Preallocazione aumenta efficienza del codice
x(:,1) = [0; 0; 1]; %rand(3,1); x(:,1) = x(:,1)/norm(x(:,1)); % Stato iniziale
rif = [1; 0; 0]; % Riferimento per l'oscillatore
v(:,1) = [-1; 1.5; 0]; % Velocità iniziale

Kx(1) = 0.5*norm(v(:,1))^2; % Energia cinetica iniziale
Vx(1) = -(1/2)*dist(x(:,1),rif)^2 + (1/4)*kappa*dist(x(:,1),rif)^4; % Energia potenziale iniziale
d(1) = dist(x(:,1),rif); % Distanza iniziale
Hx(1) = Kx(1) + Vx(1); % Energia totale iniziale

% Metodo numerico (Eulero)
for k=1:n
    x(:,k+1)= exph(x(:,k),v(:,k),h);
    v(:,k+1)= trasp(x(:,k),x(:,k+1),v(:,k) ...
                                   -h*mu*(norm(v(:,k))^(2*(e-1)))*v(:,k)...
                                   +h*(-1+kappa*dist(x(:,k),rif)^2)*logh(x(:,k),rif)...
                   );
    d(k+1)= dist(x(:,k+1),rif);
    Kx(k+1)=(1/2)*norm(v(:,k+1))^2;
    Vx(k+1)=-(1/2)*d(k+1)^2+(1/4)*kappa*d(k+1)^4; % <- Qui erano invertiti i segni, ecco perché H veniva variabile!
    Hx(k+1) = Kx(k+1) + Vx(k+1);
end

% Grafici sulla sfera unitaria nella figura 1
figure('Name','Double-Well Duffing on the unit sphere');
sphere(50);alpha(0.2); % Sfera trasparente
hold on;
FS = 11; % Font size
set(gca,'FontSize',FS);xlabel('$$x_1$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$x_2$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$x_3$$','interpreter','latex');
plot3(x(1,:),x(2,:),x(3,:),'r','LineWidth',2);
plot3(rif(1),rif(2),rif(3),'bo','LineWidth',2);
drawnow;

% Grafici delle energie cinetica, potenziale e totale nella figura 2 
t = 0:h:n*h; 
FS = 11; % Font size
figure('Name','Double-Well Duffing on the unit sphere');
subplot(2,2,1); plot(t,Kx,'b-',t,Kx(1)*ones(size(t)),'r--'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Kinetic energy $$K$$','interpreter','latex');
subplot(2,2,2); plot(t,Vx,'b-',t,Vx(1)*ones(size(t)),'r--'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Potential energy $$V$$','interpreter','latex');
subplot(2,2,3); plot(t,Hx,'b-',t,(Hx(1))*ones(size(t)),'r--'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Total energy $$H$$','interpreter','latex');
subplot(2,2,4); plot(t,d,'b-',t,d(1)*ones(size(t)),'r--'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(x,r)$$','interpreter','latex');




