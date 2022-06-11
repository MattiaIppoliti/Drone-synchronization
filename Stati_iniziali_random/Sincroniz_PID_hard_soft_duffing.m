% Programma per la simulazione numerica per effettuare una sincronizzazione
% tra due oscillatori non lineare di Hard Duffing su
% sfera unitaria di I. Cervigni, M. Ippoliti, C. Menotta (UnivPM), 
% Ultima revisione: 9 Maggio 2019
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

n = 100000; % Numero di punti della simulazione
h = 0.0002; % Passo del metodo numerico

kp=10;
kd=kp;
ki_hat=kp/2;

err=zeros(3,n);
eps=zeros(3,n);
omega=zeros(3,n);
ki=zeros(1,n);
u=zeros(3,n);

%oscillatore leader
el = 1.3; % Esponenziale per il nonlinear damping
kappal = 0.5; % Coefficiente nel potenziale di Duffing
mul = 0; %0.5 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
z = zeros(3,n); w = zeros(3,n); % Preallocazione aumenta efficienza del codice
z(:,1) =randn(3,1); z(:,1) = z(:,1)/norm(z(:,1)); % Stato iniziale
rifl = randn(3,1); rifl = rifl/norm(rifl); % Riferimento per l'oscillatore
w(:,1) =0.1*(eye(3)-z(:,1)*z(:,1)')*randn(3,1); % Velocità iniziale


Kz(1) = 0.5*norm(w(:,1))^2; % Energia cinetica iniziale
Vz(1) = +(1/2)*dist(z(:,1),rifl)^2 +(1/4)*kappal*dist(z(:,1),rifl)^4; % Energia potenziale iniziale
dl(1) = dist(z(:,1),rifl); % Distanza iniziale
Hz(1) = Kz(1) + Vz(1); % Energia totale iniziale


%oscillatore follower
ef = 1.3; % Esponenziale per il nonlinear damping
kappaf = 0.5; % Coefficiente nel potenziale di Duffing
muf = 0; %0.5 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
x = zeros(3,n); v = zeros(3,n); % Preallocazione aumenta efficienza del codice
x(:,1) =randn(3,1); x(:,1) = x(:,1)/norm(x(:,1)); % Stato iniziale
riff = randn(3,1); riff = riff/norm(riff); % Riferimento per l'oscillatore
v(:,1) =0.1*(eye(3)-x(:,1)*x(:,1)')*randn(3,1); % Velocità iniziale


Kx(1) = 0.5*norm(v(:,1))^2; % Energia cinetica iniziale
Vx(1) = +(1/2)*dist(x(:,1),riff)^2 -(1/4)*kappaf*dist(x(:,1),riff)^4; % Energia potenziale iniziale
df(1) = dist(x(:,1),riff); % Distanza iniziale
Hx(1) = Kx(1) + Vx(1); % Energia totale iniziale
err(:,1)= logh(x(:,1),z(:,1));
omega(:,1)= err(1)*h;
eps(:,1)=trasp(z(:,1),x(:,1),w(:,1))- v(:,1);
ki(1)=ki_hat*(omega(:,1)'* eps(:,1));
u(:,1)= kp*err(:,1) + ki(1)* omega(:,1) + kd*eps(:,1);
sigma(1)= 1/2*norm(u(:,1))^2;

d(1)=dist(z(:,1),x(:,1)); %distanza tra i due oscillatori

% Metodo numerico (Eulero)
for k=1:n
    %leader
    z(:,k+1)= exph(z(:,k),w(:,k),h);
    w(:,k+1)= trasp(z(:,k),z(:,k+1),w(:,k) ...
                                   -h*mul*(norm(w(:,k))^(2*(el-1)))*w(:,k)...
                                   +h*(1+kappal*dist(z(:,k),rifl)^2)*logh(z(:,k),rifl)...
                                  );
    Kz(k+1)=(1/2)*norm(w(:,k+1))^2;
    Vz(k+1)=(1/2)*dist(z(:,k+1),rifl)^2+(1/4)*kappal*dist(z(:,k+1),rifl)^4;
    Hz(k+1) = Kz(k+1) + Vz(k+1);
    dl(k+1)= dist(z(:,k+1),rifl);
    
    
    
    %follower
    x(:,k+1)= exph(x(:,k),v(:,k),h);
    v(:,k+1)= trasp(x(:,k),x(:,k+1),v(:,k) ...
                                   -h*muf*(norm(v(:,k))^(2*(ef-1)))*v(:,k)...
                                   -h*(-1+kappaf*dist(x(:,k),riff)^2)*logh(x(:,k),riff)...
                                   +u(:,k)*h... 
                               );
    Kx(k+1)=(1/2)*norm(v(:,k+1))^2;
    Vx(k+1)=(1/2)*dist(x(:,k+1),riff)^2-(1/4)*kappaf*dist(x(:,k+1),riff)^4;
    Hx(k+1) = Kx(k+1) + Vx(k+1);
    df(k+1)= dist(x(:,k+1),riff);
    
    d(k+1)= dist(z(:,k+1),x(:,k+1));
    
    
    err(:,k+1)= logh(x(:,k+1),z(:,k+1));
    eps(:,k+1)=trasp(z(:,k+1),x(:,k+1),w(:,k+1))- v(:,k+1);
   
    omega(:,k+1)= trasp(x(:,k),x(:,k+1),omega(:,k)) + err(:,k+1)*h;
    %omega(:,k+1)'*x(:,k+1)
    ki(k+1)= ki_hat*(omega(:,k+1)'*eps(:,k+1));
                                    
    u(:,k+1)= kp*err(:,k+1) + ki(k+1)*omega(:,k+1) + kd*eps(:,k+1);
    sigma(k+1)= 1/2*norm(u(:,k+1))^2;
    
end

% Grafici sulla sfera unitaria nella figura 1, 
% la treaiettoria del leader è in rosso e quella del follower in verde
figure('Name','Synchronization of hard Duffing on the unit sphere');
sphere(50);alpha(0.2); % Sfera trasparente
hold on;
FS = 11; % Font size
set(gca,'FontSize',FS);xlabel('$$x_1$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$x_2$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$x_3$$','interpreter','latex');
plot3(z(1,:),z(2,:),z(3,:),'r',x(1,:),x(2,:),x(3,:),'k','LineWidth',2);
plot3(rifl(1),rifl(2),rifl(3),'bo',riff(1),riff(2),riff(3),'yo','LineWidth',2);
quiver3(z(1),z(2),z(3),w(1),w(2),w(3),'wo','LineWidth',2,'ShowArrowHead','on','AutoScale','on');
quiver3(x(1),x(2),x(3),v(1),v(2),v(3),'go','LineWidth',2,'ShowArrowHead','on','AutoScale','on');
drawnow;

%grafico energia alta qualità
print Sincr_HardSoft_random -dpdf -r500;

% Grafici delle energie cinetica, potenziale e totale nella figura 2 
t = 0:h:n*h; 
FS = 11; % Font size
figure('Name','Synchronization of hard Duffing on the unit sphere');
subplot(3,2,1); plot(t,Kz,'r-',t,Kx,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Kinetic energy $$K$$','interpreter','latex');
subplot(3,2,2); plot(t,Vz,'r-',t,Vx,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Potential energy $$V$$','interpreter','latex');
subplot(3,2,3); plot(t,Hz,'r-',t,Hx,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Total energy $$H$$','interpreter','latex');
subplot(3,2,4); plot(t,dl,'r-',t,df,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distances $$d(x,r)$$, $$d(z,r)$$','interpreter','latex');
subplot(3,2,5); semilogy(t,d,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(z,x)$$','interpreter','latex');
subplot(3,2,6); semilogy(t,sigma,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Control effort $$\sigma$$','interpreter','latex');

%grafico energia alta qualità
print Sincr_HardSoft_random_energy -dpdf -r500;
