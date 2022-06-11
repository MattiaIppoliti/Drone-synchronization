% Programma per la simulazione numerica di un oscillatore Hard Duffing su
% Manifold delle iperrotazioni di I. Cervigni, M. Ippoliti, C. Menotta (UnivPM), 
% Ultima revisione: 16 Aprile 2019
clear all; close all; clc; 

% Mappa esponenziale su SO(3)
exph = @(x,v,h) x*expm(h*x*v); 
% Logarithm
logh = @(x,y) x*realLogSO(x'*y);
% Distanza Riemanniana
dist = @(x,y) norm(realLogSO(x'*y),'fro');

%opengl software % Per le trasparenze

Spigoli = [ 0 1 1 0 0 0 1 1 0 0 0 0 0 0 0 1 1 1 1 
            0 0 0 0 0 1 1 0 0 1 1 0 0 1 1 1 0 1 1 
            0 0 1 1 0 0 0 0 0 0 1 1 0 0 1 1 1 1 0];
        
Polo = [0.5 0.5 0.5]'; % Il cubo ruoterà intorno a questo punto, può essere sia interno che esterno al cubo

n = 50000; % Numero di punti della simulazione
h = 0.0002; % Passo del metodo numerico

kp=10;
kd=kp;
ki_hat=kp/2;

err=zeros(3,3,n);
eps=zeros(3,3,n);
omega=zeros(3,3,n);
ki=zeros(1,n);
u=zeros(3,3,n);
uc=zeros(3,3,n+1);
sigma=zeros(1,n);
sigma_c=zeros(1,n);
sigma_PID=zeros(1,n);

%leader
el = 1.3; % Esponenziale per il nonlinear damping
kappal = 0.5; % Coefficiente nel potenziale di Duffing
mul = 0; %0.5 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
Z = zeros(3,3,n); W = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
Z(:,:,1)= rand(3,3); Z(:,:,1)=expm(Z(:,:,1)'-Z(:,:,1)); % Stato iniziale
W(1,2,1)=randn; W(2,1,1)=-W(1,2,1); W(1,3,1)=randn; W(3,1,1)=-W(1,3,1);
W(2,3,1)=randn; W(3,2,1)=-W(2,3,1);  %velocità iniziale
rifl = rand(3,3); rifl(:,:,1)=expm(rifl(:,:,1)'-rifl(:,:,1));
%{
q=zeros(3,n); p=randn(3,1); p=p/norm(p); %p indica il punto iniziale 
rifP=rif*p; % Riferimento per l'oscillatore
%}
KZ(1) = (1/2)*norm(Z(:,:,1)*W(:,:,1),'fro')^2; % Energia cinetica iniziale
VZ(1) = (1/2)*dist(Z(:,:,1),rifl)^2+(1/4)*kappal*dist(Z(:,:,1),rifl)^4; % Energia potenziale iniziale
dl(1) = dist(Z(:,:,1),rifl); % Distanza iniziale
HZ(1) = KZ(1) + VZ(1); % Energia totale iniziale



%follower
ef = 1.3; % Esponenziale per il nonlinear damping
kappaf = 0.5; % Coefficiente nel potenziale di Duffing
muf = 0; %0.5 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
X = zeros(3,3,n); V = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
X(:,:,1)= rand(3,3); X(:,:,1)=expm(X(:,:,1)'-X(:,:,1)); % Stato iniziale
V(1,2,1)=randn; V(2,1,1)=-V(1,2,1); V(1,3,1)=randn; V(3,1,1)=-V(1,3,1);
V(2,3,1)=randn; V(3,2,1)=-V(2,3,1);  %velocità iniziale
riff = rand(3,3); riff(:,:,1)=expm(riff(:,:,1)'-riff(:,:,1));
%{
q=zeros(3,n); p=randn(3,1); p=p/norm(p); %p indica il punto iniziale 
rifP=rif*p; % Riferimento per l'oscillatore
%}
KX(1) = (1/2)*norm(X(:,:,1)*V(:,:,1),'fro')^2; % Energia cinetica iniziale
VX(1) = (1/2)*dist(X(:,:,1),riff)^2+(1/4)*kappaf*dist(X(:,:,1),riff)^4; % Energia potenziale iniziale
df(1) = dist(X(:,:,1),riff); % Distanza iniziale
HX(1) = KX(1) + VX(1); % Energia totale iniziale

err(:,:,1)= realLogSO(X(:,:,1)'*Z(:,:,1));
omega(:,:,1)= err(:,:,1)*h;
eps(:,:,1)=W(:,:,1)- V(:,:,1);
ki(1)=ki_hat*trace(omega(:,:,1)'* eps(:,:,1));
uc(:,:,1)= -(-muf*norm(V(:,:,1),'fro')^(2*(ef-1))*V(:,:,1)...
           +(1+kappaf*dist(X(:,:,1),riff)^2)*realLogSO(X(:,:,1)'*riff));
u(:,:,1)= kp*err(:,:,1) + ki(1)* omega(:,:,1) + kd*eps(:,:,1) + uc(:,:,1);
sigma(1)= 1/2*norm(u(:,1))^2;
sigma_c(1)= 1/2*norm(uc(:,1))^2;
sigma_PID(1) = 1/2*norm(u(:,1)-uc(:,1))^2;

d(1) = dist(Z(:,:,1),X(:,:,1));

% Metodo numerico (Eulero)
disp('Simulazione numerica...')
for k=1:n
    %leader
    Z(:,:,k+1) = Z(:,:,k)*expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) - h*mul*norm(W(:,:,k),'fro')^(2*(el-1))*W(:,:,k)...
                          + h*(1+kappal*dist(Z(:,:,k),rifl)^2)*realLogSO(Z(:,:,k)'*rifl);
    
    %q(:,k)= X(:,:,k)*p;
    
    KZ(k+1)=(1/2)*norm(Z(:,:,k+1)*W(:,:,k+1),'fro')^2;
    VZ(k+1)=(1/2)*dist(Z(:,:,k+1),rifl)^2+(1/4)*kappal*dist(Z(:,:,k+1),rifl)^4;
    HZ(k+1) = KZ(k+1) + VZ(k+1);
    dl(k+1) = dist(Z(:,:,k+1),rifl);
    
    %follower
    X(:,:,k+1) = X(:,:,k)*expm(h*V(:,:,k));
    V(:,:,k+1) = V(:,:,k) - h*muf*norm(V(:,:,k),'fro')^(2*(ef-1))*V(:,:,k)...
                          + h*(1+kappaf*dist(X(:,:,k),riff)^2)*realLogSO(X(:,:,k)'*riff)...
                          + h*u(:,:,k);
    
    %q(:,k)= X(:,:,k)*p;
    
    KX(k+1)=(1/2)*norm(X(:,:,k+1)*V(:,:,k+1),'fro')^2;
    VX(k+1)=(1/2)*dist(X(:,:,k+1),riff)^2+(1/4)*kappaf*dist(X(:,:,k+1),riff)^4;
    HX(k+1) = KX(k+1) + VX(k+1);
    df(k+1) = dist(X(:,:,k+1),riff);
    
    err(:,:,k+1)= realLogSO(X(:,:,k+1)'*Z(:,:,k+1));
    omega(:,:,k+1)=omega(:,:,k) + err(:,:,k+1)*h;
    eps(:,:,k+1)=W(:,:,k+1)- V(:,:,k+1);
    ki(k+1)=ki_hat*trace(omega(:,:,k+1)'* eps(:,:,k+1));
    uc(:,:,k+1)= (W(:,:,k+1)-W(:,:,k))/h...
                 -(-muf*norm(V(:,:,k+1),'fro')^(2*(ef-1))*V(:,:,k+1)...
                 +(1+kappaf*dist(X(:,:,k+1),riff)^2)*realLogSO(X(:,:,k+1)'*riff));
             
    u(:,:,k+1)= kp*err(:,:,k+1) + ki(k+1)* omega(:,:,k+1) + kd*eps(:,:,k+1) + uc(:,:,k+1);
    sigma(k+1)= 1/2*norm(u(:,:,k+1),'fro')^2;
    sigma_c(k+1)=1/2*norm(uc(:,k+1))^2;
    sigma_PID(k+1) = 1/2*norm(u(:,k+1)-uc(:,k+1))^2;
    
    d(k+1)= dist(Z(:,:,k+1),X(:,:,k+1));
    
end

%La rappresentazione del cubo viene messa in un altro ciclo, poichè non si
%prendono tutti i campioni ma se ne prende solo 1 ogni 20, in maniera da rendere la
%rappresentazione più veloce
disp('Animazione grafica 3D...')
FS=11;
figure('Name','hard duffing oscillator applied to a cube');
set(gca,'FontSize',FS);xlabel('$$x$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$y$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$z$$','interpreter','latex');
%leader
SpigoliRuotati_l = Z(:,:,1)*(Spigoli - Polo) + Polo;
    
%follower
SpigoliRuotati_f = X(:,:,1)*(Spigoli - Polo) + Polo;
      
hh = plot3(SpigoliRuotati_l(1,:),SpigoliRuotati_l(2,:),SpigoliRuotati_l(3,:),'ro-',...
          SpigoliRuotati_f(1,:),SpigoliRuotati_f(2,:),SpigoliRuotati_f(3,:),'ko-',...           
          Polo(1),Polo(2),Polo(3),'bo');
grid on; axis([-2 2 -2 2 -2 2]);
hold on

%Video
v=VideoWriter('cubosfera.mov');
open(v);

for k=1:20:n
    %leader
    SpigoliRuotati_l = Z(:,:,k+1)*(Spigoli - Polo) + Polo;
    
    %follower
    SpigoliRuotati_f = X(:,:,k+1)*(Spigoli - Polo) + Polo;
    
    delete(hh);
    hh =  plot3(SpigoliRuotati_l(1,:),SpigoliRuotati_l(2,:),SpigoliRuotati_l(3,:),'ro-',...
          SpigoliRuotati_f(1,:),SpigoliRuotati_f(2,:),SpigoliRuotati_f(3,:),'ko-',...           
          Polo(1),Polo(2),Polo(3),'bo');
    drawnow
    %video
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

%grafico alta qualità
print Sincr_SO(3)_uc -dpdf -r500;

%{
% Grafici sulla sfera unitaria nella figura 1
figure('Name','hard duffing oscillator on the unit sphere');
sphere(50);alpha(0.2); % Sfera trasparente
hold on;
FS = 11; % Font size
set(gca,'FontSize',FS);xlabel('$$x$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$y$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$z$$','interpreter','latex');
plot3(q(1,:),q(2,:),q(3,:),'r','LineWidth',2);  %oscillazione sulla sfera
plot3(rifP(1),rifP(2),rifP(3),'bo','LineWidth',2);  %punto di riferimento dell'oscillazione
plot3(p(1),p(2),p(3),'ro','LineWidth',2); %punto iniziale dell'oscillazione
drawnow;
%}

% Grafici delle energie cinetica, potenziale e totale nella figura 2 
t = 0:h:n*h; 
FS = 11; % Font size
figure('Name','hard duffing oscillator on the unit sphere');
subplot(3,2,1); plot(t,KZ,'r-',t,KX,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Kinetic energy $$K$$','interpreter','latex');
subplot(3,2,2); plot(t,VZ,'r-',t,VX,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Potential energy $$V$$','interpreter','latex');
subplot(3,2,3); plot(t,HZ,'r-',t,HX,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Total energy $$H$$','interpreter','latex');
subplot(3,2,4); plot(t,dl,'r-',t,df,'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distances $$d(Z,R)$$ et $$d(X,R)$$' ,'interpreter','latex');
subplot(3,2,5); semilogy(t,d,'b-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(Z,X)$$','interpreter','latex');
subplot(3,2,6);  semilogy(t,sigma,'b-',t,sigma_c,'g-',t,sigma_PID,'m-'); axis tight
                 legend('\sigma','\sigma_C', '\sigma_{PID}','Location','northoutside','Orientation','horizontal');
                 legend('boxoff');
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Control effort $$\sigma$$','interpreter','latex');

%grafico energia alta qualità
print Sincr_SO(3)_uc_energy -dpdf -r500;

