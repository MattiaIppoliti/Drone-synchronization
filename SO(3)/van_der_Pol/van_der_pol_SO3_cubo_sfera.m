% Programma per la simulazione numerica di un oscillatore nonlineare di van der Pol su
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

n = 100000; % Numero di punti della simulazione
h = 0.0005; % Passo del metodo numerico
e = 1.3; %(1.3 non-damped) (3 damped) Esponenziale per il nonlinear damping
kappa = 1; % Coefficiente nel potenziale di van der Pol
mu = 0; %0.5 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
X = zeros(3,3,n); W = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
q=zeros(3,n); p=randn(3,1); p=p/norm(p); %p indica il punto iniziale 
X(:,:,1)=eye(3); % Stato iniziale
rif = rand(3); X(:,:,1)=eye(3); rif=expm(rif'-rif); rifP=rif*p; % Riferimento per l'oscillatore
W(:,:,1) = [0 -0.5 -1; 0.5 0 0.7; 1 -0.7 0]; 

Kx(1) = (1/2)*norm(X(:,:,1)*W(:,:,1),'fro')^2; % Energia cinetica iniziale
Vx(1) = (1/2)*kappa*dist(X(:,:,1),rif)^2; % Energia potenziale iniziale
d(1) = dist(X(:,1),rif); % Distanza iniziale
Hx(1) = Kx(1) + Vx(1); % Energia totale iniziale

% Metodo numerico (Eulero)
for k=1:n
    X(:,:,k+1) = X(:,:,k)*expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) - h*mu*norm(W(:,:,k),'fro')^(2*(e-1))*W(:,:,k) ...
                          + h*kappa*realLogSO(X(:,:,k)'*rif);
    
    
    q(:,k)= X(:,:,k)*p;
    
    Kx(k+1)=(1/2)*norm(X(:,:,k+1)*W(:,:,k+1),'fro')^2;
    Vx(k+1)=(1/2)*kappa*dist(X(:,:,k+1),rif)^2;
    Hx(k+1) = Kx(k+1) + Vx(k+1);
    d(k+1) = dist(X(:,:,k+1),rif);
end

%La rappresentazione del cubo viene messa in un altro ciclo, poichè non si
%prendono tutti i campioni ma se ne prende solo 1 ogni 20, in maniera da rendere la
%rappresentazione più veloce

for k=1:20:n
    SpigoliRuotati = X(:,:,k+1)*(Spigoli - Polo) + Polo;
    plot3(SpigoliRuotati(1,:),SpigoliRuotati(2,:),SpigoliRuotati(3,:),'ro-',Polo(1),Polo(2),Polo(3),'bo')
    grid on; axis([-2 2 -2 2 -2 2]) 
    xlabel('x');ylabel('y');zlabel('z')
    drawnow;

end

% Grafici sulla sfera unitaria nella figura 1
figure('Name','Van der Pol oscillator on the unit sphere');
sphere(50);alpha(0.2); % Sfera trasparente
hold on;
FS = 11; % Font size
set(gca,'FontSize',FS);xlabel('$$x_1$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$x_2$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$x_3$$','interpreter','latex');

plot3(q(1,:),q(2,:),q(3,:),'r','LineWidth',2);  %oscillazione sulla sfera
plot3(rifP(1),rifP(2),rifP(3),'bo','LineWidth',2);  %punto di riferimento dell'oscillazione
plot3(p(1),p(2),p(3),'ro','LineWidth',2); %punto iniziale dell'oscillazione
drawnow;

% Grafici delle energie cinetica, potenziale e totale nella figura 2 
t = 0:h:n*h; 
FS = 11; % Font size
figure('Name','Van der Pol oscillator on the unit sphere');
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
set(gca,'FontSize',FS);ylabel('Distance $$d(X,R)$$','interpreter','latex');

