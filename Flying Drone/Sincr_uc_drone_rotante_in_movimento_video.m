% Programma per la simulazione numerica di un drone, per il momento capace di fare soltanto movimenti rotazionali,
% implementato sul Manifold delle iperrotazioni SO(3) di I. Cervigni, M. Ippoliti, C. Menotta (UnivPM), 
% Ultima revisione: 24/04/2019
clear all; close all; clc; 

%caricamento figura drone
load drone
figure;
Droni = [Xr';Yr';Zr'];  %le coordinate dei punti che costituiscono la figura del drone
close all
Centrorot = [-0.0125  -0.0048 0]'; %centro di massa del nostro drone

%Lie bracket
brak = @(x,y) x*y-y*x;
%matrix anti-commutator
anticom = @(x,y) x*y+y*x;
% Distanza Riemanniana SO(3)
dist_SO3 = @(x,y) norm(realLogSO(x'*y),'fro');
% Distanza R^3
dist_R3 = @(x,y) norm(x-y); 
% vee
vee = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];

%opengl software % Per le trasparenze
n = 6000; % Numero di punti della simulazione
h=0.001; % passo di campionamento


%PID SO(3)
kp=50;
kd=kp;
ki_hat=10;

err=zeros(3,3,n);
eps=zeros(3,3,n);
omega=zeros(3,3,n);
ki=zeros(1,n);
u=zeros(3,3,n);
uc=zeros(3,3,n);
sigma=zeros(1,n);
sigma_uc=zeros(1,n);
sigma_PID=zeros(1,n);
d_SO3=zeros(1,n);


%PID R^3
kp2=100;
kd2=80;
ki_hat2=10;

err2=zeros(3,n);
eps2=zeros(3,n);
omega2=zeros(3,n);
ki2=zeros(1,n);
uc2=zeros(3,n);
alfa=zeros(3,n);
sigma2=zeros(1,n);
sigma_uc2=zeros(1,n);
sigma_PID2=zeros(1,n);
d_R3=zeros(1,n);

%costanti dei droni follower e leader 
Jx=7.5*10^(-3); Jy=Jx; Jz=1.3*10^(-2);
Jq=(1/2)*diag([Jy-Jx+Jz Jx-Jy+Jz Jx+Jy-Jz]);
D=diag([sqrt((Jy*Jz)/Jx) sqrt((Jx*Jz)/Jy) sqrt((Jx*Jy)/Jz)]);
inD= inv(D);

b= 3.13*10^(-5);
r=0.23;
y=7.15*10^(-7);
p=9.81;
Mq=0.650;

Wx=[0 0 0; 0 0 -1;0 1 0];
Wy=[0 0 1; 0 0 0; -1 0 0];
Wz=[ 0 -1 0; 1 0 0; 0 0 0];

ez = [0;0;1];

nh = 3048;

G = 1/4 * eye(3);

%leader 
%Yaw & Pitch(+x)
RPM1 = 3047; RPM2 = 3048; RPM3 = 3049; RPM4 = 3048;
w1 = RPM1*pi/30; w2=RPM2*pi/30; w3=RPM3*pi/30; w4=RPM4*pi/30;

Z = zeros(3,3,n); W = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
Z(:,:,1)=eye(3); % Stato iniziale

z = zeros(3,n); w = zeros(3,n);

taul=b*r*(w4^2-w2^2)*Wx+b*r*(w3^2-w1^2)*Wy+y*(-w1^2+w2^2-w3^2+w4^2)*Wz;
Wr=-w1+w2-w3+w4;

%follower
X = zeros(3,3,n); V = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
X(:,:,1)= 0.1*[-0.8 4 3; 3 -0.9 0; 0 2 -1]; X(:,:,1)=expm(X(:,:,1)'-X(:,:,1));%rand(3,3); X(:,:,1)=expm(X(:,:,1)'-X(:,:,1)); % Stato iniziale

x = [5;3;3]; v = zeros(3,n);

%PID R^3
err2(:,1)= z(:,1)-x(:,1)-[0;0;5];
omega2(:,1)= err2(:,1)*h;
eps2(:,1)=w(:,1)- v(:,1);
%ki2(1)=ki_hat2*(omega2(:,1)'* eps2(:,1));
ki2(1)=ki_hat2;
uc2(:,1)= +p*ez+(1/Mq)*G*v(:,1);
alfa(:,1)= kp2*err2(:,1) + ki2(1)* omega2(:,1) + kd2*eps2(:,1)+uc2(:,1);
sigma2(1)= 1/2*norm(alfa(:,1))^2;
sigma_uc2(1)= 1/2*norm(uc2(:,1))^2;
sigma_PID2(1)= 1/2*norm(alfa(:,1)-uc2(:,1))^2;

d_R3(1) = dist_R3(z(:,1),x(:,1));

%PID SO(3)
err(:,:,1)= realLogSO(X(:,:,1)'*Z(:,:,1)) +  0.0001*d_R3(1)*vee(cross(X(:,:,1)*ez, alfa(:,1)));
omega(:,:,1)= err(:,:,1)*h;
eps(:,:,1)=W(:,:,1)- V(:,:,1);
ki(1)=ki_hat*trace(omega(:,:,1)'* eps(:,:,1));
uc(:,:,1)= inD*brak(Jq,V(:,:,1)^2)*inD;
u(:,:,1)= kp*err(:,:,1) + ki(1)* omega(:,:,1) + kd*eps(:,:,1) + uc(:,:,1);
sigma(1)= 1/2*norm(u(:,:,1),'fro')^2;
sigma_uc(1)= 1/2*norm(uc(:,:,1),'fro')^2;
sigma_PID(1)= 1/2*norm(u(:,:,1)-uc(:,:,1),'fro')^2;
d_SO3(1) = dist_SO3(Z(:,:,1),X(:,:,1));

% Metodo numerico (Eulero)
disp('Simulazione numerica...')
for k=1:n
    %leader
    Z(:,:,k+1) = Z(:,:,k)* expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) + h*inD*(brak(Jq,W(:,:,k)^2)+taul)*inD;
    
    z(:,k+1) = z(:,k)+ h*w(:,k);
    w(:,k+1) = w(:,k)+ h*((1/2)*(b/Mq)*(w1^2+w2^2+w3^2+w4^2)*Z(:,:,k)*ez-p*ez-(1/Mq)*G*w(:,k));
    
    %follower
    X(:,:,k+1) = X(:,:,k)* expm(h*V(:,:,k));
    V(:,:,k+1) = V(:,:,k) + h*inD*(brak(Jq,V(:,:,k)^2))*inD + h*u(:,:,k);
    
    x(:,k+1) = x(:,k)+ h*v(:,k);
    v(:,k+1) = v(:,k)+ h*((alfa(:,k)'*(X(:,:,k)*ez))*(X(:,:,k)*ez)-p*ez-(1/Mq)*G*v(:,k));
    
    
    
    %PID R^3
    
    err2(:,k+1)= z(:,k+1)-x(:,k+1)-[0;0;5];
    omega2(:,k+1)= omega(:,k)+ err2(:,k+1)*h;
    eps2(:,k+1)=w(:,k+1)- v(:,k+1);
    %ki2(k+1)=ki_hat2*(omega2(:,k+1)'* eps2(:,k+1));
    ki2(k+1) = ki_hat2;
    uc2(:,k+1)= (w(:,k+1)-w(:,k))/h+p*ez+(1/Mq)*G*v(:,k+1);
    alfa(:,k+1)= kp2*err2(:,k+1) + ki2(k+1)* omega2(:,k+1) + kd2*eps2(:,k+1)+uc2(:,k+1);
    sigma2(k+1)= 1/2*norm(alfa(:,k+1))^2;
    sigma_uc2(k+1)= 1/2*norm(uc2(:,k+1))^2;
    sigma_PID2(k+1)= 1/2*norm(alfa(:,k+1)-uc2(:,k+1))^2;

    d_R3(k+1) = dist_R3(z(:,k+1),x(:,k+1));
    
    %PID SO(3)
    
    %err(:,:,k+1)= realLogSO(X(:,:,k+1)'*Z(:,:,k+1)) - 0.00001*vee(alfa(:,k));
    %err(:,:,k+1)= realLogSO(X(:,:,k+1)'*Z(:,:,k+1)) - 0.00001*vee(cross(z(:,k+1)-x(:,k+1)-[0;0;5], alfa(:,k)));
    %err(:,:,k+1)= 0.01*vee(cross(z(:,k+1)-x(:,k+1)-[0;0;5], X(:,:,k+1)*ez));
    m(:,:,k+1) = 0.001*vee(cross(X(:,:,k+1)*ez, alfa(:,k+1)));
    err(:,:,k+1)= realLogSO(X(:,:,k+1)'*Z(:,:,k+1))+ m(:,:,k+1);
    
    omega(:,:,k+1)=omega(:,:,k) + err(:,:,k+1)*h;
    eps(:,:,k+1)=W(:,:,k+1)- V(:,:,k+1)+ (m(:,:,k+1)-m(:,:,k))/h;
    ki(k+1)=ki_hat*trace(omega(:,:,k+1)'* eps(:,:,k+1));
    uc(:,:,k+1)= (W(:,:,k+1)-W(:,:,k))/h...
                 -inD*brak(Jq,V(:,:,k)^2)*inD;
             
    u(:,:,k+1)= kp*err(:,:,k+1) + ki(k+1)* omega(:,:,k+1) + kd*eps(:,:,k+1) + uc(:,:,k+1);
    sigma(k+1)= 1/2*norm(u(:,:,k+1),'fro')^2;
    sigma_uc(k+1)= 1/2*norm(uc(:,:,k+1),'fro')^2;
    sigma_PID(k+1)= 1/2*norm(u(:,:,k+1)-uc(:,:,k+1),'fro')^2;

    d_SO3(k+1) = dist_SO3(Z(:,:,k+1),X(:,:,k+1));
    
end

figure;
disp('Animazione grafica 3D...')

%La rappresentazione del droneo viene messa in un altro ciclo, poichè non si
%prendono tutti i campioni ma se ne prende solo 1 ogni 100, in maniera da rendere la
%rappresentazione più veloce
DroniRuotati_leader = Z(:,:,1)*(Droni - Centrorot) + z(:,1);        
DroniRuotati_follower = X(:,:,1)*(Droni - Centrorot) + x(:,1);        
hh=plot3(DroniRuotati_leader(1,:),DroniRuotati_leader(2,:),DroniRuotati_leader(3,:),'ro-',...
         DroniRuotati_follower(1,:),DroniRuotati_follower(2,:),DroniRuotati_follower(3,:),'bo-');
grid on;
axis([-5 45 -10 10 -35 5])
FS = 14; % Font size
set(gca,'FontSize',FS);xlabel('$$x$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$y$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$z$$','interpreter','latex');
hold on

%Video
%v=VideoWriter('Drone_nomix.mov');
%open(v);

for k=1:10:n
    delete(hh);
    DroniRuotati_leader = Z(:,:,k+1)*(Droni - Centrorot) + z(:,k+1);        
    DroniRuotati_follower = X(:,:,k+1)*(Droni - Centrorot) + x(:,k+1);               
    hh=plot3(DroniRuotati_leader(1,:),DroniRuotati_leader(2,:),DroniRuotati_leader(3,:),'ro-',...
             DroniRuotati_follower(1,:),DroniRuotati_follower(2,:),DroniRuotati_follower(3,:),'bo-');
    drawnow
    
    %video
    %frame = getframe(gcf);
    %writeVideo(v,frame);
end

%close(v);

%grafici della distanza e sforzo di controllo
t = 0:h:n*h; 
FS = 11; % Font size
figure('Name','Quadcopter drone');
subplot(2,2,1); semilogy(t,d_SO3,'b-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(Z,X)$$','interpreter','latex');
subplot(2,2,2); semilogy(t,sigma,'b-',t,sigma_uc,'g-',t,sigma_PID,'m-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Control effort $$\sigma_{SO(3)}$$','interpreter','latex');
                    legend('\sigma','\sigma_C', '\sigma_{PD}','Location','northoutside', 'Orientation','horizontal');
                    legend('boxoff')

subplot(2,2,3); semilogy(t,d_R3,'b-',t,5*ones(size(d_R3)),'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(z,x)$$','interpreter','latex');
subplot(2,2,4); semilogy(t,sigma2,'b-',t,sigma_uc2,'g-',t,sigma_PID2,'m-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Control effort $$\sigma_{(R^3)}$$','interpreter','latex');
                    legend('\sigma','\sigma_C', '\sigma_{PD}','Location','northoutside', 'Orientation','horizontal');
                    legend('boxoff')
