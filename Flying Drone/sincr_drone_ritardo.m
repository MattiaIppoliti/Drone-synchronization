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
n = 12000; % Numero di punti della simulazione
h=0.001; % passo di campionamento
late=2000;

%PID SO(3)
kp=150;
kd=kp;
ki_hat=kp/2;

err=zeros(3,3,n);
eps=zeros(3,3,n);
omega=zeros(3,3,n);
ki=zeros(1,n);
u=zeros(3,3,n);
uc=zeros(3,3,n);
sigma=zeros(1,n);
sigma_uc=zeros(1,n);
d_SO3=zeros(1,n);


%PID R^3
kp2=50;
kd2=2*kp2;
ki_hat2=kp2/2;

err2=zeros(3,n);
eps2=zeros(3,n);
omega2=zeros(3,n);
ki2=zeros(1,n);
uc2=zeros(3,n);
alfa=zeros(3,n);
sigma2=zeros(1,n);
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
%Pitch (+x) & Roll (+y)
RPM1=3047; RPM2=3049; RPM3=3049; RPM4=3047;
w1(1:(n/3)) = RPM1*pi/30; w2(1:(n/3)) =RPM2*pi/30; w3(1:(n/3)) =RPM3*pi/30; w4(1:(n/3)) =RPM4*pi/30;

%Yaw & Pitch(+x)
RPM1 = 3047; RPM2 = 3048; RPM3 = 3049; RPM4 = 3048;
w1((n/3)+1:(2*n/3)) = RPM1*pi/30; w2((n/3)+1:(2*n/3)) =RPM2*pi/30; w3((n/3)+1:(2*n/3)) =RPM3*pi/30; w4((n/3)+1:(2*n/3)) =RPM4*pi/30;

%Yaw & Roll(-y)
RPM1=3048; RPM2=3047; RPM3=3048; RPM4=3049;
w1((2*n/3)+1:n) = RPM1*pi/30; w2((2*n/3)+1:n)  =RPM2*pi/30; w3((2*n/3)+1:n)  =RPM3*pi/30; w4((2*n/3)+1:n)=RPM4*pi/30;

Z = zeros(3,3,n); W = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
Z(:,:,1)=eye(3); % Stato iniziale

z = zeros(3,n); w = zeros(3,n);

Wr=-w1+w2-w3+w4;

%follower
X = zeros(3,3,n); V = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
X(:,:,1+late)= rand(3,3); X(:,:,1+late)=expm(X(:,:,1+late)'-X(:,:,1+late)); % Stato iniziale

x = zeros(3,n);
x(:,1+late) = [5;3;3]; v = zeros(3,n);

%PID R^3
err2(:,1+late)= z(:,1)-x(:,1+late)-[0;0;5];
omega2(:,1+late)= err2(:,1+late)*h;
eps2(:,1+late)=w(:,1)- v(:,1+late);
%ki2(1)=ki_hat2*(omega2(:,1)'* eps2(:,1));
ki2(1+late)=ki_hat2;
uc2(:,1+late)= +p*ez+(1/Mq)*G*v(:,1+late);
alfa(:,1+late)= kp2*err2(:,1+late) + ki2(1+late)* omega2(:,1+late) + kd2*eps2(:,1+late)+uc2(:,1+late);
sigma2(1+late)= 1/2*norm(alfa(:,1+late))^2;
d_R3(1+late) = dist_R3(z(:,1),x(:,1+late));

%PID SO(3)
err(:,:,1+late)= realLogSO(X(:,:,1+late)'*Z(:,:,1)) +  0.0001*d_R3(1+late)*vee(cross(X(:,:,1+late)*ez, alfa(:,1+late)));
omega(:,:,1+late)= err(:,:,1+late)*h;
eps(:,:,1+late)=W(:,:,1)- V(:,:,1+late);
ki(1+late)=ki_hat*trace(omega(:,:,1+late)'* eps(:,:,1+late));
uc(:,:,1+late)= inD*brak(Jq,V(:,:,1+late)^2)*inD;
u(:,:,1+late)= kp*err(:,:,1+late) + ki(1+late)* omega(:,:,1+late) + kd*eps(:,:,1+late) + uc(:,:,1+late);
sigma(1+late)= 1/2*norm(u(:,:,1+late),'fro')^2;
sigma_uc(1+late)= 1/2*norm(uc(:,:,1+late),'fro')^2;
d_SO3(1+late) = dist_SO3(Z(:,:,1),X(:,:,1+late));

% Metodo numerico (Eulero)
disp('Simulazione numerica...')
for k=1:n
    %leader
    taul(:,:,k)=b*r*(w4(k)^2-w2(k)^2)*Wx+b*r*(w3(k)^2-w1(k)^2)*Wy+y*(-w1(k)^2+w2(k)^2-w3(k)^2+w4(k)^2)*Wz;
    Z(:,:,k+1) = Z(:,:,k)* expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) + h*inD*(brak(Jq,W(:,:,k)^2)+taul(:,:,k))*inD;
    
    z(:,k+1) = z(:,k)+ h*w(:,k);
    w(:,k+1) = w(:,k)+ h*((1/2)*(b/Mq)*(w1(k)^2+w2(k)^2+w3(k)^2+w4(k)^2)*Z(:,:,k)*ez-p*ez-(1/Mq)*G*w(:,k));
    
    %follower
    X(:,:,k+1+late) = X(:,:,k+late)* expm(h*V(:,:,k+late));
    V(:,:,k+1+late) = V(:,:,k+late) + h*inD*(brak(Jq,V(:,:,k+late)^2))*inD + h*u(:,:,k+late);
    
    x(:,k+1+late) = x(:,k+late)+ h*v(:,k+late);
    v(:,k+1+late) = v(:,k+late)+ h*((alfa(:,k+late)'*(X(:,:,k+late)*ez))*(X(:,:,k+late)*ez)-p*ez-(1/Mq)*G*v(:,k+late));
    
    
    
    %PID R^3
    err2(:,k+1+late)= z(:,k+1)-x(:,k+1+late)-[0;0;5];
    omega2(:,k+1+late)= omega(:,k+late)+ err2(:,k+1+late)*h;
    eps2(:,k+1+late)=w(:,k+1)- v(:,k+1+late);
    %ki2(k+1)=ki_hat2*(omega2(:,k+1)'* eps2(:,k+1));
    ki2(k+1+late) = ki_hat2;
    uc2(:,k+1+late)= (w(:,k+1)-w(:,k))/h+p*ez+(1/Mq)*G*v(:,k+1+late);
    alfa(:,k+1+late)= kp2*err2(:,k+1+late) + ki2(k+1+late)* omega2(:,k+1+late) + kd2*eps2(:,k+1+late)+uc2(:,k+1+late);
    sigma2(k+1+late)= 1/2*norm(alfa(:,k+1+late))^2;

    d_R3(k+1) = dist_R3(z(:,k+1),x(:,k+1+late));
    
    %PID SO(3)
    %err(:,:,k+1)= realLogSO(X(:,:,k+1)'*Z(:,:,k+1)) - 0.00001*vee(alfa(:,k));
    %err(:,:,k+1)= realLogSO(X(:,:,k+1)'*Z(:,:,k+1)) - 0.00001*vee(cross(z(:,k+1)-x(:,k+1)-[0;0;5], alfa(:,k)));
    %err(:,:,k+1)= 0.01*vee(cross(z(:,k+1)-x(:,k+1)-[0;0;5], X(:,:,k+1)*ez));
    m(:,:,k+1+late) = 0.001*vee(cross(X(:,:,k+1+late)*ez, alfa(:,k+1+late)));
    err(:,:,k+1+late)= realLogSO(X(:,:,k+1+late)'*Z(:,:,k+1))+ m(:,:,k+1+late);
    
    omega(:,:,k+1+late)=omega(:,:,k+late) + err(:,:,k+1+late)*h;
    eps(:,:,k+1+late)=W(:,:,k+1)- V(:,:,k+1+late)+ (m(:,:,k+1+late)-m(:,:,k+late))/h;
    ki(k+1+late)=ki_hat*trace(omega(:,:,k+1+late)'* eps(:,:,k+1+late));
    uc(:,:,k+1+late)= (W(:,:,k+1)-W(:,:,k))/h...
                 -inD*brak(Jq,V(:,:,k+late)^2)*inD;
             
    u(:,:,k+1+late)= kp*err(:,:,k+1+late) + ki(k+1+late)* omega(:,:,k+1+late) + kd*eps(:,:,k+1+late) + uc(:,:,k+1+late);
    sigma(k+1+late)= 1/2*norm(u(:,:,k+1+late),'fro')^2;
    sigma_uc(k+1+late)= 1/2*norm(uc(:,:,k+1+late),'fro')^2;

    d_SO3(k+1) = dist_SO3(Z(:,:,k+1),X(:,:,k+1+late));
    
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
axis([-5 50 -5 50 -180 5])
FS = 11; % Font size
set(gca,'FontSize',FS);xlabel('$$x_1$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$x_2$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$x_3$$','interpreter','latex');
hold on

%Video
v=VideoWriter('Drone_nomix.mov');
open(v);

for k=1:10:n
    delete(hh);
    DroniRuotati_leader = Z(:,:,k+1)*(Droni - Centrorot) + z(:,k+1);        
    DroniRuotati_follower = X(:,:,k+1)*(Droni - Centrorot) + x(:,k+1);               
    hh=plot3(DroniRuotati_leader(1,:),DroniRuotati_leader(2,:),DroniRuotati_leader(3,:),'ro-',...
             DroniRuotati_follower(1,:),DroniRuotati_follower(2,:),DroniRuotati_follower(3,:),'bo-');
    drawnow;
    %video
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

%grafici della distanza e sforzo di controllo
t = 0:h:n*h; 
t1 = 0:h:(n+late)*h;
FS = 11; % Font size
figure('Name','Quadcopter drone');
subplot(2,2,1); semilogy(t,d_SO3,'b-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(Z,X)$$','interpreter','latex');
subplot(2,2,2); semilogy(t1,sigma,'b-',t1,sigma_uc,'r-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Control effort $$\sigma_{SO(3)}$$','interpreter','latex');
subplot(2,2,3); semilogy(t,d_R3,'b-',t,5*ones(size(d_R3)),'k-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Distance $$d(z,x)$$','interpreter','latex');
subplot(2,2,4); semilogy(t1,sigma2,'b-'); axis tight
set(gca,'FontSize',FS);xlabel('Time $$t$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('Control effort $$\sigma_{(R^3)}$$','interpreter','latex');