% Programma per la simulazione numerica di un drone, capace di fare movimenti rotazionali e traslazionali,
% implementato sul Manifold delle iperrotazioni SO(3)(per le rotazioni)
% e su R^3 (per le traslazioni) di I. Cervigni, M. Ippoliti, C. Menotta (UnivPM), 
% Ultima revisione: 24/04/2019
clear all; close all; clc; 

%caricamento figura drone
load drone
Droni = [Xr';Yr';Zr'];  %le coordinate dei punti che costituiscono la figura del drone
close all
Centrorot = [-0.0125  -0.0048 0]'; %centro di massa del nosro drone


%Lie bracket
brak = @(x,y) x*y-y*x;
%matrix anti-commutator
anticom = @(x,y) x*y+y*x;

%opengl software % Per le trasparenze
n = 5000; % Numero di punti della simulazione
h=0.001; % passo di campionamento

%costanti del drone 
Jx=7.5*10^(-3); Jy=Jx; Jz=1.3*10^(-2);
Jq=(1/2)*diag([Jy-Jx+Jz Jx-Jy+Jz Jx+Jy-Jz]);
D=diag([sqrt((Jy*Jz)/Jx) sqrt((Jx*Jz)/Jy) sqrt((Jx*Jy)/Jz)]);
inD= inv(D);

b= 3.13*10^(-5);
r=0.23;
y=7.5*10^(-7);
p=9.81;
Mq=0.650;

nh=3048;

ez = [0;0;1];

R = zeros(3,3,n); W = zeros(3,3,n); % Preallocazione aumenta efficienza del codice
q = zeros(3,n); v = zeros(3,n); % Preallocazione aumenta efficienza del codice, Q0= origine degli assi
R(:,:,1)=eye(3); % Stato iniziale

Wx=[0 0 0; 0 0 -1;0 1 0];
Wy=[0 0 1; 0 0 0; -1 0 0];
Wz=[ 0 -1 0; 1 0 0; 0 0 0];

%velocità dei 4 motori
%per avere il drone fermo tutte le velocità devono essere 320

wh=3048/30*pi;

choise=menu('Choose the maneuver to perform','Hover','Yaw & Roll(-y)', 'Yaw & Pitch(+x)', 'Pitch (+x) & Roll (+y)');

switch(choise)
    case 1
        RPM1=3048; RPM2=3048; RPM3=3048; RPM4=3048;
    case 2
        RPM1=3048; RPM2=3047; RPM3=3048; RPM4=3049;
    case 3
        RPM1=3047; RPM2=3048; RPM3=3049; RPM4=3048;
    case 4
        RPM1=3047; RPM2=3049; RPM3=3049; RPM4=3047;
end

w1=RPM1*pi/30; w2=RPM2*pi/30; w3=RPM3*pi/30; w4=RPM4*pi/30;

tau=b*r*(w4^2-w2^2)*Wx+b*r*(w3^2-w1^2)*Wy+y*(-w1^2+w2^2-w3^2+w4^2)*Wz;
Wr=-w1+w2-w3+w4;

G=1/4 * eye(3);


% Metodo numerico (Eulero)
disp('Simulazione numerica...')
for k=1:n
    R(:,:,k+1) = R(:,:,k)* expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) + h*inD*(brak(Jq,W(:,:,k)^2)+tau)*inD;
    
    q(:,k+1) = q(:,k)+ h*v(:,k);
    v(:,k+1) = v(:,k)+ h*((1/2)*(b/Mq)*(w1^2+w2^2+w3^2+w4^2)*R(:,:,k)*ez-p*ez-(1/Mq)*G*v(:,k));
    
end


disp('Animazione grafica 3D...')

%La rappresentazione del droneo viene messa in un altro ciclo, poichè non si
%prendono tutti i campioni ma se ne prende solo 1 ogni 100, in maniera da rendere la
%rappresentazione più veloce
figure;

DroniRuotati = R(:,:,1)*(Droni - Centrorot) + q(:,1); 
hh =  plot3(DroniRuotati(1,:),DroniRuotati(2,:),DroniRuotati(3,:),'ro-',q(1,1),q(2,1),q(3,1),'b.');
grid on; 
axis([-30 30 -30 30 -30 30]) 
%limito la finestra grafica alla porzione che mi interessa per vedere bene il drone
FS = 15; % Font size 
set(gca,'FontSize',FS);xlabel('$$x$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$y$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$z$$','interpreter','latex');
hold on

%Video
%v=VideoWriter('Drone.avi', 'Uncompressed AVI');
%open(v);

for k=1:10:n
    delete(hh);
    DroniRuotati = R(:,:,k+1)*(Droni - Centrorot) + q(:,k);        
    hh = plot3(DroniRuotati(1,:),DroniRuotati(2,:),DroniRuotati(3,:),'ro-',q(1,1:k),q(2,1:k),q(3,1:k),'b.');
    drawnow
    %frame = getframe(gcf);
    %writeVideo(v,frame);
end

%close(v);