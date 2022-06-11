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

%opengl software % Per le trasparenze
n = 100000; % Numero di punti della simulazione
h=0.001; % passo di campionamento

%costanti del drone 
Jx=7.5*10^(-3); Jy=Jx; Jz=1.3*10^(-2);
Jq=(1/2)*diag([Jy-Jx+Jz Jx-Jy+Jz Jx+Jy-Jz]);
D=diag([sqrt((Jy*Jz)/Jx) sqrt((Jx*Jz)/Jy) sqrt((Jx*Jy)/Jz)]);
inD= inv(D);

b= 3.13*10^(-5);
r=0.23;
y=7.15*10^(-7);
p=9.81;
Mq=0.650;

nh=3048;

w1=10.05; w2=10; w3=10; w4=10;

% Preallocazione aumenta efficienza del codice
R = zeros(3,3,n);
W = zeros(3,3,n); 
R(:,:,1)=eye(3); % Stato iniziale

tau=zeros(3,3,n); % Preallocazione aumenta efficienza del codice

Wx=[0 0 0; 0 0 -1;0 1 0];
Wy=[0 0 1; 0 0 0; -1 0 0];
Wz=[ 0 -1 0; 1 0 0; 0 0 0];
tau=b*r*(w4^2-w2^2)*Wx+b*r*(w3^2-w1^2)*Wy+y*(-w1^2+w2^2-w3^2+w4^2)*Wz;
Wr=-w1+w2-w3+w4;


% Metodo numerico (Eulero)
disp('Simulazione numerica...')
for k=1:n
   
    R(:,:,k+1) = R(:,:,k)* expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) + h*inD*(brak(Jq,W(:,:,k)^2)+tau)*inD;
    
end

figure;
disp('Animazione grafica 3D...')

%La rappresentazione del droneo viene messa in un altro ciclo, poichè non si
%prendono tutti i campioni ma se ne prende solo 1 ogni 100, in maniera da rendere la
%rappresentazione più veloce
DroniRuotati = R(:,:,1)*(Droni - Centrorot) + Centrorot;        
hh=plot3(DroniRuotati(1,:),DroniRuotati(2,:),DroniRuotati(3,:),'ro-');
grid on;
axis([-3 3 -3 3 -3 3])
FS = 16; % Font size
set(gca,'FontSize',FS);xlabel('$$x$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$y$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$z$$','interpreter','latex');
hold on

for k=1:100:n
    delete(hh);
    DroniRuotati = R(:,:,k+1)*(Droni - Centrorot) + Centrorot;        
    hh=plot3(DroniRuotati(1,:),DroniRuotati(2,:),DroniRuotati(3,:),'ro-');
    drawnow;
end
