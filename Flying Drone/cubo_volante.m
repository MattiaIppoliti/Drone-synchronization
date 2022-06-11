% Programma per la simulazione numerica di un "cubo volante" con un oscillatore di tipo pendolo classico :)
% di I. Cervigni, M. Ippoliti, C. Menotta (UnivPM), 
% Ultima revisione: 18 Aprile 2019
clear all; close all; clc;
% Mappa esponenziale sulla sfera
exphS2 = @(x,v,h) x*cos(h*norm(v)) + h*v*sinc(h*norm(v)/pi); 
% Trasporto parallelo sulla sfera
trasp = @(x,y,w) w - ( y*(y'*w) - x*(x'*y)*(y'*w) ) / (1+x'*y) - x*(y'*w);
% distS2anza Riemanniana
distS2 = @(x,y) abs(acos(x'*y));
% Logarithm
loghS2 = @(x,y) (y - x*(x'*y)) / sinc(distS2(x,y)/pi);

% Logarithm
loghSO3 = @(x,y) x*realLogSO(x'*y);
% Distanza Riemanniana
distSO3 = @(x,y) norm(realLogSO(x'*y),'fro');

%opengl software % Per le trasparenze

Spigoli = [ 0 0.3 0.3 0 0 0 0.3 0.3 0 0 0 0 0 0 0 0.3 0.3 0.3 0.3 
            0 0 0 0 0 0.3 0.3 0 0 0.3 0.3 0 0 0.3 0.3 0.3 0 0.3 0.3 
            0 0 0.3 0.3 0 0 0 0 0 0 0.3 0.3 0 0 0.3 0.3 0.3 0.3 0];
        
Polo = [0.15 0.15 0.15]'; % Il cubo ruoterà intorno a questo punto, può essere sia interno che esterno al cubo


n = 15001; % Numero di punti della simulazione
h = 0.001; % Passo del metodo numerico
e = 1.3; % Esponenziale per il nonlinear damping
kappa = 0.5; % Coefficiente nel potenziale del pendolo
mu = 0;%0.5; % 0.5 per il caso smorzato per avere belle traiettorie % Coefficiente di smorzamento
x = zeros(3,n); v = zeros(3,n); % Preallocazione aumenta efficienza del codice
x(:,1) = [0; 0; 1]; % randn(3,1); x(:,1) = x(:,1)/norm(x(:,1)); % Stato iniziale
rifS2 = [1; 0; 0];  % riferimento per l'oscillatore
v(:,1) = [0.5; -0.9; 0]; % Velocità iniziale
X = zeros(3,3,n); W = zeros(3,3,n); % Preallocazione aclcumenta efficienza del codice
X(:,:,1)=eye(3); % Stato iniziale
rifSO3 = [-1 0 0;0 1 0;0 0 -1]; % Riferimento per l'oscillatore
W(:,:,1) = [0 -0.5 -1; 0.5 0 0.7; 1 -0.7 0]; 


 
% Metodo numerico (Eulero)
for k=1:n
    x(:,k+1)= exphS2(x(:,k),v(:,k),h);
    v(:,k+1)= trasp(x(:,k),x(:,k+1),v(:,k)...
                                    -h*mu*(norm(v(:,k))^(2*(e-1)))*v(:,k)...
                                    +h*kappa*sinc(distS2(x(:,k),rifS2)/pi)*loghS2(x(:,k),rifS2)...
                                    );
                                
    X(:,:,k+1) = X(:,:,k)*expm(h*W(:,:,k));
    W(:,:,k+1) = W(:,:,k) - h*mu*norm(W(:,:,k),'fro')^(2*(e-1))*W(:,:,k)...
                        + h*kappa*sinc(distSO3(X(:,:,k),rifSO3)/pi)*realLogSO(X(:,:,k)'*rifSO3);
           
end
 
%La rappresentazione del cubo viene messa in un altro ciclo, poichè non si
%prendono tutti i campioni ma se ne prende solo 1 ogni 20, in maniera da rendere la
%rappresentazione più veloce
figure;
sphere(50);alpha(0.1); % Sfera trasparente
grid on; 
FS = 15; % Font size
set(gca,'FontSize',FS);xlabel('$$x$$','interpreter','latex');
set(gca,'FontSize',FS);ylabel('$$y$$','interpreter','latex');
set(gca,'FontSize',FS);zlabel('$$z$$','interpreter','latex');
hold on;
axis([-1.3 1.3 -1.3 1.3 -1.3 1.3]) 
plot3(rifS2(1),rifS2(2),rifS2(3),'bo','LineWidth',2);
SpigoliTraslati= Spigoli - Polo;
SpigoliRuotati = X(:,:,1)*SpigoliTraslati + x(:,1);
hh = plot3(SpigoliRuotati(1,:),SpigoliRuotati(2,:),SpigoliRuotati(3,:),'ro-',x(1,1),x(2,1),x(3,1),'b.');
for k=1:20:n
    delete(hh);
    SpigoliTraslati= Spigoli - Polo;
    SpigoliRuotati = X(:,:,k)*SpigoliTraslati + x(:,k);
    hh = plot3(SpigoliRuotati(1,:),SpigoliRuotati(2,:),SpigoliRuotati(3,:),'ro-',x(1,1:k),x(2,1:k),x(3,1:k),'b.');
    drawnow; 
    %pause(0.01)
end
