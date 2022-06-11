function L = realLogSO(R)
% Funzione per il calcolo del logaritmo reale di una matrice in SO(n)
% basata sull'articolo di N. Sherif e E. Morsy (Computing real logarithm of
% a real matrix) di Simone Fiori (DII-UnivPM, Ottobre 2018)
[V,T] = schur(R);
% R è una matrice di rotazione quindi se è di ordine dispari ha un
% autovalore pari a 1. Essendo una matrice "normale", la sua decomposizione 
% di Schur è diagonale a blocchi 1x1 e 2x2.
D = []; % Inizializza D ad una matrice vuota
[n,~]=size(R);
if mod(n,2) % Se n è dispari, T contiene un blocco "1" che occorre trovare
    precisione = 10^(-10);
    [~,pos] = min(abs(diag(T)-1-precisione)); % Posizione stimata del blocco "1": 
    % il valore non è esattamente 1 quindi si cerca il valore più vicino;
    % una precisione non adeguata potrebbe causare il malfunzionamento della
    % funzione
    T(pos,:)=[]; T(:,pos)=[]; % Elimina la riga e la colonna pos-ma
    app = V(:,pos); V(:,pos)=[]; V=[app V]; % Porta la colonna pos-ma in prima posizione
    n = n-1; % Riduce la dimensione di 1
    D(1,1) = 0; % Logaritmo del blocco 1 = 0
end
for b = 1:n/2
    B = T(2*b-1:2*b,2*b-1:2*b); % Estrae i blocchi 2x2
    ll = log(B(1,1) + B(1,2)*1i); % real(ll) = 0 per una matrice di rotazione
    [s,~] = size(D);
    D(s+1:s+2,s+1:s+2) = [0 imag(ll);-imag(ll) 0]; % Aggiunge un blocco-logaritmo a D
end
L = V*D*V'; % Risultato finale
return
