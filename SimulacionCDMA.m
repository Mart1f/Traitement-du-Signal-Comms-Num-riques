%Simulacion CDMA, Info N, TX
T = [  1 -1  1 -1  1  1 -1 -1  1; 
      -1 -1  1  1  1 -1 -1  1 -1;
       1  1 -1 -1 -1  1  1 -1  1;
       1  1  1  1 -1 -1 -1 -1  1];

%Codigos únicos para cada Tx (Conjunto de Walsh-Hadamard)
H = [ 1   1  1  1 ;
      1  -1  1 -1 ;
      1   1 -1 -1 ;
      1  -1 -1  1 ];

HH=hadamard(4);

nb = length(H);       %numero de bits
Y = size(T);          %Tamaño matriz Tx
N = Y(1);             % Flujo Tx/bit
nf = Y(2);            % numero de bits/flujo

% bits codificados en los Tx M(nb) I(nf)
G = zeros(nf,nb);
for n = 1:N
    Z = zeros(nf,nb);
    for i = 1:nf
        for m = 1:nb
            Z(i,m) = [T(n,i)*H(n,m)];        
        end
    end
    G = G + Z;
end

%Valores a transmitir
Tx=reshape(G,1,N*nf);
figure, plot(Tx,'LineWidth',3)

% Señal CDMA a traves canal AWGN
ch=awgn(Tx,30);
figure, plot(ch,'LineWidth',3)

%Conversion Serie/paralelo
chx=(reshape(ch,nf,N));

% Decodificacion en los Rx
Rx=[];          % vector de bits reconstruidos bits
for n=1:N
    TOT=zeros(1,nf);
    R=zeros(nf,nb);
    for i = 1:nf
        for m = 1:nb
            R(i,m) = chx(i,m) * H(n,m);
            TOT(i) = TOT(i) + R(i,m);
         end
    end
    Rx = [Rx;TOT/nb];
end

figure, plot(T(3,:),'LineWidth',3)
hold on
plot(round(Rx(3,:)),'r','LineWidth',2)