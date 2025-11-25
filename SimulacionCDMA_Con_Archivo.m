
channels=4;
xn=(audioread('E.wav'))';
for i=1:channels
    m(i,:)=xn;
end
figure, plot(m(1,:),'LineWidth',3)

% Parametro A-law
A=87.6;
% Niveles de cuantizacion
n=4;

xnew=zeros(channels,length(xn));
for j=1:channels
   for i=1:length(xn)
       if m(j,i)==0
          xnew(j,i)=0.01;
       else xnew(j,i)=m(j,i);
       end
   end
end
m=xnew;

%%%%%%%%
%Ley-A
for i=1:channels
    y=((1+log(A*abs(m(i,:))))./(1+log(A))).*sign(m(i,:));
    y_quan(i,:)=y;
    b_quan=y_quan;
end
figure, plot(y_quan(1,:))

d=2/n;
%Distribución de los niveles (0,2)
q=d.*[0:n-1];
%Distribución de los niveles (-1,1)
q=q-((n-1)/2)*d;


%Barrido para los valores de salida
for i=1:n
  y_quan(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2)))=...
  q(i).*ones(1,length(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2))));
  b_quan(find(y_quan==q(i)))=(i-1).*ones(1,length(find(y_quan==q(i))));
end
figure, plot(y_quan(1,:))
figure, plot(b_quan(1,:))

% Para Valores de palabras de codigo
bits=ceil(log2(n));

code1=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(1,i)/(2^j)) == 1)
	code1(i,(bits-j)) = 1;
	b_quan(1,i) = b_quan(1,i) - 2^j;
          end   
    end    
end

code2=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(2,i)/(2^j)) == 1)
	code2(i,(bits-j)) = 1;
	b_quan(2,i) = b_quan(2,i) - 2^j;
          end   
    end    
end
code3=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(3,i)/(2^j)) == 1)
	code3(i,(bits-j)) = 1;
	b_quan(3,i) = b_quan(3,i) - 2^j;
          end   
    end    
end
code4=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(4,i)/(2^j)) == 1)
	code4(i,(bits-j)) = 1;
	b_quan(4,i) = b_quan(4,i) - 2^j;
          end   
    end    
end

%Paralelo a Serie
Data1=reshape(code1,[length(code1)*bits,1]);
figure, plot(Data1)
Data1(find(Data1 == 0))=-1;
figure, plot(Data1)
Data2=reshape(code2,[length(code2)*bits,1]);
Data2(find(Data2 == 0))=-1;
Data3=reshape(code3,[length(code3)*bits,1]);
Data3(find(Data3 == 0))=-1;
Data4=reshape(code4,[length(code4)*bits,1]);
Data4(find(Data4 == 0))=-1;

%Generación de matriz de transmisión
T = [(Data1)';
     (Data2)';
     (Data3)';
     (Data4)'];

figure, plot(T(1,:))

%Codigos únicos para cada Tx (Conjunto de Walsh-Hadamard)
H = [ 1   1  1  1 ;
      1  -1  1 -1 ;
      1   1 -1 -1 ;
      1  -1 -1  1 ];

H=hadamard(4);

nb = length(H);       %numero de bits
Y = size(T);              %Tamaño matriz Tx
N = Y(1);                  % Flujo Tx/bit
nf = Y(2);                   % numero de bits/flujo

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

%Valores transmitidos
Tx=reshape(G,1,N*nf);
figure, plot(Tx,'LineWidth',3)

% Señal CDMA a traves canal AWGN
ch=awgn(Tx,9);
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
            R(i,m) = chx(i,m) * H (n,m);
            TOT(i) = TOT(i) + R (i,m);
         end
    end
    Rx = [Rx;TOT/nb];
end

figure, plot(T(1,:),'LineWidth',1)
hold on
plot(round(Rx(1,:)),'r','LineWidth',1)

Rxx=round(Rx);
Rxx(find(Rxx == -1))=0;
Rx=Rxx;

DataRx1=reshape(Rx(1,:),[length(Rx)/bits,bits]);
DataRx2=reshape(Rx(2,:),[length(Rx)/bits,bits]);
DataRx3=reshape(Rx(3,:),[length(Rx)/bits,bits]);
DataRx4=reshape(Rx(4,:),[length(Rx)/bits,bits]);

%Valores de amplitud a partir de la codificacion
[row,col]=size(DataRx1);
for i=1:length(DataRx1)
    mag1=0;
    k=col;
    for j=k:-1:1
       mag1=mag1+DataRx1(i,j)*2^(k-j);        
    end
    amp1(i)=mag1;   
end
figure, plot(amp1)

[row,col]=size(DataRx2);
for i=1:length(DataRx2)
    mag2=0;
    k=col;
    for j=k:-1:1
       mag2=mag2+DataRx2(i,j)*2^(k-j);        
    end
    amp2(i)=mag2;   
end
figure, plot(amp2)

[row,col]=size(DataRx3);
for i=1:length(DataRx3)
    mag3=0;
    k=col;
    for j=k:-1:1
       mag3=mag3+DataRx3(i,j)*2^(k-j);        
    end
    amp3(i)=mag3;   
end
figure, plot(amp3)

[row,col]=size(DataRx4);
for i=1:length(DataRx4)
    mag4=0;
    k=col;
    for j=k:-1:1
       mag4=mag4+DataRx4(i,j)*2^(k-j);        
    end
    amp4(i)=mag4;   
end
figure, plot(amp4)

%%%%%%%%%%%%
%Valores de Amplitud [-1,+1]
val1=amp1-mean(amp1);
vmax1=max(abs(val1));
val_quan1=val1/vmax1;
figure, plot(val_quan1)

val2=amp1-mean(amp2);
vmax2=max(abs(val2));
val_quan2=val2/vmax2;
figure, plot(val_quan2)

val3=amp3-mean(amp3);
vmax3=max(abs(val3));
val_quan3=val3/vmax3;
figure, plot(val_quan3)

val4=amp4-mean(amp4);
vmax4=max(abs(val4));
val_quan4=val4/vmax4;
figure, plot(val_quan4)

%%%%%%%%%%%%%%%%%%
%Para proceso inverso
x_new1=((exp((abs(val_quan1)*(1+log(A)))-1))./A).*sign(val_quan1);
figure, plot(x_new1)
xnew_max1=max(abs(x_new1));
%Normalizar [-1,+1]
xnew_norm1=x_new1/xnew_max1;
figure, plot(xnew_norm1)

x_new2=((exp((abs(val_quan2)*(1+log(A)))-1))./A).*sign(val_quan2);
figure, plot(x_new2)
xnew_max2=max(abs(x_new2));
%Normalizar [-1,+1]
xnew_norm2=x_new2/xnew_max2;
figure, plot(xnew_norm2)

x_new3=((exp((abs(val_quan3)*(1+log(A)))-1))./A).*sign(val_quan3);
figure, plot(x_new3)
xnew_max3=max(abs(x_new3));
%Normalizar [-1,+1]
xnew_norm3=x_new3/xnew_max3;
figure, plot(xnew_norm3)

x_new4=((exp((abs(val_quan4)*(1+log(A)))-1))./A).*sign(val_quan4);
figure, plot(x_new4)
xnew_max4=max(abs(x_new4));
%Normalizar [-1,+1]
xnew_norm4=x_new4/xnew_max4;
figure, plot(xnew_norm4)

%%
%Simulacion CDMA, Info N, TX
clear;
clc;
channels=4;
[xn, fs] = audioread('E.wav');
plot(xn);
title('Señal de voz Original Mensaje de Voz')
sound(xn);
for i=1:channels
    m(i,:)=xn;
end
figure, plot(m(1,:),'LineWidth',3)

% Parametro A-law
A=87.6;
% Niveles de cuantizacion
n=4;

xnew=zeros(channels,length(xn));
for j=1:channels
   for i=1:length(xn)
       if m(j,i)==0
          xnew(j,i)=0.01;
       else xnew(j,i)=m(j,i);
       end
   end
end
m=xnew;

%%%%%%%%
%Ley-A
for i=1:channels
    y=((1+log(A*abs(m(i,:))))./(1+log(A))).*sign(m(i,:));
    y_quan(i,:)=y;
    b_quan=y_quan;
end
figure, plot(y_quan(1,:))

d=2/n;
%Distribución de los niveles (0,2)
q=d.*[0:n-1];
%Distribución de los niveles (-1,1)
q=q-((n-1)/2)*d;


%Barrido para los valores de salida
for i=1:n
  y_quan(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2)))=...
  q(i).*ones(1,length(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2))));
  b_quan(find(y_quan==q(i)))=(i-1).*ones(1,length(find(y_quan==q(i))));
end
figure, plot(y_quan(1,:))
figure, plot(b_quan(1,:))

% Para Valores de palabras de codigo
bits=ceil(log2(n));

code1=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(1,i)/(2^j)) == 1)
	code1(i,(bits-j)) = 1;
	b_quan(1,i) = b_quan(1,i) - 2^j;
          end   
    end    
end

code2=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(2,i)/(2^j)) == 1)
	code2(i,(bits-j)) = 1;
	b_quan(2,i) = b_quan(2,i) - 2^j;
          end   
    end    
end
code3=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(3,i)/(2^j)) == 1)
	code3(i,(bits-j)) = 1;
	b_quan(3,i) = b_quan(3,i) - 2^j;
          end   
    end    
end
code4=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(4,i)/(2^j)) == 1)
	code4(i,(bits-j)) = 1;
	b_quan(4,i) = b_quan(4,i) - 2^j;
          end   
    end    
end

%Paralelo a Serie
Data1=reshape(code1,[length(code1)*bits,1]);
figure, plot(Data1)
Data1(find(Data1 == 0))=-1;
figure, plot(Data1)
Data2=reshape(code2,[length(code2)*bits,1]);
Data2(find(Data2 == 0))=-1;
Data3=reshape(code3,[length(code3)*bits,1]);
Data3(find(Data3 == 0))=-1;
Data4=reshape(code4,[length(code4)*bits,1]);
Data4(find(Data4 == 0))=-1;

%Generación de matriz de transmisión
T = [(Data1)';
     (Data2)';
     (Data3)';
     (Data4)'];

figure, plot(T(1,:))
title('Señal Transmitida del mensaje de audio')

%Codigos únicos para cada Tx (Conjunto de Walsh-Hadamard)
H = [ 1   1  1  1 ;
      1  -1  1 -1 ;
      1   1 -1 -1 ;
      1  -1 -1  1 ];

H=hadamard(4);

nb = length(H);       %numero de bits
Y = size(T);              %Tamaño matriz Tx
N = Y(1);                  % Flujo Tx/bit
nf = Y(2);                   % numero de bits/flujo

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

%Valores transmitidos
Tx=reshape(G,1,N*nf);
figure, plot(Tx,'LineWidth',3)


% Señal CDMA a traves canal AWGN
ch=awgn(Tx,15);
figure, plot(ch,'LineWidth',2)

%Conversion Serie/paralelo
chx=(reshape(ch,nf,N));

% Decodificacion en los Rx
Rx=[];          % vector de bits reconstruidos bits
for n=1:N
    TOT=zeros(1,nf);
    R=zeros(nf,nb);
    for i = 1:nf
        for m = 1:nb
            R(i,m) = chx(i,m) * H (n,m);
            TOT(i) = TOT(i) + R (i,m);
         end
    end
    Rx = [Rx;TOT/nb];
end

figure, plot(T(1,:),'LineWidth',3)
title('Señal Transmitida Mensaje de Audio')
hold on
plot(round(Rx(1,:)),'r','LineWidth',2)
title('Señal Recibida Mensaje de Audio')

Rxx=round(Rx);
Rxx(find(Rxx == -1))=0;
Rx=Rxx;

DataRx1=reshape(Rx(1,:),[length(Rx)/bits,bits]);
DataRx2=reshape(Rx(2,:),[length(Rx)/bits,bits]);
DataRx3=reshape(Rx(3,:),[length(Rx)/bits,bits]);
DataRx4=reshape(Rx(4,:),[length(Rx)/bits,bits]);

%Valores de amplitud a partir de la codificacion
[row,col]=size(DataRx1);
for i=1:length(DataRx1)
    mag1=0;
    k=col;
    for j=k:-1:1
       mag1=mag1+DataRx1(i,j)*2^(k-j);        
    end
    amp1(i)=mag1;   
end
figure, plot(amp1)

[row,col]=size(DataRx2);
for i=1:length(DataRx2)
    mag2=0;
    k=col;
    for j=k:-1:1
       mag2=mag2+DataRx2(i,j)*2^(k-j);        
    end
    amp2(i)=mag2;   
end
figure, plot(amp2)

[row,col]=size(DataRx3);
for i=1:length(DataRx3)
    mag3=0;
    k=col;
    for j=k:-1:1
       mag3=mag3+DataRx3(i,j)*2^(k-j);        
    end
    amp3(i)=mag3;   
end
figure, plot(amp3)

[row,col]=size(DataRx4);
for i=1:length(DataRx4)
    mag4=0;
    k=col;
    for j=k:-1:1
       mag4=mag4+DataRx4(i,j)*2^(k-j);        
    end
    amp4(i)=mag4;   
end
figure, plot(amp4)

%%%%%%%%%%%%
%Valores de Amplitud [-1,+1]
val1=amp1-mean(amp1);
vmax1=max(abs(val1));
val_quan1=val1/vmax1;
figure, plot(val_quan1)

val2=amp1-mean(amp2);
vmax2=max(abs(val2));
val_quan2=val2/vmax2;
figure, plot(val_quan2)

val3=amp3-mean(amp3);
vmax3=max(abs(val3));
val_quan3=val3/vmax3;
figure, plot(val_quan3)

val4=amp4-mean(amp4);
vmax4=max(abs(val4));
val_quan4=val4/vmax4;
figure, plot(val_quan4)

%%%%%%%%%%%%%%%%%%
%Para proceso inverso
x_new1=((exp((abs(val_quan1)*(1+log(A)))-1))./A).*sign(val_quan1);
figure, plot(x_new1)
xnew_max1=max(abs(x_new1));
%Normalizar [-1,+1]
xnew_norm1=x_new1/xnew_max1;
figure, plot(xnew_norm1)

x_new2=((exp((abs(val_quan2)*(1+log(A)))-1))./A).*sign(val_quan2);
figure, plot(x_new2)
xnew_max2=max(abs(x_new2));
%Normalizar [-1,+1]
xnew_norm2=x_new2/xnew_max2;
figure, plot(xnew_norm2)

x_new3=((exp((abs(val_quan3)*(1+log(A)))-1))./A).*sign(val_quan3);
figure, plot(x_new3)
xnew_max3=max(abs(x_new3));
%Normalizar [-1,+1]
xnew_norm3=x_new3/xnew_max3;
figure, plot(xnew_norm3)

x_new4=((exp((abs(val_quan4)*(1+log(A)))-1))./A).*sign(val_quan4);
figure, plot(x_new4)
xnew_max4=max(abs(x_new4));
%Normalizar [-1,+1]
xnew_norm4=x_new4/xnew_max4;
figure, plot(xnew_norm4)
title('Señal final a partir de la multiplexación CDMA')
sound(xnew_norm3);
audiowrite('CDMA_VocalE_15dB.wav',xnew_norm3, fs);