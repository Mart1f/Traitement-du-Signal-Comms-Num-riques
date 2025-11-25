%%
% Modulación QAM
M=32;
k=log2(M);
n=2400;
datos=randi([0 1],n,1);
data=reshape(datos,length(datos)/k,k);
symbol=bi2de(data);
QAM=qammod(symbol,M);
scatterplot(QAM,1,0,'y.');
snr=3;
Rx=awgn(QAM,snr,'measured');
SRx=scatterplot(Rx,1,0,'r.');
hold on
scatterplot(QAM,1,0,'y*',SRx)
Demod=qamdemod(Rx,M);
bits=de2bi(Demod,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Valores de entrada
xn=(audioread('MARTINBERMUDEZ3.wav'))';
sound(xn)
plot(xn)
% Parametro A-law
A=87.6;
% Niveles de cuantizacion
n=4;

xnew=zeros(1,length(xn));
for i=1:length(xn)
    if xn(i)==0
     xnew(i)=0.01;
    else xnew(i)=xn(i);
    end
end

xn=xnew;
%%%%%%%%
%Ley-A
y=((1+log(A*abs(xn)))./(1+log(A))).*sign(xn);
figure, plot(y)

y_quan=y;
b_quan=y_quan;

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
 
figure, plot(y_quan)
figure, plot(b_quan)

MSE=(sum((xn-y_quan).^2))/length(xn)*100

x_orden=sort(xn);
quan_orden=sort(y_quan);
figure, plot(x_orden,quan_orden,'LineWidth',3)

% Para Valores de palabras de codigo
bits=ceil(log2(n));
code=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(i)/(2^j)) == 1)
	code(i,(bits-j)) = 1;
	b_quan(i) = b_quan(i) - 2^j;
          end   
    end    
end

codeserial=reshape(code,[1,length(code)*bits]);
serial=zeros(1,length(codeserial)+1);
serial(1,2:length(codeserial)+1)=codeserial;

% Termina Tx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comienza Rx
Rx=(out.DataRx);

codeparal=reshape(Rx,[length(Rx)/bits,bits]);

%Valores de amplitud a partir de la codificacion
[row,col]=size(codeparal);
for i=1:length(codeparal)
    mag=0;
    k=col;
    for j=k:-1:1
       mag=mag+codeparal(i,j)*2^(k-j);        
    end
    amp(i)=mag;   
end
figure, plot(amp)

amp_max=max(abs(amp));
new_amp=amp-(amp_max/2);
%Normalizar [-1,+1]
new_amp_max=max(abs(new_amp));
amp_norm=new_amp/new_amp_max;
figure, plot(amp_norm)

%%%%%%%%%%%%%%%%%%
%Para proceso inverso
x_new=((exp((abs(amp_norm)*(1+log(A)))-1))./A).*sign(amp_norm);
audiowrite('16QAM_3dB_LeyA.wav', x_new,fs);
figure, plot(x_new)

sound(x_new)