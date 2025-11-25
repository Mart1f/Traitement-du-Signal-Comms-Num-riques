%Codigo de Cuantizacion No Uniforme
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

%Paralelo a Serie
DataTx=reshape(code,[length(code)*bits,1]);
figure, plot(DataTx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solo por codigo
trellis=poly2trellis(3,{'x2 + 1','x + x2 + 1'});
codedData=convenc(DataTx,trellis);
figure, plot(codedData);
%Para BPSK
%M=2;
%Para QPSK y 4QAM
%M=4;

%Para 16QAM
M=16;

%Modulacion para BPSK y QPSK 
%ModData=pskmod(codedData,M,InputType='bit');

%Modulacion para 4QAM 
ModData=qammod(codedData,M,InputType='bit');

scatterplot(ModData);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proceso Paso por Canal AWGN
%%Canal AWGN por codigo
Channel_Rx=awgn(ModData,11,'measured');
scatterplot(Channel_Rx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proceso de Recepcion
%%Recepcion por codigo
% Demodulacion para BPSK y QPSK
%demodData=pskdemod(Channel_Rx,M,OutputType='bit');

% Demodulacion para 4QAM y 16QAM
demodData=qamdemod(Channel_Rx,M,OutputType='bit');
figure,plot(demodData);

Rx=vitdec(demodData,trellis,1,'cont','hard');
%Rx=abs(out.Rx);
figure, plot(Rx);
%%%%%%%%%%%%%%%%%%%%%&
DataRx=reshape(Rx,[length(Rx)/bits,bits]);

%Valores de amplitud a partir de la codificacion
[row,col]=size(DataRx);
for i=1:length(DataRx)
    mag=0;
    k=col;
    for j=k:-1:1
       mag=mag+DataRx(i,j)*2^(k-j);        
    end
    amp(i)=mag;   
end
figure, plot(amp)

%%%%%%%%%%%%
%Valores de Amplitud [-1,+1]
val=amp-mean(amp);
vmax=max(abs(val));
val_quan=val/vmax;
figure, plot(val_quan)

%%%%%%%%%%%%%%%%%%
%Para proceso inverso
x_new=((exp((abs(val_quan)*(1+log(A)))-1))./A).*sign(val_quan);
figure, plot(x_new)

xnew_max=max(abs(x_new));
%Normalizar [-1,+1]
xnew_norm=x_new/xnew_max;
figure, plot(xnew_norm)

MSE=(sum((xn-xnew_norm).^2))/length(xn)*100
audiowrite('MensajeSConvol_11dB_16QAM.wav', xnew_norm, fs_original);
sound(xnew_norm)