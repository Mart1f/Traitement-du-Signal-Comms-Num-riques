%OFDM para QPSK y 4QAM
clc;
clear;

%Cargar archivos de voz
[ch1, fs] = audioread('A.wav');
[ch2, ~] = audioread('E.wav');
[ch3, ~] = audioread('U.wav');
[ch4, ~] = audioread('MARTINBERMUDEZ3.wav');

% Ahora sí, realiza la transposición si es necesario
ch1 = ch1';
ch2 = ch2';
ch3 = ch3';
ch4 = ch4';
%Unir variables
xn=cat(2,ch1,ch2,ch3,ch4);

figure;
plot(xn);
title('Mesajes Originales');
xlabel('Muestra');
ylabel('Amplitud');
grid on;

%% proceso de compresion y cuantizacion
A=87.6; %constante Ley-A
n=4;% Niveles de cuantizacion

xnew=zeros(1,length(xn));      %Reemplazo de ceros
for i=1:length(xn)
    if xn(i)==0
     xnew(i)=0.001;
    else xnew(i)=xn(i);
    end
end
xn=xnew;

y=((1+log(A*abs(xn)))./(1+log(A))).*sign(xn); %ecuacion ley-A

figure;
plot(y);
title('Ley-A');
xlabel('Muestra');
ylabel('Amplitud');

y_quan=y;
b_quan=y_quan;

d=2/n;
q=d.*[0:n-1];      %Distribución de los niveles (0,2)
q=q-((n-1)/2)*d;   %Distribución de los niveles (-1,1)
%Barrido para los valores de salida
for i=1:n
  y_quan(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2)))=...
  q(i).*ones(1,length(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2))));
  b_quan(find(y_quan==q(i)))=(i-1).*ones(1,length(find(y_quan==q(i))));
end

figure,
plot(y_quan)
title('y_  quan')
xlabel('Muestra');
ylabel('Bit');
%sound(y_quan)
figure,
plot(b_quan)
title('b_  quan')
xlabel('Muestra');
ylabel('Bit');

x_orden=sort(xn);
quan_orden=sort(y_quan);

bits=ceil(log2(n));            % Para Valores de palabras de codigo
code=zeros(length(xn),bits);
for i=1:length(xn)
   for j=bits:-1:0
        if ( fix(b_quan(i)/(2^j)) == 1)
	code(i,(bits-j)) = 1;
	b_quan(i) = b_quan(i) - 2^j;
          end   
    end    
end

%% Transmision OFDM
%%%Para QPSK y 4QAM
data=code;
% Modulacion
mod_data=pskmod(data',4,InputType='bit'); %QPSK
%mod_data=qammod(data',4,InputType='bit'); %4QAM

data_bits=length(mod_data); % Num. de bits 
M=4; % Num. de subportadoras
bloque=data_bits/M; % bloque OFDM 
Pre_len=floor(0.1*bloque); %Largo del cyclic prefix

figure,
stem(y_quan); 
grid on; 
xlabel('Puntos de Datos'); 
ylabel('Amplitud')
title('Datos Originales');

figure,
stem(abs(mod_data));
title('Modulación QPSK/4QAM');
grid on;
scatterplot(mod_data);%constelacion

Flujo=reshape(mod_data,data_bits/M,M);
Sc1 = Flujo(:,1);
Sc2 = Flujo(:,2);
Sc3 = Flujo(:,3);
Sc4 = Flujo(:,4);

figure, 
subplot(4,1,1),stem(real(Sc1))%,title('Sub-portadora_1'),grid on;
title('SUBPORTADORAS ANTES DEL CANAL')
subplot(4,1,2),stem(real(Sc2))%,title('Sub-portadora_2'),grid on;
subplot(4,1,3),stem(real(Sc3))%,title('Sub-portadora_3'),grid on;
subplot(4,1,4),stem(real(Sc4))%,title('Sub-portadora_4'),grid on;

% IFFT de las subportdoras
n_Sc=4;
Pre_start=bloque-Pre_len;
f_Sc1 = ifft(Sc1);
f_Sc2 = ifft(Sc2);
f_Sc3 = ifft(Sc3);
f_Sc4 = ifft(Sc4);

figure,
subplot(4,1,1),stem(real(f_Sc1),'r'),
title('IFFT de las Sub-Portadores')
subplot(4,1,2),stem(real(f_Sc2),'c')
subplot(4,1,3),stem(real(f_Sc3),'b')
subplot(4,1,4),stem(real(f_Sc4),'g')

% Adicionar el CYCLIC PREFIX 
for i=1:n_Sc,
      ifft_Sc(:,i)=ifft((Flujo(:,i)),bloque);
      for j=1:Pre_len,
          cyclic_prefix(j,i)=ifft_Sc(j+Pre_start,i);
      end
      Append_prefix(:,i)=vertcat(cyclic_prefix(:,i), ifft_Sc(:,i));
end
A1=Append_prefix(:,1);
A2=Append_prefix(:,2);
A3=Append_prefix(:,3);
A4=Append_prefix(:,4);

figure,
subplot(4,1,1),plot(real(A1),'r','LineWidth',2),
title('Prefijo Cíclico + Sub-portadoras')
subplot(4,1,2),plot(real(A2),'c','LineWidth',2)
subplot(4,1,3),plot(real(A3),'b','LineWidth',2)
subplot(4,1,4),plot(real(A4),'g','LineWidth',2)

figure,
plot((real(A1)),'r','LineWidth',2),
title('Superposicion de Espectros de los 4 canales'),
hold on ,
plot((real(A2)),'c','LineWidth',2),
hold on ,
plot((real(A3)),'b','LineWidth',2),
hold on ,
plot((real(A4)),'g','LineWidth',2),
hold on,
grid on 

% Conversion paralelo a serie para Tx
[filas_prefix cols_prefix]=size(Append_prefix);
L_ofdm = filas_prefix*cols_prefix;
% Señal OFDM para ser Tx
ofdm_signal = reshape(Append_prefix, 1, L_ofdm);

figure,
plot(real(ofdm_signal),'LineWidth',2);
xlabel('Tiempo');
ylabel('Amplitud');
title('Señal OFDM antes del canal AWGN');
grid on;

% Señal pasa por canal AWGN
rx_signal=awgn(ofdm_signal,50);

figure,
plot(real(rx_signal),'LineWidth',2);
xlabel('Tiempo');
ylabel('Amplitud');
title('Señal OFDM despúes del Canal AWGN');grid on;

%INICIA BLOQUE RECEPCION 
% Parte de RX OFDM
Rx_Ch = reshape(rx_signal,filas_prefix, cols_prefix);

% Remover el cyclic Prefix
Rx_Ch(1:Pre_len,:)=[];
R1= Rx_Ch (:,1);
R2= Rx_Ch (:,2);
R3= Rx_Ch (:,3);
R4= Rx_Ch (:,4);

figure,
subplot(4,1,1),plot(real(R1),'r','LineWidth',2),
title('Prefijo Cíclico removido de Sub-Portadoras')
subplot(4,1,2),plot(real(R2),'c','LineWidth',2)
subplot(4,1,3),plot(real(R3),'b','LineWidth',2)
subplot(4,1,4),plot(real(R4),'g','LineWidth',2)

% Realizar la FFT a señal Rx
for i=1:n_Sc,
      fft_data(:,i)=fft(Rx_Ch(:,i),bloque);
end
F1=fft_data(:,1);
F2=fft_data(:,2);
F3=fft_data(:,3);
F4=fft_data(:,4);

figure,
subplot(4,1,1),stem(real(F1),'r'),
title('SUBPORTADORAS DESPUES DEL CANAL')
subplot(4,1,2),stem(real(F2),'c')
subplot(4,1,3),stem(real(F3),'b')
subplot(4,1,4),stem(real(F4),'g')

% Señal Reconstruida
% Conversion a serial y demodulación
rx_serial = reshape(fft_data, 1,(bloque*M));
scatterplot(rx_serial);

demod_data = pskdemod(rx_serial,4,OutputType='bit');%QPSK
%demod_data = qamdemod(rx_serial,4,OutputType='bit');%4QAM

%% Proceso inverso ley A y Grabacion de los mensajes modificados
codeparal=(demod_data)';

%%%%%%%%%%%%%%
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

amp_max=max(abs(amp));
new_amp=amp-(amp_max/2);
%Normalizar [-1,+1]
new_amp_max=max(abs(new_amp));
amp_norm=new_amp/new_amp_max;
%%%%%%%%%%%%%%%%%%
%Para proceso inverso ley A
x_new=((exp((abs(amp_norm)*(1+log(A)))-1))./A).*sign(amp_norm);

figure,
plot(x_new);
title('Mesajes Recuperado');
xlabel('Muestra');
ylabel('Amplitud');
grid on;

%separar canales modificados
msg1=x_new(1:8000);
msg2=x_new(8001:16000);
msg3=x_new(16001:24000);
msg4=x_new(24001:48000);

sound(msg1);
pause(2);
sound(msg2);
pause(2);
sound(msg3);
pause(3);
sound(msg4);


audiowrite('OFDM_QPSK_LetraA.wav', msg1, fs);
audiowrite('OFDM_QPSK_LetraE.wav', msg2, fs);
audiowrite('OFDM_QPSK_LetraU.wav', msg3, fs);
audiowrite('OFDM_QPSK_MensajeVoz.wav', msg4, fs);

% audiowrite('OFDM_4QAM_LetraA.wav', msg1, fs);
% audiowrite('OFDM_4QAM_LetraE.wav', msg2, fs);
% audiowrite('OFDM_4QAM_LetraU.wav', msg3, fs);
% audiowrite('OFDM_4QAM_MensajeVoz.wav', msg4, fs);



%MSE
MSE1=(sum((ch1-msg1).^2))/length(ch1)*100
MSE2=(sum((ch2-msg2).^2))/length(ch2)*100
MSE3=(sum((ch3-msg3).^2))/length(ch3)*100
MSE4=(sum((ch4-msg4).^2))/length(ch4)*100