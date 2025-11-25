%Simulacion FDMA
muestras=2000;
% numero de usuarios
channels=8;
% Frecuencia Señal moduladora en Hz
fm=[30 40 50 60 70 80 90 100]; 
% Frecuencia de portadora asignada a los diferentes usuarios en Hz
fc=[300 600 900 1200 1500 1800 2100 2400];
% Desviacion de frecuencia
fd=10;

% Generar la señal moduladora 
t=linspace(0,muestras,muestras);
for i=1:channels
    m(i,:)=sin(2*pi*fm(1,i)*t)+2*sin(pi*8*t);
end
figure, plot(m(1,:),'LineWidth',3)

% Generar Señal modulada
for i=1: channels
    y(i,:)=fmmod(m(i,:),fc(1,i),4*fc(1,i),fd);
end
figure, plot(y(1,:) , 'LineWidth',3)



% Señal modulada a traves canal AWGN 
for i=1:channels
        chn(i,:)=awgn(y(i,:),10,'measured');
end
figure, plot(chn(1,:),'LineWidth',3)



%Parte Recepción
% Demodulacion señal recibida en la estacion base 
for i=1:channels
    Rx(i,:)=fmdemod(chn(i,:),fc(1,i),4*fc(1,i),fd);
end
figure, plot(Rx(1,:),'LineWidth',3)

figure, plot(Rx(1,:),'LineWidth',4)
hold on
plot(m(1,:),'r','LineWidth',2)

