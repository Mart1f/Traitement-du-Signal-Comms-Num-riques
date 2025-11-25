%Codigo de Cuantizacion Uniforme VOCAL U 32
% Valores de entrada
%xn=(audioread('U.wav'))';
[xn, fs_original] = audioread('U.wav');
sound(xn)
plot(xn)
% Parametro mu-law
mu=255;


% Niveles de cuantizacion
n=32;
%%%%%%%%
%Ley-mu
y=(log(1+mu*abs(xn))./log(1+mu)).*sign(xn);
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
bb_quan=b_quan;

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
     if(fix(b_quan(i)/(2^j)) == 1)
	code(i,(bits-j)) = 1;
	b_quan(i) = b_quan(i) - 2^j;
     end   
 end    
end

%%%%%%%%%%%%%%
%Valores de amplitud a partir de la codificacion
[row,col]=size(code);
for i=1:length(code)
    mag=0;
    k=col;
    for j=k:-1:1
       mag=mag+code(i,j)*2^(k-j);        
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
x_new=(((1+mu).^(abs(val_quan))-1)./mu).*sign(y_quan);
figure, plot(x_new)

xnew_max=max(abs(x_new));
%Normalizar [-1,+1]
xnew_norm=x_new/xnew_max;
figure, plot(xnew_norm)

MSE=(sum((xn-x_new).^2))/length(xn)*100

audiowrite('U32LeyMu.wav', xnew_norm, fs_original);

% Reproducir el archivo grabado para verificar
sound(xnew_norm, fs_original);