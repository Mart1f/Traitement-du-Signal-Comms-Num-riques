%Cuantizacion Uniforme vocal U 32 niveles

%Vocal U 32 niveles
%xn=audioread('U.wav');
[xn, fs_original] = audioread('U.wav');
sound(xn)
plot(xn)


%Niveles de cuantizacion
n=32;
a_quan=xn;
b_quan=a_quan;
%Niveles-Regiones
d=2/n;
q=d.*[0:n-1];
q=q-((n-1)/2)*d;

%Barrido valores a cuantizar
for i=1:n
 a_quan(find((q(i)-d/2<=a_quan)&(a_quan<=q(i)+d/2)))=...
 q(i).*ones(1,length(find((q(i)-d/2<=a_quan)&(a_quan<=q(i)+d/2))));
 b_quan(find(a_quan==q(i)))=(i-1).*ones(1,length(find(a_quan==q(i))));
end
figure, plot(a_quan)
figure, plot(b_quan)
sound(a_quan)
bb_quan=b_quan;

%error cuadratico medio
MSE=((sum((xn-a_quan).^2))/length(xn))*100

%Ordenar valores
xn_orden=sort(xn);
quan_orden=sort(a_quan);
figure, plot(xn_orden,quan_orden,'LineWidth',3)

%Palabras de Codigo
bits=ceil(log2(n));
code=zeros(length(xn),bits);
for i=1:length(xn)
 for j=bits:-1:0
     if(fix(b_quan(i)/(2^j))==1)
       code(i,(bits-j))=1;
       b_quan(i)=b_quan(i)-2^j;
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

audiowrite('U32Unif.wav', xnew_norm, fs_original);

% Reproducir el archivo grabado para verificar
sound(xnew_norm, fs_original);