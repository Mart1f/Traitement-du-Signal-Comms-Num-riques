%%
% Modulación QAM utilizando Ley Mu
% Valores de entrada
[xn, fs] = audioread('MARTINBERMUDEZ3.wav');
sound(xn)
figure, plot(xn)

% Parámetro Mu-law
mu = 255;

% Niveles de cuantización
n = 4;

% Procesamiento inicial para evitar ceros
xnew = zeros(1, length(xn));
for i = 1:length(xn)
    if xn(i) == 0
        xnew(i) = 0.01;
    else
        xnew(i) = xn(i);
    end
end
xn = xnew;

%%%%%%%%
% Ley-Mu
y = (log(1 + mu * abs(xn)) ./ log(1 + mu)) .* sign(xn);
figure, plot(y)

% Cuantización de la señal
y_quan = y;
b_quan = y_quan;

d = 2/n;
% Distribución de los niveles (0,2)
q = d .* [0:n-1];
% Distribución de los niveles (-1,1)
q = q - ((n-1)/2) * d;

% Barrido para los valores de salida
for i = 1:n
    y_quan(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2))) = ...
    q(i) .* ones(1, length(find((q(i)-d/2 <= y_quan) & (y_quan <= q(i)+d/2))));
    b_quan(find(y_quan == q(i))) = (i-1) .* ones(1, length(find(y_quan == q(i))));
end

figure, plot(y_quan)
figure, plot(b_quan)

MSE = (sum((xn - y_quan).^2)) / length(xn) * 100;

x_orden = sort(xn);
quan_orden = sort(y_quan);
figure, plot(x_orden, quan_orden, 'LineWidth', 3)

% Codificación a palabras de código
bits = ceil(log2(n));
code = zeros(length(xn), bits);
for i = 1:length(xn)
    for j = bits:-1:0
        if (fix(b_quan(i)/(2^j)) == 1)
            code(i, (bits-j)) = 1;
            b_quan(i) = b_quan(i) - 2^j;
        end
    end    
end

codeserial = reshape(code, [1, length(code)*bits]);
serial = zeros(1, length(codeserial) + 1);
serial(1, 2:length(codeserial) + 1) = codeserial;

% Termina Tx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comienza Rx

Rx = (out.DataRx);

codeparal = reshape(Rx, [length(Rx)/bits, bits]);

% Valores de amplitud a partir de la codificación
[row, col] = size(codeparal);
for i = 1:length(codeparal)
    mag = 0;
    k = col;
    for j = k:-1:1
        mag = mag + codeparal(i,j) * 2^(k-j);        
    end
    amp(i) = mag;   
end
figure, plot(amp)

amp_max = max(abs(amp));
new_amp = amp - (amp_max/2);
% Normalizar [-1,+1]
new_amp_max = max(abs(new_amp));
amp_norm = new_amp / new_amp_max;
figure, plot(amp_norm)

% Para proceso inverso
x_new = (((1 + mu).^(abs(amp_norm)) - 1) ./ mu) .* sign(amp_norm);
figure, plot(x_new)

audiowrite('16QAM_3dB_LeyMu.wav', x_new, fs);
sound(x_new)
