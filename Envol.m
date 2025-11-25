function [Y, H, w] = Envolvente(y, WS, L_FFT, p, r)

    y = filter([1; -0.625],1,y); 
    y = y.*hamming(WS);

    for k = 0 : p                        
        r(k+1) = y((k+1):WS)'*y(1:(WS-k)); 
    end

    R = r(1:p);                        
    r = r(2:(p+1)); 

    for k1 = 2 : p                         
        R(:,k1) = [R(k1,1); R(1:p-1,k1-1)];
    end

    a = -inv(R)*r;

    Y = fft(y,L_FFT);
    Y = Y(1:(L_FFT/2));

    [H, w] = freqz(1, [1; a]', (L_FFT/2));

    H = H * (max(abs(Y)))/(max(abs(H)));  

    Y = 20*log10(abs(Y));   
    H = 20*log10(abs(H));  

end