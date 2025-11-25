function yn = NormalizarSenal(in)
    
    % Máximo absoluto de in
    max_entrada = max(abs(in));

    % Dividir in por su valor maximo -> [0,1]
    yn = in / max_entrada;
end
