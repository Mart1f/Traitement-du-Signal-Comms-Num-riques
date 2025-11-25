function AudioRecord(fs, bitsPerSample, numChannels, duration, filename)

    % Crear el objeto de grabación
    recObj = audiorecorder(fs, bitsPerSample, numChannels);

    % Iniciar la grabación
    disp('Comenzando la grabación...');
    record(recObj, duration);
    
    % Pausa para permitir la grabación durante la duración especificada
    pause(duration);

    % Detener la grabación
    stop(recObj);
    disp('Grabación finalizada.');

    % Obtener los datos de audio grabados
    audioData = getaudiodata(recObj);

    % Guardar los datos de audio en un archivo .wav
    audiowrite(filename, audioData, fs);

    disp(['Audio grabado guardado en ' filename ' exitosamente.']);
end
