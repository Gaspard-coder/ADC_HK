function y_out = CAN(x, N_bits, PE, DSR)

    % 1. Sous-échantillonnage (Downsampling)
    % On ne garde qu'un échantillon tous les 'DSR' échantillons
    if nargin < 4
        DSR = 1;
    end
    
    if DSR > 1
        x_sampled = x(1:DSR:end);
    else
        x_sampled = x;
    end

    quantizedInput = floor((x_sampled+(PE/2))*2^(N_bits-1)/(PE/2)); % Quantification
    quantizedInput(quantizedInput<0) = 0; % Clipping Down
    quantizedInput(quantizedInput>2^N_bits-1) = 2^N_bits-1; % Clipping Up
    y_out = (quantizedInput-2^(N_bits-1))/2^(N_bits-1)*(PE/2)+(PE/2)/2^N_bits;

end