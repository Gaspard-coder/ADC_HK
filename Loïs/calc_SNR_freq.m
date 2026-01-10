function SNR = calc_SNR_freq(x, Fe, f_min, f_max, Fsig)
% Entrées :
%   x     : Signal temporel (sortie du CAN)
%   Fe    : Fréquence d'échantillonnage
%   f_min : Début de la bande d'intégration
%   f_max : Fin de la bande d'intégration (Bw = f_max - f_min)

    N = length(x);
    yDSP = calc_spectre(x);

    yDSP = yDSP/N; % Parseval
    sig_bin = fix(Fsig/Fe * N) + 1;
    
    sig_bin_win = sig_bin + (-2:2);

    % Conversion des fréquences limites (Hz) en indices de bins
    idx_min = round(f_min/Fe*N) + 1;
    idx_max = round(f_max/Fe*N) + 1;
    
    % Sécurisation des indices (pour ne pas dépasser la taille du vecteur)
    idx_min = max(1, idx_min);
    idx_max = min(N, idx_max);
    
    % Calcul des bins de bruit : Bande totale moins les bins du signal
    err_bin = setdiff(idx_min:idx_max, sig_bin_win);
    

    % Calcul du SNR
    SNR = 10*log10(sum(yDSP(sig_bin_win))/sum(yDSP(err_bin)));

end