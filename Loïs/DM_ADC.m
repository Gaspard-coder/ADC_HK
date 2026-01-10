close all;
clc;

%%% LELARGE Loïs & VIBERT Gaspard

% Paramètres du CAN 
Nbits   =   10;        % Nombre de bits pour le CAN
PE      =   2;         % Pleine Echelle

% --- Paramètres ---

Fe      =   30.72e6;            % Fréquence d'échantillonnage (Hz)
Fcarr   =   1e6;                % Fréquence de porteuse à modifier pour être sur un bin
N       =   2^16;               % Nombre de sample
t       =   0:1/Fe:(N-1)/Fe;    % temps de simulation avec discretisation     
A       =   PE/2;               % Amplitude du sinus

Fsig    =   round(Fcarr/Fe*length(t))/length(t)*Fe; % fréquence de porteuse adaptée pour être sur le bin

x       =   A*sin(2*pi*Fsig*t);   % Création du signal d'entrée (Amplitude 1, donc de -1 à 1)


%% Question 1 : Proposer une fonction qui modélise un CAN
y_can = CAN(x,Nbits,PE);
figure('Name', 'Question 2.1 - Modélisation du CAN');
plot(t(1:500), y_can(1:500), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Signal quantifié');
grid on;
xlabel('Temps (s)');
ylabel('Amplitude (V)');
title(sprintf('Signal après quantification - CAN %d bits (PE = %.1f V)', Nbits, PE));
legend('Location', 'best');
xlim([0 t(500)]);

% On peut zoomer pour mieux voir la quantification :
T_sig = 1/Fsig;
nb_points_par_periode = round(Fe * T_sig);
nb_pts_zoom = round(nb_points_par_periode * 1.5); % 1.5 période pour bien voir
Q = PE / 2^Nbits; % Quantum de quantification

figure('Name', 'Question 2.1 - zoom avec escalier', 'Position', [150 150 1200 800]);
subplot(2,1,1);
stairs(t(1:nb_pts_zoom)*1e6, y_can(1:nb_pts_zoom), 'r-', 'LineWidth', 2, 'DisplayName', 'Signal quantifié (escalier)');
hold on;
plot(t(1:nb_pts_zoom)*1e6, x(1:nb_pts_zoom), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Signal original');
plot(t(1:nb_pts_zoom)*1e6, y_can(1:nb_pts_zoom), 'ko', 'MarkerSize', 4, 'HandleVisibility', 'off');
grid on;
xlabel('Temps (µs)');
ylabel('Amplitude (V)');
title(sprintf('ZOOM - zoom avec escalier CAN %d bits (Q = %.4f V)', Nbits, Q));
legend('Location', 'best');

%% Question 2 : Traçage du spectre de sortie d'un CAN 10 bits
bin_freq_val = 1:1:N;

figure('Name', 'Question 2.2 - Calcul du spectre de la fonction');
Xpsd = calc_spectre(y_can);
plot(bin_freq_val,10*log10(Xpsd))
xlabel('Bin');
ylabel('Densité Spectrale de Puissance (dB)');
title('DSP');

%% Question 3 : Calcul du SNR
SNR = calc_SNR_freq(Xpsd, Fe, 0, 10e6, Fsig);

%% Question 4 : Calcul du SQNR et comparaison 
Nbits_vecteur = 1:10;     % Vecteur de bits de 1 à 10
SNR_mesure = zeros(size(Nbits_vecteur)); % Pour stocker les résultats

for i = 1:length(Nbits_vecteur)
    n_bits = Nbits_vecteur(i);    
    y_can = CAN(x,n_bits,PE);
    % --- Calcul du SNR ---
    SNR_mesure(i) = calc_SNR_freq(y_can,Fe,0,10e6,Fsig);
end
% --- 3. Comparaison Théorique ---
SNR_theorie = 6.02 * Nbits_vecteur + 1.76 + 20*log10(2*A/PE) + 10*log10(Fe/(2*10e6));

% --- 4. Affichage des résultats ---
figure('Color', 'w');
plot(Nbits_vecteur, SNR_mesure, '-o', 'LineWidth', 2, 'DisplayName', 'Mesure (calc\_SNR\_freq)');
hold on;
plot(Nbits_vecteur, SNR_theorie, '--r', 'LineWidth', 2, 'DisplayName', 'Théorie SQNR');
grid on;
xlabel('Nombre de bits (N)');
ylabel('SNR (dB)');
title('Évolution du SNR en fonction de la résolution du CAN');
legend('Location', 'best');


%% Question 5
A_vec = 0.01:0.01:PE/2*1.5; 

% Allocation mémoire
SNR_simu    = zeros(size(A_vec));

for i = 1:length(A_vec)
    A = A_vec(i);
    
    % A. Génération du signal
    x_in = A * sin(2*pi*Fsig*t);
    
    y_can = CAN(x_in, Nbits, PE);

    SNR_simu(i) = calc_SNR_freq(y_can, Fe, 0, 10e6, Fsig);
end
% --- Comparaison Théorique ---
SNR_theorie = 6.02 * Nbits + 1.76 + 20*log10(2*A_vec/PE) + 10*log10(Fe/(2*10e6));

% --- Affichage des résultats ---
figure('Color', 'w', 'Name', 'Evolution SNR vs Amplitude');
plot(A_vec, SNR_theorie, '--r', 'LineWidth', 2, 'DisplayName', 'Théorie (Eq. 1)');
hold on;
plot(A_vec, SNR_simu, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Simulation');
grid on;
xlabel('Amplitude d''entrée A_{in} (V)');
ylabel('SNR (dB)');
title(['SNR en sortie du CAN ' num2str(Nbits) ' bits (BW = 10 MHz)']);
legend('Location', 'SouthEast');
