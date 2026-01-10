function Xpsd = calc_spectre(y_can)

win = blackman(length(y_can),'periodic');   % Création de la fenêtre
y_windowed = y_can(:).*win(:);              % Fenetrage du signal
Xpsd = abs(fft(y_windowed)).^2;