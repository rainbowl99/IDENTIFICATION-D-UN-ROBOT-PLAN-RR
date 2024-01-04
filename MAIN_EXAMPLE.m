%close all 
clear all
clc

% Define a trajectory
TR.T = 2; %s
TR.N = 3; % nb periods
TR.Q = [pi/2; 0];
TR.C = [0 1 0 0; 1 0 0 0];
TR.S = [0 0 0 0; 0 1 0 2]; 

% Run script
option = 'only_acquisition_hardware' ; % 'full_robot' or  'only_acquisition_hardware'
[q,tau] = myrobot(TR,option);

% Plot results
a = figure(1) ;
subplot(211); plot(q); xlabel('Time Sample'); ylabel('q'); axis tight 
subplot(212); plot(tau); xlabel('Time Sample'); ylabel('\tau'); axis tight

%% Traitement des signaux
% 3.1 Proposer un ou des filtres permettant d'enlever au maximum les bruits parasites
% Pour enlever le bruit parasite du reseau electrique 50 Hz, on peut utiliser un filtre passe-bas

% Perform FFT on the signals
Q_fft_1 = fft(q(:,1));
Q_fft_2 = fft(q(:,2));
Tau_fft_1 = fft(tau(:,1));
Tau_fft_2 = fft(tau(:,2));
% Frequency vector
Fs = 1000; % Sampling frequency
L_q = length(q(:,1)); % Length of signal
L_tau = length(tau(:,1));
f_q = Fs*(0:(L_q-1)/2)/L_q;
f_tau = Fs*(0:(L_tau-1)/2)/L_tau;
% Two-sided spectrum P2 and single-sided spectrum P1
P2_Q_1 = abs(Q_fft_1/L_q);
P1_Q_1 = P2_Q_1(1:(L_q-1)/2+1);
P1_Q_1(2:end-1) = 2*P1_Q_1(2:end-1);

P2_Q_2 = abs(Q_fft_2/L_q);
P1_Q_2 = P2_Q_2(1:(L_q-1)/2+1);
P1_Q_2(2:end-1) = 2*P1_Q_2(2:end-1);

P2_Tau_1 = abs(Tau_fft_1/L_tau);
P1_Tau_1 = P2_Tau_1(1:(L_tau-1)/2+1);
P1_Tau_1(2:end-1) = 2*P1_Tau_1(2:end-1);

P2_Tau_2 = abs(Tau_fft_2/L_tau);
P1_Tau_2 = P2_Tau_2(1:(L_tau-1)/2+1);
P1_Tau_2(2:end-1) = 2*P1_Tau_2(2:end-1);

% Plot FFT results
b = figure(2);
subplot(211); plot(f_q,P1_Q_1) 
title('Single-Sided Amplitude Spectrum of Q1(f)')
xlabel('f (Hz)')
ylabel('|P1\_Q(f)|')

subplot(212); plot(f_q,P1_Q_1) 
title('Single-Sided Amplitude Spectrum of Q1(f)')
xlabel('f (Hz)')
ylabel('|P1\_Q(f)|')
xlim([0,10]);

c = figure(3);
subplot(211); plot(f_q,P1_Q_2) 
title('Single-Sided Amplitude Spectrum of Q2(f)')
xlabel('f (Hz)')
ylabel('|P2\_Q(f)|')

subplot(212); plot(f_q,P1_Q_2) 
title('Single-Sided Amplitude Spectrum of Q2(f)')
xlabel('f (Hz)')
ylabel('|P2\_Q(f)|')
xlim([0,10]);

d = figure(4);
subplot(211); plot(f_tau,P1_Tau_1) 
title('Single-Sided Amplitude Spectrum of Tau1(f)')
xlabel('f (Hz)')
ylabel('|P1\_Tau(f)|')

subplot(212); plot(f_tau,P1_Tau_1) 
title('Single-Sided Amplitude Spectrum of Tau1(f)')
xlabel('f (Hz)')
ylabel('|P1\_Tau(f)|')
xlim([0,10]);

e = figure(5);
subplot(211); plot(f_tau,P1_Tau_2) 
title('Single-Sided Amplitude Spectrum of Tau2(f)')
xlabel('f (Hz)')
ylabel('|P2\_Tau(f)|')

subplot(212); plot(f_tau,P1_Tau_2) 
title('Single-Sided Amplitude Spectrum of Tau2(f)')
xlabel('f (Hz)')
ylabel('|P2\_Tau(f)|')
xlim([0,10]);

f = figure(6);
plot(f_q,P1_Q_1)
title('Single-Sided Amplitude Spectrum of bruit blanc Gaussien(f)')
xlabel('f (Hz)')
ylabel('|bruit blanc Gaussien(f)|')
ylim([-0.02,0.02]);
% les parametres
F_cutoff = 3; % Frequence de coupure en Hz
Fs_q = 1000 ; % frequence d'echantillonnage pour l'angle en Hz
Fs_tau = 2500 ; % frequence d'echantillonnage pour la couple en Hz
Fs_noise = 50;

% Conception de filtres à trappe pour atténuer le bruit de 50 Hz
notchFilt = designfilt('bandstopiir', 'FilterOrder', 2, ...
                       'HalfPowerFrequency1', Fs_noise-1, 'HalfPowerFrequency2', Fs_noise+1, ...
                       'SampleRate', Fs_q);

% Application des filtres de trappe
q_notchfiltered = filtfilt(notchFilt, q);
tau_notchfiltered = filtfilt(notchFilt, tau);

% definition un filtre passe-bas a l'orde 2 avec une frequence de coupure a 3 Hz
[B_q,A_q] = butter(2,F_cutoff/(Fs_q/2),'low'); 
[B_tau,A_tau] = butter(2,F_cutoff/(Fs_tau/2),'low');
% Appliquer ce filtre passe-bas aux signaux
angle_position_filtree = filtfilt(B_q,A_q,q_notchfiltered);
tau_filtree = filter(B_tau,A_tau,tau_notchfiltered);

g = figure(7) ;
subplot(211); plot(angle_position_filtree); xlabel('Time Sample'); ylabel('q'); axis tight 
subplot(212); plot(tau_filtree); xlabel('Time Sample'); ylabel('\tau'); axis tight

%% 3.2 Proposer une méthode pour que les signaux de couple et d angle soient à la même fréquence d$échantillonnage.
Fs_cible = 1000; % la frequence d'echantionnage cible
tau_resample = resample(tau_filtree,Fs_cible,Fs_tau);
angle_position_resample = resample(angle_position_filtree,Fs_cible,Fs_q);

% Plot results
h = figure(8) ;
subplot(211); plot(angle_position_resample); xlabel('Time Sample'); ylabel('q'); axis tight 
subplot(212); plot(tau_resample); xlabel('Time Sample'); ylabel('\tau'); axis tight

%% 3.3 Estimer le décalage temporel entre les deux interfaces d acquisition et proposer une méthode pour le compenser.
% Estimer le decalage temporel entre les signaux de angle position et de
% couple
[correction1,lags1] = xcorr(tau_resample(:,1),angle_position_resample(:,1)); % q1 et tau1
[correction2,lags2] = xcorr(tau_resample(:,2),angle_position_resample(:,2)); % q2 et tau2

[~, I1] = max(abs(correction1));
Delay_time1 = lags1(I1);% 

[~, I2] = max(abs(correction2));
Delay_time2 = lags2(I2);

% Compensation du délai
if Delay_time1 > 0
    % tau_resample est en retard par rapport à angle_position_resample
    tau_resample_adjusted_1 = [tau_resample(Delay_time1+1:end, 1); zeros(Delay_time1, 1)];
    angle_position_adjusted_1 = angle_position_resample(:, 1);
elseif Delay_time1 < 0
    % angle_position_resample est en retard par rapport à tau_resample
    angle_position_adjusted_1 = [angle_position_resample(-Delay_time1+1:end, 1); zeros(-Delay_time1, 1)];
    tau_resample_adjusted_1 = tau_resample(:, 1);
else
    % Pas de délai
    tau_resample_adjusted_1 = tau_resample(:, 1);
    angle_position_adjusted_1 = angle_position_resample(:, 1);
end

if Delay_time2 > 0
    % tau_resample est en retard par rapport à angle_position_resample
    tau_resample_adjusted_2 = [tau_resample(Delay_time2+1:end, 2); zeros(Delay_time2, 1)];
    angle_position_adjusted_2 = angle_position_resample(:, 2);
elseif Delay_time2 < 0
    % angle_position_resample est en retard par rapport à tau_resample
    angle_position_adjusted_2 = [angle_position_resample(-Delay_time2+1:end, 2); zeros(-Delay_time2, 1)];
    tau_resample_adjusted_2 = tau_resample(:, 2);
else
    % Pas de délai
    tau_resample_adjusted_2 = tau_resample(:, 2);
    angle_position_adjusted_2 = angle_position_resample(:, 2);
end

% Affichage des signaux ajustés
i = figure(9);
subplot(211);
plot(angle_position_adjusted_1);
hold on;
plot(tau_resample_adjusted_1);
title('Signaux Alignés dans le Temps');
legend('Position Angulaire', 'Tau');
xlim([0,6000]);

% Répétition pour le second signal
subplot(212);
plot(angle_position_adjusted_2);
hold on;
plot(tau_resample_adjusted_2);
title('Time Aligned Signals');
legend('Angle Position', 'Tau');
xlim([0,6000]);
%% 3.4 Supposons que angle_position_adjusted_1 soit le signal d'angle déjà aligné
% Calcul de l'intervalle de temps
dt = 1 / Fs_cible;  % Intervalle de temps d'échantillonnage

% Calcul de la vitesse angulaire - méthode de différence centrale
angular_velocity_1 = diff(angle_position_adjusted_1) / dt;
% Comme la différence réduit le nombre de points de données, nous ajoutons des zéros au début et à la fin
angular_velocity_1 = [0; angular_velocity_1; 0];

% Calcul de l'accélération angulaire - méthode de différence centrale
angular_acceleration_1 = diff(angular_velocity_1) / dt;
% Ajout de zéros de la même manière
angular_acceleration_1 = [0; angular_acceleration_1; 0];

% Calcul de la vitesse angulaire - méthode de différence centrale
angular_velocity_2 = diff(angle_position_adjusted_2) / dt;
% Comme la différence réduit le nombre de points de données, nous ajoutons des zéros au début et à la fin
angular_velocity_2 = [0; angular_velocity_2; 0];

% Calcul de l'accélération angulaire - méthode de différence centrale
angular_acceleration_2 = diff(angular_velocity_2) / dt;
% Ajout de zéros de la même manière
angular_acceleration_2 = [0; angular_acceleration_2; 0];


% Affichage des vitesses et accélérations angulaires
samples_range = 250:5750;% En raison des énormes fluctuations au point limite, nous choisissons les données pour l'intervalle [250,5750].

j = figure(10);
subplot(211);
plot(samples_range, angular_velocity_1(samples_range));
title('Vitesse Angulaire Q1 (différence)');
xlabel('Échantillon de Temps');
ylabel('Vitesse (rad/s)');

subplot(212);
plot(samples_range, angular_acceleration_1(samples_range));
title('Accélération Angulaire Q1 (différence)');
xlabel('Échantillon de Temps');
ylabel('Accélération (rad/s^2)');

k = figure(11);
subplot(211);
plot(samples_range, angular_velocity_2(samples_range));
title('Vitesse Angulaire Q2 (différence)');
xlabel('Échantillon de Temps');
ylabel('Vitesse (rad/s)');

subplot(212);
plot(samples_range, angular_acceleration_2(samples_range));
title('Accélération Angulaire Q2 (différence)');
xlabel('Échantillon de Temps');
ylabel('Accélération (rad/s^2)');

% Supposons que angle_position_adjusted_1 soit le signal d'angle déjà aligné
% Calcul de l'intervalle de temps
dt = 1 / Fs_cible;  % Intervalle de temps d'échantillonnage

% Initialisation des tableaux
n = length(angle_position_adjusted_1);
angular_velocity_1 = zeros(n, 1);
angular_acceleration_1 = zeros(n, 1);

% Calcul des vitesses et accélérations angulaires pour Q1
for i = 2:n-1
    % Calcul de la vitesse angulaire en utilisant le développement en série de Taylor
    angular_velocity_1(i) = (angle_position_adjusted_1(i+1) - angle_position_adjusted_1(i-1)) / (2*dt);

    % Calcul de l'accélération angulaire en utilisant le développement en série de Taylor
    angular_acceleration_1(i) = (angle_position_adjusted_1(i+1) - 2*angle_position_adjusted_1(i) + angle_position_adjusted_1(i-1)) / (dt^2);
end

% Affichage des vitesses et accélérations angulaires
l = figure(12);
subplot(211);
plot(angular_velocity_1);
title('Vitesse Angulaire Q1 (Expansion de Taylor)');
xlabel('Échantillon de Temps');
ylabel('Vitesse (rad/s)');

subplot(212);
plot(angular_acceleration_1);
title('Accélération Angulaire Q1 (Expansion de Taylor)');
xlabel('Échantillon de Temps');
ylabel('Accélération (rad/s^2)');

% Initialisation des tableaux
n = length(angle_position_adjusted_2);
angular_velocity_2 = zeros(n, 1);
angular_acceleration_2 = zeros(n, 1);


% Calcul des vitesses et accélérations angulaires pour Q2
for i = 2:n-1
    % Calcul de la vitesse angulaire en utilisant le développement en série de Taylor
    angular_velocity_2(i) = (angle_position_adjusted_2(i+1) - angle_position_adjusted_2(i-1)) / (2*dt);

    % Calcul de l'accélération angulaire en utilisant le développement en série de Taylor
    angular_acceleration_2(i) = (angle_position_adjusted_2(i+1) - 2*angle_position_adjusted_2(i) + angle_position_adjusted_2(i-1)) / (dt^2);
end

% Affichage des vitesses et accélérations angulaires
m = figure(13);
subplot(211);
plot(angular_velocity_2);
title('Angular Velocity Q2(Taylor Expansion)');
xlabel('Time Sample');
ylabel('Velocity (rad/s)');

subplot(212);
plot(angular_acceleration_2);
title('Angular Acceleration Q2(Taylor Expansion)');
xlabel('Time Sample');
ylabel('Acceleration (rad/s^2)');