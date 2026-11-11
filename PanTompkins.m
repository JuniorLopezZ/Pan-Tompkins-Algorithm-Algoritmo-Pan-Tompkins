clear;
axis equal;

%% ==================== LECTURA DE DATOS ======================
opts = detectImportOptions('Lola.txt', ...
    'Delimiter','\t', ...
    'DecimalSeparator',',');                      

data = readtable('Lola.txt', opts);
clf;

foto = data.Var1;
ecg  = data.Var2;

fs = 1000;
N  = length(ecg);
t_ecg = (0:N-1)/fs;

%% ==================== ECG LIMPIA ============================
t_ecg_clean = t_ecg(290001:320001)';
ecg_clean   = ecg(290001:320001)';

figure(1);
subplot(2,1,1);
hold on;
plot(t_ecg_clean, ecg_clean);
title('Señal ECG 290-320s');
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');

% Detección de picos
umbral_y = 3*mean(abs(ecg_clean));
umbral_x = 0.4*fs;
[pks_ecg, locs_ecg] = findpeaks(ecg_clean, 'MinPeakHeight', umbral_y, 'MinPeakDistance', umbral_x);
scatter(t_ecg_clean(locs_ecg), pks_ecg);

% Intervalos R-R y BPM
HVR = diff(t_ecg_clean(locs_ecg));
HVR_BPM = 60./HVR;

% Interpolación
new_t = downsample(t_ecg_clean,1);
HRV_f = interp1(t_ecg_clean(locs_ecg(1:end-1)), HVR, new_t);
BPM = (1./HRV_f)*60;

subplot(2,1,2);
plot(new_t,BPM);
title('Frecuencia cardiaca instantánea (ECG limpia)');
xlabel('Tiempo (s)');
ylabel('BPM');
hold off;

locs_ecg=locs_ecg(1:end-1)';
figure(10);
plot(locs_ecg,HVR);
BPM_TEST=(1./HVR)*60;
plot(locs_ecg,BPM_TEST);

% Número de latidos y BPM medio
t_total = t_ecg_clean(end) - t_ecg_clean(1);
num_latidos = length(pks_ecg);
BPM_MEDIO = (num_latidos / t_total)*60;
fprintf('El valor de la tasa cardiaca media en la señal ECG sin ruido es de %.2f\n', BPM_MEDIO);

%% ==================== ECG CON RUIDO =========================
figure(2);
subplot(2,1,1);
hold on;
plot(t_ecg, ecg);
title('Señal ECG traza completa');
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');

umbral_y_2 = mean(abs(ecg));
umbral_x_2 = 0.3*fs;
[pks_ecg_d, locs_ecg] = findpeaks(ecg, 'MinPeakHeight', umbral_y_2, 'MinPeakDistance', umbral_x_2);
scatter(t_ecg(locs_ecg), pks_ecg_d);
t_t=t_ecg(locs_ecg)';


HVR = diff(t_ecg(locs_ecg));
new_t_2 = downsample(t_ecg,1);
HRV_f_2 = interp1(t_ecg(locs_ecg(1:end-1)), HVR, new_t_2);
BPM_2 = (1./HRV_f_2)*60;

subplot(2,1,2);
bpm_2=(1./HVR)*60;
t_bpm=t_t(2:end);
plot(t_t,pks_ecg_d,'g-',t_bpm,bpm_2,'b-',new_t_2,BPM_2,'r-')
title('Frecuencia cardiaca instantánea (ECG con ruido)');
xlabel('Tiempo (s)');
ylabel('BPM');

t_total_2 = t_ecg(end) - t_ecg(1);
num_latidos_2 = length(pks_ecg_d);
BPM_MEDIO_2 = (num_latidos_2 / t_total_2)*60;
fprintf('El valor de la tasa cardiaca media en la señal ECG con ruido es de %.2f\n', BPM_MEDIO_2);
hold off;

%% ==================== FOTOPLETISMOGRAFÍA CON RUIDO =====================
n = length(foto);
t_foto = (0:n-1)/fs;

% Señal completa con ruido
figure(3);
subplot(2,1,1);
hold on;
plot(t_foto, foto);
title('Señal FOTO completa con ruido');
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');

umbral_y_foto = mean(abs(foto));
umbral_x_foto = 0.1*fs;
[pks_ecg, locs_ecg] = findpeaks(foto, 'MinPeakHeight', umbral_y_foto, 'MinPeakDistance', umbral_x_foto);


scatter(t_foto(locs_ecg), pks_ecg);

diff_entre_picos = diff(t_foto(locs_ecg));
new_t_foto = downsample(t_foto,1);
HVR_f_foto = interp1(t_foto(locs_ecg(1:end-1)), diff_entre_picos, new_t_foto);
BPM_foto = (1./HVR_f_foto)*60;

subplot(2,1,2);
plot(new_t_foto, BPM_foto);
title('Frecuencia cardiaca instantánea (FOTO con ruido)');
xlabel('Tiempo (s)');
ylabel('BPM');

t_total = t_foto(end) - t_foto(1);
n_latidos_foto = length(pks_ecg);
BPM_MEDIO_FOTO = (n_latidos_foto / t_total)*60;
fprintf('El valor de la tasa cardiaca media en la señal FOTO con ruido es de %.2f\n', BPM_MEDIO_FOTO);

%% ==================== FOTOPLETISMOGRAFÍA LIMPIA =====================
t_foto_clean = t_foto(190001:210001)';
foto_clean = foto(190001:210001)';

figure(4);
subplot(2,1,1);
hold on;
plot(t_foto_clean, foto_clean);
title('Señal FOTO 150-180s');
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');

umbral_y_foto_clean = mean(abs(foto_clean));
umbral_x_foto_clean = 0.3*fs;
[pks_ecg, locs_ecg] = findpeaks(foto_clean, 'MinPeakHeight', umbral_y_foto_clean, 'MinPeakDistance', umbral_x_foto_clean);
scatter(t_foto_clean(locs_ecg), pks_ecg);

diff_entre_picos = diff(t_foto_clean(locs_ecg));
new_t_foto_clean = downsample(t_foto_clean,1);
HVR_f_foto_clean = interp1(t_foto_clean(locs_ecg(1:end-1)), diff_entre_picos, new_t_foto_clean);
BPM_foto_clean = (1./HVR_f_foto_clean)*60;

subplot(2,1,2);
plot(new_t_foto_clean, BPM_foto_clean);
title('Frecuencia cardiaca instantánea (FOTO limpia)');
xlabel('Tiempo (s)');
ylabel('BPM');

t_total_clean = t_foto_clean(end) - t_foto_clean(1);
n_latidos_foto = length(pks_ecg);
BPM_MEDIO_FOTO_clean = (n_latidos_foto / t_total_clean)*60;
fprintf('El valor de la tasa cardiaca media en la señal FOTO sin ruido es de %.2f\n', BPM_MEDIO_FOTO_clean);

%% ==================== FUSIÓN ECG + FOTO =====================
figure(5);
hold on;
plot(t_foto, foto, 'green', 'LineWidth', 1.5);
plot(t_ecg, ecg, 'red', 'LineWidth', 1.5);
title('Señales combinadas ECG y FOTO');
xlabel('Tiempo (s)');
ylabel('Amplitud (mV)');

alpha = 0.7;
HR_fusion = alpha * BPM_2 + (1-alpha) * BPM_foto;

w_ecg = 0.5;
w_ppg = 0.8;
media_fusion = (w_ppg * BPM_MEDIO_FOTO + w_ecg * BPM_MEDIO_2) / (w_ecg + w_ppg);

fprintf('El valor de la tasa cardiaca media en la señal conjunta con ruido es de %.2f\n', media_fusion);

figure(7);
subplot(3,1,1)
plot(t_ecg,foto);
plot(t_ecg,ecg);
plot(t_foto, HR_fusion, 'k', 'LineWidth', 1.5);
title('Frecuencia cardiaca fusionada (ECG + FOTO)');
xlabel('Tiempo (s)');
ylabel('BPM ');
subplot(3,1,2)
plot(t_ecg,BPM_foto);
title('Frecuencia cardiaca foto');
xlabel('Tiempo (s)');
ylabel('BPM');
subplot(3,1,3)
plot(t_ecg,BPM_2);
title('Frecuencia cardiaca ecg');
xlabel('Tiempo (s)');
ylabel('BPM');



