%%
clc; clear; 
baseDataPath = 'data';
fileName = 'ECG_MIT_01.mat';
data = load(fullfile(baseDataPath, fileName));
ecg=data.ECG_MIT_1;
sample_plot=2000;

%% Pan Tompkins BW

% Coeficientes do filtro (equação dada)
b = [1/16, 0, -1/8, 0, 1/16];
a = [1, -3.176, 3.808, -2.076, 0.444];

% x: seu sinal de entrada (coluna ou linha)
bw_out = filter(b, a, ecg);

b = [1/4, 1/8, 0, -1/8, -1/4];
a = 1;

diff = filter(b, a, bw_out);
square = diff.^2;

N = 30;                 % tamanho da janela (ex.: 30 amostras)
b = (1/N) * ones(1, N+1);
a = 1;

MWI = filter(b, a, square);    % x = sinal ao quadrado

%
fs=360;
gr=1;
plot_win = [0 2000];   % amostras
beat = Thresholding_BP_MWI(bw_out, MWI, fs, gr, plot_win);

%% Plot
sample_plot = [1 2000];   % EM AMOSTRAS
% ou
%sample_plot = 0;           % sinal inteiro
plot_signal_window(ecg, sample_plot, 1, 'ECG MIT BIH 01');
plot_signal_window(bw_out,          sample_plot, 2, 'Bandpass BW');
plot_signal_window(diff,            sample_plot, 3, 'Differentiator');
plot_signal_window(square,          sample_plot, 4, 'Square');
plot_signal_window(MWI,             sample_plot, 5, 'Moving Window Integration');

%% Pan Tompkins Original
fs=360;
gr=1;
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecg, fs, 1, 6000, 1);

%%
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecg, fs, 1, 6000, 20000);

