%% Init
clc; clear; 
baseDataPath = ("../MITBIH");
fileName = 'MITBIH_ALL.mat';
data = load(fullfile(baseDataPath, fileName)).data;
ecg=data.ECG_MIT_01;
sample_plot=2000;

filt = panTompkinsFilterBW(ecg);
fs=360;
ecg_h= filt.bandpass;
ecg_m = filt.mwi;
N  = length(ecg_m);

%% Parâmetros (PTA clássico vs modificado)

% Escolha rápida: defina PROFILE = 'classic' ou 'modified'
PROFILE = 'classic'; % 'modified' usa janela integradora menor e thresholds mais "apertados"
%PROFILE = 'modified';
%[a, missed_count, recovered_count, SIG_LEV, NOISE_LEV, THR_SIG, SIG_LEV1, NOISE_LEV1, THR_SIG1, r_idx]=full_conta_0(fs, N, PROFILE, ecg_h, ecg_m);
[a, missed_count, recovered_count, r_idx, trace]=full_conta_15(fs, N, PROFILE, ecg_h, ecg_m);
fprintf('Detecções: %d | Timeouts: %d | Recuperados: %d\n', a, missed_count, recovered_count);

fprintf("PROCESSO FINALIZADO!\n");

%% Tabela com #QRS reais por record (MIT-BIH) — ajuste conforme seu conjunto
dados = [ ...
    01 2273; 02 1865; 03 2187; 04 2084; 05 2229; 06 2572; 07 2027; 08 2137; 09 1763; 10 2532; ... 
    11 2124; 12 2539; 13 1795; 14 1879; 15 1953; 16 2412; 17 1535; 18 2278; 19 1987; 20 1863; ...
    21 2476; 22 1518; 23 1619; 24 2601; 25 1963; 26 2136; 27 2980; 28 2656; 29 2332; 30 2955; ...
    31 3005; 32 2650; 33 2748; 34 3251; 35 2262; 36 3363; 37 2208; 38 2154; 39 2048; 40 2427; ...
    41 2483; 42 2605; 43 2053; 44 2256; 45 1571; 46 1780; 47 3079; 48 2753];

valor = dados(dados(:,1) == 1, 2);
%disp(valor);

[ahw2, miss2, recov2, r_idx_hw2, trace_hw2] = full_conta_15(fs, N, PROFILE, ecg_h, ecg_m);

fprintf('Total Real = %d batimentos Identificados= %d \n (missed=%d, recovered=%d)', valor, ahw2, miss2, recov2);

%plot_out_3(fs, N, ecg_m_29, ecg_h_29, r_idx_hw2, trace_hw2);
plot_out_3(fs, N, ecg_m, ecg_h, r_idx_hw2, trace_hw2, 2, 430000, 460000);

%% Métricas simples por contagem (sem alinhamento batimento-a-batimento)
% ---------------------------------------------------------------
% Tabela com #QRS reais por record (MIT-BIH) — ajuste conforme seu conjunto
dados = [ ...
    01 2273; 02 1865; 03 2187; 04 2084; 05 2229; 06 2572; 07 2027; 08 2137; 09 1763; 10 2532; ... 
    11 2124; 12 2539; 13 1795; 14 1879; 15 1953; 16 2412; 17 1535; 18 2278; 19 1987; 20 1863; ...
    21 2476; 22 1518; 23 1619; 24 2601; 25 1963; 26 2136; 27 2980; 28 2656; 29 2332; 30 2955; ...
    31 3005; 32 2650; 33 2748; 34 3251; 35 2262; 36 3363; 37 2208; 38 2154; 39 2048; 40 2427; ...
    41 2483; 42 2605; 43 2053; 44 2256; 45 1571; 46 1780; 47 3079; 48 2753];

T = array2table(dados, 'VariableNames', {'Record','QRS_Real'});

% Número de QRS detectados
qrs_detect = numel(r_idx);

% Entrada do record (GUI ou prompt)
use_gui = usejava('desktop');
if use_gui
    answ = inputdlg({'Record (ex.: 100, 101, ...):'}, 'Entrada', [1 40], {'100'});
    if isempty(answ), rec = 100; else, rec = str2double(answ{1}); end
else
    rec = 100; % defina aqui se não houver GUI
end

idx = find(T.Record == rec, 1);
if isempty(idx)
    warning('Record %d não está na tabela. Métricas por contagem não serão exibidas.', rec);
    metrics = struct('Se',NaN,'PPV',NaN,'Acc_like',NaN,'DER',NaN,...
        'QRS_real',NaN,'QRS_detect',qrs_detect);
else
    qrs_real = T.QRS_Real(idx);
    
    % Aproximação por contagem (sem matching batimento-a-batimento)
    TP = min(qrs_detect, qrs_real);
    FP = max(qrs_detect - qrs_real, 0);
    FN = max(qrs_real   - qrs_detect, 0);
    
    Se       = 100 * TP / qrs_real;                 % Sensibilidade
    PPV      = 100 * TP / max(qrs_detect,1);        % Preditividade positiva
    Acc_like = 100 * TP / max(TP + FP + FN,1);      % "Acurácia" sem TN
    DER      = 100 * (FP + FN) / max(qrs_real,1);   % Detection Error Rate (%)
    
    fprintf('\n=== Record %d ===\n', rec);
    fprintf('QRS Real      : %d\n', qrs_real);
    fprintf('QRS Detectado : %d\n', qrs_detect);
    fprintf('TP (aprox)    : %d | FP: %d | FN: %d\n', TP, FP, FN);
    fprintf('Se            : %.2f %%\n', Se);
    fprintf('+P (PPV)      : %.2f %%\n', PPV);
    fprintf('Acc-like      : %.2f %%\n', Acc_like);
    fprintf('DER           : %.2f %%\n', DER);
    
    metrics = struct('Se',Se,'PPV',PPV,'Acc_like',Acc_like,'DER',DER,...
        'QRS_real',qrs_real,'QRS_detect',qrs_detect);
end

%% %% Tabela QRS reais (1–48)
i=1;
PROFILE = 'classic';
%PROFILE = 'modified';

beat = zeros(1,48);
missed_count = zeros(1,48);
recovered_count = zeros(1,48);
r_idx = cell(1,48);
trace = cell(1,48);

for i = 1:48
    ecg = data.(sprintf('ECG_MIT_%02d', i));

    filt = panTompkinsFilterBW(ecg);
    ecg_h = filt.bandpass(:);
    ecg_m = filt.mwi(:);

    N = length(ecg_m);

    [beat(i), missed_count(i), recovered_count(i), r_idx{i}, trace{i}] = ...
        full_conta_15(fs, N, PROFILE, ecg_h, ecg_m);

    fprintf('Record %02d processado | count=%d | missed=%d | recovered=%d\n', ...
        i, beat(i), missed_count(i), recovered_count(i));
end

dados = [ ...
    01 2273; 02 1865; 03 2187; 04 2084; 05 2229; 06 2572; 07 2027; 08 2137; 09 1763; 10 2532; ... 
    11 2124; 12 2539; 13 1795; 14 1879; 15 1953; 16 2412; 17 1535; 18 2278; 19 1987; 20 1863; ...
    21 2476; 22 1518; 23 1619; 24 2601; 25 1963; 26 2136; 27 2980; 28 2656; 29 2332; 30 2955; ...
    31 3005; 32 2650; 33 2748; 34 3251; 35 2262; 36 3363; 37 2208; 38 2154; 39 2048; 40 2427; ...
    41 2483; 42 2605; 43 2053; 44 2256; 45 1571; 46 1780; 47 3079; 48 2753];

T = array2table(dados, 'VariableNames', {'Record','QRS_Real'});

numRecords = height(T);

% ========================================
% 2) Pré-alocação dos vetores de métricas
% ========================================
QRS_detect = zeros(numRecords,1);
TP         = zeros(numRecords,1);
FP         = zeros(numRecords,1);
FN         = zeros(numRecords,1);
Se         = zeros(numRecords,1);
PPV        = zeros(numRecords,1);
Acc_like   = zeros(numRecords,1);
DER        = zeros(numRecords,1);

% ==========================================
% 3) Loop sobre as 48 amostras (beat 1x48)
% ==========================================
for i = 1:numRecords
    
    % ---- número de QRS detectados ----
    qrs_detect = beat(i);   % ajusta se for struct/cell
    
    QRS_detect(i) = qrs_detect;
    
    % ---- QRS real da tabela ----
    qrs_real = T.QRS_Real(i);
    
    % ---- Métricas aproximadas por contagem ----
    TP(i) = min(qrs_detect, qrs_real);
    FP(i) = max(qrs_detect - qrs_real, 0);
    FN(i) = max(qrs_real   - qrs_detect, 0);
    
    Se(i)       = 100 * TP(i) / qrs_real;                      % Sensibilidade
    PPV(i)      = 100 * TP(i) / max(qrs_detect,1);             % +P (PPV)
    Acc_like(i) = 100 * TP(i) / max(TP(i) + FP(i) + FN(i),1);  % "Acurácia" sem TN
    DER(i)      = 100 * (FP(i) + FN(i)) / max(qrs_real,1);     % Detection Error Rate (%)
end

% ==================================
% 4) Monta tabela final
% ==================================
MetricsTable = table( ...
    T.Record, T.QRS_Real, QRS_detect, TP, FP, FN, ...
    Se, PPV, Acc_like, DER, ...
    'VariableNames', {'Record','QRS_Real','QRS_Detect', ...
                      'TP','FP','FN','Se','PPV','Acc_like','DER'});

% ==================================
% 5) Adiciona linha com as MÉDIAS
% ==================================
meanSe   = mean(Se);
meanPPV  = mean(PPV);
meanAcc  = mean(Acc_like);
meanDER  = mean(DER);

% Linha de médias: outros campos como NaN ou vazio
MeanRow = {NaN, NaN, NaN, NaN, NaN, NaN, ...
           meanSe, meanPPV, meanAcc, meanDER};

MeanTable = cell2table(MeanRow, ...
    'VariableNames', MetricsTable.Properties.VariableNames);

% Junta a tabela original com a linha de médias
MetricsTable = [MetricsTable; MeanTable];

% ==================================
% 6) Exporta para CSV
% ==================================
writetable(MetricsTable, 'metricas_QRS_soft_and_fullconta.csv');

disp('Arquivo metricas_QRS_1a48.csv gerado com sucesso.');
open('metricas_QRS_soft_and_fullconta.csv');
fprintf("PROCESSO FINALIZADO!\n");

%% Plot
sample_plot = 2000;

figure;
subplot(5,1,1);
plot(filt.ecg(1:sample_plot));
title('ECG original');

subplot(5,1,2);
plot(filt.bandpass(1:sample_plot));
title('Bandpass');

subplot(5,1,3);
plot(filt.derivative(1:sample_plot));
title('Derivative');

subplot(5,1,4);
plot(filt.square(1:sample_plot));
title('Square');

subplot(5,1,5);
plot(filt.mwi(1:sample_plot));
title('MWI');

%% Pan Tompkins Original
fs=360;
gr=1;
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecg, fs, 1, 6000, 1);

%%
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecg, fs, 1, 6000, 20000);

