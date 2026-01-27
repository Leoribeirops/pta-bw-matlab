function [Beat_C1] = Thresholding_BP_MWI(ecg_h, ecg_m, fs, gr, plot_win)
%% Thresholding_BP_MWI
% Versão do Thresholding do Pan-Tompkins que recebe:
%   ecg_h : sinal já filtrado em bandpass (mesmo comprimento do ECG original)
%   ecg_m : sinal já processado até o MWI (moving window integration)
%   fs    : frequência de amostragem (Hz)
%   gr    : flag de plot (1 = plota, 0 = não plota)
%   plot_win : janela de plotagem
%             - 0 ou []  : plota o sinal inteiro (se gr=1)
%             - [ini fim]: plota somente o intervalo em amostras (se gr=1)
%
% Saída:
%   Beat_C1 : número de batimentos (QRS) confirmados no bandpass

%% ---------------------- Defaults e validações ---------------------- %%
if nargin < 4 || isempty(gr)
    gr = 1;
end
if nargin < 5
    plot_win = 0;
end

if ~isvector(ecg_h) || ~isvector(ecg_m)
    error('ecg_h e ecg_m devem ser vetores (linha ou coluna).');
end

ecg_h = ecg_h(:);
ecg_m = ecg_m(:);

if length(ecg_h) ~= length(ecg_m)
    error('ecg_h e ecg_m devem ter o mesmo comprimento.');
end

if ~isscalar(fs) || fs <= 0
    error('fs deve ser um escalar positivo.');
end

% Regra solicitada: se gr==0 e plot_win==0 => não gera gráfico.
% (Na prática: se gr==0 não plota nunca.)
doPlot = (gr ~= 0);

N = length(ecg_h);

% Normalização leve (opcional, mas mantém comportamento próximo ao original)
if max(abs(ecg_h)) ~= 0
    ecg_h = ecg_h / max(abs(ecg_h));
end
if max(abs(ecg_m)) ~= 0
    ecg_m = ecg_m / max(abs(ecg_m));
end

%% ---------------------- Janela de plotagem ------------------------- %%
% plot_win:
%   0 ou []      => plota tudo (se doPlot)
%   [ini fim]    => plota somente [ini:fim] em amostras
if isempty(plot_win)
    plot_win = 0;
end

useWin = false;
if isnumeric(plot_win) && isscalar(plot_win) && plot_win == 0
    useWin = false;
elseif isnumeric(plot_win) && numel(plot_win) == 2
    plot_win = round(plot_win(:)');
    plot_win(1) = max(1, plot_win(1));
    plot_win(2) = min(N, plot_win(2));
    if plot_win(1) >= plot_win(2)
        error('plot_win inválida: inicio >= fim.');
    end
    useWin = true;
else
    error('plot_win deve ser 0, [] ou [inicio fim] em amostras.');
end

%% -------------------- Parâmetros/Inicialização --------------------- %%
skip = 0;                 % vira 1 quando detecta onda T
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0;

% "Delay" não é saída aqui, mas mantemos a mesma convenção (janela ~150ms)
Nwin = round(0.150 * fs);

%% --------------------------- Plots iniciais ------------------------ %%
if doPlot
    if useWin
        idx = plot_win(1):plot_win(2);
    else
        idx = 1:N;
    end

    figure;
    ax1 = subplot(211);
    plot(idx, ecg_h(idx));
    axis tight;
    title('Entrada: Bandpass (ecg\_h)');

    ax2 = subplot(212);
    plot(idx, ecg_m(idx));
    axis tight;
    title('Entrada: MWI (ecg\_m)');

    linkaxes([ax1 ax2],'x');
    zoom on;
end

%% ===================== Fiducial Marks (MWI) ======================== %%
[pks, locs] = findpeaks(ecg_m, 'MINPEAKDISTANCE', round(0.2*fs));
LLp = length(pks);

if LLp == 0
    Beat_C1 = 0;
    return;
end

%% ================= Buffers e variáveis do algoritmo ================= %%
% Detecção no MWI (interno)
qrs_c = zeros(1,LLp);
qrs_i = zeros(1,LLp);

% Ruído (interno)
nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);

% Tracking de níveis e thresholds (para plot/depuração)
SIGL_buf   = zeros(1,LLp);
NOISL_buf  = zeros(1,LLp);
THRS_buf   = zeros(1,LLp);
SIGL_buf1  = zeros(1,LLp);
NOISL_buf1 = zeros(1,LLp);
THRS_buf1  = zeros(1,LLp);

%% ================= Training phase (2 s) ============================ %%
L2 = min(N, 2*fs);

THR_SIG   = max(ecg_m(1:L2)) * (1/3);
THR_NOISE = mean(ecg_m(1:L2)) * (1/2);
SIG_LEV   = THR_SIG;
NOISE_LEV = THR_NOISE;

THR_SIG1   = max(ecg_h(1:L2)) * (1/3);
THR_NOISE1 = mean(ecg_h(1:L2)) * (1/2);
SIG_LEV1   = THR_SIG1;
NOISE_LEV1 = THR_NOISE1;

%% =================== Thresholding and decision ====================== %%
Beat_C  = 0;   % QRS no MWI (interno)
Beat_C1 = 0;   % QRS confirmados no bandpass (SAÍDA)
Noise_Count = 0;

for i = 1:LLp

    %% localizar pico correspondente no bandpass em torno do pico do MWI
    win150 = round(0.150*fs);
    if (locs(i) - win150) >= 1 && locs(i) <= N
        seg = ecg_h(locs(i)-win150:locs(i));
        [y_i, x_i] = max(seg);
    else
        if i == 1
            seg = ecg_h(1:locs(i));
            [y_i, x_i] = max(seg);
            ser_back = 1;
        elseif locs(i) >= N
            seg = ecg_h(max(1,locs(i)-win150):end);
            [y_i, x_i] = max(seg);
        else
            seg = ecg_h(1:locs(i));
            [y_i, x_i] = max(seg);
        end
    end

    %% update RR (igual ao original)
    if Beat_C >= 9
        diffRR  = diff(qrs_i(Beat_C-8:Beat_C));
        mean_RR = mean(diffRR);
        comp    = qrs_i(Beat_C) - qrs_i(Beat_C-1);

        if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
            THR_SIG  = 0.5 * THR_SIG;
            THR_SIG1 = 0.5 * THR_SIG1;
        else
            m_selected_RR = mean_RR;
        end
    end

    if m_selected_RR
        test_m = m_selected_RR;
    elseif mean_RR && m_selected_RR == 0
        test_m = mean_RR;
    else
        test_m = 0;
    end

    %% Search-back (igual ao original)
    if test_m
        if Beat_C > 0 && (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)
            lo = qrs_i(Beat_C) + round(0.200*fs);
            hi = locs(i) - round(0.200*fs);

            if lo >= 1 && hi > lo && hi <= N
                [pks_temp, rel] = max(ecg_m(lo:hi));
                locs_temp = lo + rel - 1;

                if pks_temp > THR_NOISE
                    Beat_C = Beat_C + 1;
                    qrs_c(Beat_C) = pks_temp;
                    qrs_i(Beat_C) = locs_temp;

                    % localizar também no bandpass
                    seg_lo = max(1, locs_temp - win150);
                    seg_hi = min(N, locs_temp);
                    [y_i_t, x_i_t] = max(ecg_h(seg_lo:seg_hi));

                    if y_i_t > THR_NOISE1
                        Beat_C1 = Beat_C1 + 1;
                        % índice do R no bandpass
                        qrs_idx = seg_lo + (x_i_t - 1);
                        % atualiza nível bandpass
                        SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;
                    end

                    SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV;
                end
            end
        end
    end

    %% Classificação pelo MWI
    if pks(i) >= THR_SIG

        % checagem de onda T (igual ao original)
        if Beat_C >= 3
            if (locs(i) - qrs_i(Beat_C)) <= round(0.3600*fs)
                s1_lo = max(1, locs(i)-round(0.075*fs));
                s1_hi = locs(i);
                s2_lo = max(1, qrs_i(Beat_C)-round(0.075*fs));
                s2_hi = qrs_i(Beat_C);

                Slope1 = mean(diff(ecg_m(s1_lo:s1_hi)));
                Slope2 = mean(diff(ecg_m(s2_lo:s2_hi)));

                if abs(Slope1) <= abs(0.5*Slope2)
                    Noise_Count = Noise_Count + 1;
                    nois_c(Noise_Count) = pks(i);
                    nois_i(Noise_Count) = locs(i);
                    skip = 1;

                    NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                    NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;
                else
                    skip = 0;
                end
            end
        end

        if skip == 0
            Beat_C = Beat_C + 1;
            qrs_c(Beat_C) = pks(i);
            qrs_i(Beat_C) = locs(i);

            % Confirmação no bandpass
            if y_i >= THR_SIG1
                Beat_C1 = Beat_C1 + 1;

                % atualiza nível no bandpass
                SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;
            end

            % atualiza nível no MWI
            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV;
        end

    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;

    else
        Noise_Count = Noise_Count + 1;
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);

        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;
    end

    %% Ajuste de thresholds (SNR)
    if (NOISE_LEV ~= 0) || (SIG_LEV ~= 0)
        THR_SIG   = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*THR_SIG;
    end
    if (NOISE_LEV1 ~= 0) || (SIG_LEV1 ~= 0)
        THR_SIG1   = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*THR_SIG1;
    end

    % buffers p/ visualização
    SIGL_buf(i)   = SIG_LEV;
    NOISL_buf(i)  = NOISE_LEV;
    THRS_buf(i)   = THR_SIG;

    SIGL_buf1(i)  = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i)  = THR_SIG1;

    % reset
    skip = 0;
    ser_back = 0;

end

%% ----------------------------- Plot final --------------------------- %%
if doPlot
    if useWin
        idx = plot_win(1):plot_win(2);
    else
        idx = 1:N;
    end

    % Para marcar detecções no plot, reconstruímos um conjunto aproximado
    % de instantes QRS no MWI (qrs_i) e mostramos o threshold track
    qrs_i_plot = qrs_i(1:Beat_C);
    qrs_c_plot = qrs_c(1:Beat_C);

    % restringe à janela (se houver)
    q_in = (qrs_i_plot >= idx(1)) & (qrs_i_plot <= idx(end));
    qrs_i_plot = qrs_i_plot(q_in);
    qrs_c_plot = qrs_c_plot(q_in);

    locs_in = (locs >= idx(1)) & (locs <= idx(end));
    locs_w = locs(locs_in);

    % buffers têm tamanho LLp; alinhar por "locs"
    SIGL_w   = SIGL_buf(locs_in);
    NOISL_w  = NOISL_buf(locs_in);
    THRS_w   = THRS_buf(locs_in);

    SIGL1_w  = SIGL_buf1(locs_in);
    NOISL1_w = NOISL_buf1(locs_in);
    THRS1_w  = THRS_buf1(locs_in);

    figure;

    subplot(211);
    plot(idx, ecg_h(idx)); axis tight;
    title('Bandpass (ecg\_h) + níveis/limiares');
    hold on;
    if ~isempty(locs_w)
        plot(locs_w, NOISL1_w, '--k', 'LineWidth', 2);
        plot(locs_w, SIGL1_w, '-.r', 'LineWidth', 2);
        plot(locs_w, THRS1_w, '-.g', 'LineWidth', 2);
    end

    subplot(212);
    plot(idx, ecg_m(idx)); axis tight;
    title('MWI (ecg\_m) + QRS (interno) + níveis/limiares');
    hold on;
    if ~isempty(qrs_i_plot)
        scatter(qrs_i_plot, qrs_c_plot, 'm');
    end
    if ~isempty(locs_w)
        plot(locs_w, NOISL_w, '--k', 'LineWidth', 2);
        plot(locs_w, SIGL_w, '-.r', 'LineWidth', 2);
        plot(locs_w, THRS_w, '-.g', 'LineWidth', 2);
    end
end

end
