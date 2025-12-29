function [qrs_amp_raw, qrs_i_raw, delay] = Thresholding_BP_MWI(ecg_h, ecg_m, fs, gr)
%% Thresholding_BP_MWI
% Entrada:
%   ecg_h : sinal já filtrado em bandpass (mesmo comprimento do ECG original)
%   ecg_m : sinal já processado até o MWI (moving window integration) usado para findpeaks/limiar
%   fs    : frequência de amostragem
%   gr    : (opcional) flag de plot (1 plota, 0 não plota)
%
% Saída:
%   qrs_amp_raw : amplitudes R no sinal bandpass (ecg_h)
%   qrs_i_raw   : índices dos picos R no sinal bandpass (ecg_h)
%   delay       : atraso (amostras) associado ao MWI (se você quiser usar para alinhar)

if nargin < 4
    gr = 1;
end

if ~isvector(ecg_h) || ~isvector(ecg_m)
    error('ecg_h e ecg_m devem ser vetores (linha ou coluna).');
end

ecg_h = ecg_h(:);
ecg_m = ecg_m(:);

if length(ecg_h) ~= length(ecg_m)
    error('ecg_h e ecg_m devem ter o mesmo comprimento.');
end

%% ======================= Initialize =============================== %
delay = 0;
skip = 0;                 % vira 1 quando detecta onda T
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0;
ax = zeros(1,2);

%% (Opcional) normalização leve (mantém comportamento similar ao original)
if max(abs(ecg_h)) ~= 0
    ecg_h = ecg_h / max(abs(ecg_h));
end
if max(abs(ecg_m)) ~= 0
    ecg_m = ecg_m / max(abs(ecg_m));
end

%% Se seu ecg_m foi feito com janela ~150 ms, o atraso de integração é ~N/2
Nwin = round(0.150 * fs);
delay = delay + Nwin/2;   % mantém convenção do código original

%% ===================== (Opcional) Plot de entrada ================== %
if gr
    figure;
    ax(1) = subplot(211); plot(ecg_h); axis tight; title('Entrada: Bandpass (ecg\_h)');
    ax(2) = subplot(212); plot(ecg_m); axis tight; title('Entrada: MWI (ecg\_m)');
    linkaxes(ax,'x'); zoom on;
end

%% ===================== Fiducial Marks ============================== %
[pks, locs] = findpeaks(ecg_m, 'MINPEAKDISTANCE', round(0.2*fs));

%% =================== Initialize Some Other Parameters =============== %
LLp = length(pks);

qrs_c       = zeros(1,LLp);   % amplitude no MWI (apenas interno)
qrs_i       = zeros(1,LLp);   % índice no MWI (apenas interno)
qrs_i_raw   = zeros(1,LLp);   % índice no bandpass (saída)
qrs_amp_raw = zeros(1,LLp);   % amplitude no bandpass (saída)

nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);

SIGL_buf  = zeros(1,LLp);
NOISL_buf = zeros(1,LLp);
SIGL_buf1  = zeros(1,LLp);
NOISL_buf1 = zeros(1,LLp);
THRS_buf  = zeros(1,LLp);
THRS_buf1 = zeros(1,LLp);

%% Training phase (2 s) para limiares
L2 = min(length(ecg_m), 2*fs);

THR_SIG   = max(ecg_m(1:L2)) * (1/3);
THR_NOISE = mean(ecg_m(1:L2)) * (1/2);
SIG_LEV   = THR_SIG;
NOISE_LEV = THR_NOISE;

THR_SIG1   = max(ecg_h(1:L2)) * (1/3);
THR_NOISE1 = mean(ecg_h(1:L2)) * (1/2);
SIG_LEV1   = THR_SIG1;
NOISE_LEV1 = THR_NOISE1;

%% ============ Thresholding and decision rule ============= %%
Beat_C  = 0;    % QRS detectados no MWI
Beat_C1 = 0;    % QRS confirmados no bandpass (saída efetiva)
Noise_Count = 0;

for i = 1:LLp

    %% localizar pico correspondente no bandpass ao redor de locs(i)
    if locs(i)-round(0.150*fs) >= 1 && locs(i) <= length(ecg_h)
        [y_i, x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
    else
        if i == 1
            [y_i, x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
        elseif locs(i) >= length(ecg_h)
            [y_i, x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
        else
            [y_i, x_i] = max(ecg_h(1:locs(i)));
        end
    end

    %% update RR / ajuste de threshold (igual ao original)
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

    %% search-back (igual ao original)
    if test_m
        if Beat_C > 0 && (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)
            lo = qrs_i(Beat_C) + round(0.200*fs);
            hi = locs(i) - round(0.200*fs);

            if lo >= 1 && hi > lo && hi <= length(ecg_m)
                [pks_temp, locs_temp_rel] = max(ecg_m(lo:hi));
                locs_temp = lo + locs_temp_rel - 1;

                if pks_temp > THR_NOISE
                    Beat_C = Beat_C + 1;
                    qrs_c(Beat_C) = pks_temp;
                    qrs_i(Beat_C) = locs_temp;

                    if locs_temp <= length(ecg_h)
                        seg_lo = max(1, locs_temp - round(0.150*fs));
                        seg_hi = locs_temp;
                        [y_i_t, x_i_t] = max(ecg_h(seg_lo:seg_hi));
                    else
                        seg_lo = max(1, locs_temp - round(0.150*fs));
                        [y_i_t, x_i_t] = max(ecg_h(seg_lo:end));
                    end

                    if y_i_t > THR_NOISE1
                        Beat_C1 = Beat_C1 + 1;
                        seg_lo = max(1, locs_temp - round(0.150*fs));
                        qrs_i_raw(Beat_C1)   = seg_lo + (x_i_t - 1);
                        qrs_amp_raw(Beat_C1) = y_i_t;
                        SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;
                    end

                    SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV;
                end
            end
        end
    end

    %% classificação de pico (QRS vs ruído) baseada no MWI
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

            if y_i >= THR_SIG1
                Beat_C1 = Beat_C1 + 1;
                if ser_back
                    qrs_i_raw(Beat_C1) = x_i;
                else
                    seg_lo = max(1, locs(i) - round(0.150*fs));
                    qrs_i_raw(Beat_C1) = seg_lo + (x_i - 1);
                end
                qrs_amp_raw(Beat_C1) = y_i;
                SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;
            end

            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV;
        end

    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;

    else % pks(i) < THR_NOISE
        Noise_Count = Noise_Count + 1;
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;
    end

    %% reajuste de thresholds (igual ao original)
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG   = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*THR_SIG;
    end
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1   = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*THR_SIG1;
    end

    SIGL_buf(i)   = SIG_LEV;
    NOISL_buf(i)  = NOISE_LEV;
    THRS_buf(i)   = THR_SIG;

    SIGL_buf1(i)  = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i)  = THR_SIG1;

    skip = 0;
    ser_back = 0;
end

%% Ajustar comprimentos finais (igual ao original)
qrs_i_raw   = qrs_i_raw(1:Beat_C1);
qrs_amp_raw = qrs_amp_raw(1:Beat_C1);

%% Plots de overlay (mantido, mas agora só com ecg_h/ecg_m)
if gr
    figure;

    subplot(211);
    plot(ecg_h); axis tight;
    title('QRS no Bandpass (ecg\_h)');
    hold on; scatter(qrs_i_raw, qrs_amp_raw, 'm');
    hold on; plot(locs, NOISL_buf1, '--k', 'LineWidth', 2);
    hold on; plot(locs, SIGL_buf1, '-.r', 'LineWidth', 2);
    hold on; plot(locs, THRS_buf1, '-.g', 'LineWidth', 2);

    subplot(212);
    plot(ecg_m); axis tight;
    title('QRS no MWI (ecg\_m) + níveis/limiar adaptativo');
    hold on; scatter(locs(pks >= THR_SIG), pks(pks >= THR_SIG), 'm'); % visual simples
    hold on; plot(locs, NOISL_buf, '--k', 'LineWidth', 2);
    hold on; plot(locs, SIGL_buf, '-.r', 'LineWidth', 2);
    hold on; plot(locs, THRS_buf, '-.g', 'LineWidth', 2);
end

end
