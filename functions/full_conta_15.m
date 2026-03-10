function [a, missed_count, recovered_count, r_idx, trace] = full_conta_15(fs, N, PROFILE, ecg_h, ecg_m)
% FULL_CONTA_10
% Versão derivada da full_conta_08/09 com:
%   - Piso suave para SIG_LEV/SIG_LEV1
%   - Suavização de THR_SIG1 (h) em 0.5/0.5
%   - Suavização mais forte de THR_SIG (m) em 0.8/0.2
%
% Entradas:
%   fs, N, PROFILE, ecg_h, ecg_m (mesmo formato das versões anteriores)

% =====================================================================
% 1) Perfil
% =====================================================================
switch PROFILE
    case 'classic'
        ALPHA_THR_Q       = 0.25;
        SEARCH_REL        = 0.50;
        RR_TIMEOUT_FACTOR = 1.66;
    case 'modified'
        ALPHA_THR_Q       = 0.189;
        SEARCH_REL        = 0.30;
        RR_TIMEOUT_FACTOR = 1.50;
    otherwise
        error('PROFILE inválido.');
end

% =====================================================================
% 2) Constantes
% =====================================================================
DELAY  = 54;
RR_MIN = round(0.200*fs);

HYST_SAMPLES = round(0.05*fs);
HYST_GAIN    = 0.125;

EMA_COEF = 1/8;

EMA_NOISE_UP   = 1/16;
EMA_NOISE_DOWN = EMA_COEF;

MAX_NOISE_STEP_FRAC = 1/16;

NOISE_MAX_FRACTION  = 0.75;
NOISE_HARD_FRACTION = 0.50;

FALLBACK_RR_MULT = 2.0;
FALLBACK_GAIN    = 0.50;
MIN_THR_FRAC     = 0.15;

thr_scale_m   = 1.0;
thr_scale_h   = 1.0;
MIN_THR_SCALE = 0.20;
MAX_THR_SCALE = 1.80;

% =====================================================================
% 3) Treinamento inicial (0–2 s)
% =====================================================================
ini_tr = 1;
fim_tr = 2*fs;
assert(fim_tr <= N,'Sinal curto demais para treino');

trace.SIG_LEV   = zeros(N,1);
trace.NOISE_LEV = zeros(N,1);
trace.THR_SIG   = zeros(N,1);

trace.SIG_LEV1   = zeros(N,1);
trace.NOISE_LEV1 = zeros(N,1);
trace.THR_SIG1   = zeros(N,1);

% SIG/NOISE (m)
[vmax_m, ~] = max_in_range(ecg_m, ini_tr, fim_tr);
SIG_LEV = 0.25*vmax_m;

sum_m = 0;
len_tr = fim_tr - ini_tr + 1;
for k = ini_tr:fim_tr
    sum_m = sum_m + ecg_m(k);
end
NOISE_LEV = 0.25*(sum_m/len_tr);

% SIG/NOISE (h)
[vmax_h, ~] = max_in_range(ecg_h, ini_tr, fim_tr);
SIG_LEV1 = 0.25*vmax_h;

sum_h = 0;
for k = ini_tr:fim_tr
    sum_h = sum_h + ecg_h(k);
end
NOISE_LEV1 = 0.25*(sum_h/len_tr);

% ---------- PATCH 1: pisos suaves de SIG ----------
SIG_INIT  = SIG_LEV;
SIG_INIT1 = SIG_LEV1;

SIG_FLOOR  = SIG_INIT  * (1/64);   % piso ~1,56% do inicial
SIG1_FLOOR = SIG_INIT1 * (1/64);

% ---------- PATCH 2: registradores para suavização de THR ----------
THR_SIG_prev  = 0;
THR_SIG1_prev = 0;

THR_SIG = 0; THR_SIG1 = 0; THR_NOISE = 0; THR_NOISE1 = 0;
recalc_thresholds();

% =====================================================================
% 4) Estados / buffers
% =====================================================================
IDLE = 0; ARM = 1; EVAL = 2;
state = IDLE;

valor_arm   = 0; pos_m = 0;
valor_arm_h = 0; pos_h = 0;

hyst_count = 0;

r_idx   = zeros(N,1);
r_amp_h = zeros(N,1);
a = 2;

RR_buf   = zeros(8,1);
rr_count = 0;
last_R   = 0;

missed_count    = 0;
recovered_count = 0;
timeout_processed = false;

% =====================================================================
% 5) Loop principal
% =====================================================================
i_start = fim_tr + 1;

for i = i_start:N

    % ----------------- FSM principal -----------------
    switch state
        case IDLE
            if ecg_m(i) > NOISE_LEV
                state     = ARM;
                valor_arm = ecg_m(i);
                pos_m     = i;

                bs_from = max(i-DELAY+1,1);
                [valor_arm_h,pos_h] = max_in_range(ecg_h,bs_from,i);
            end

        case ARM
            if ecg_m(i) > valor_arm
                valor_arm = ecg_m(i);
                pos_m     = i;

                bs_from = max(i-DELAY+1,1);
                [valor_arm_h,pos_h] = max_in_range(ecg_h,bs_from,i);
            else
                state = EVAL;
            end

        case EVAL
            accept = false;
            rr_last_for_thresh = 0;

            if (last_R==0) || ((pos_h-last_R) >= RR_MIN)

                pass_t = pass_slope_and_width(ecg_m,pos_m,fs,last_R);

                if pass_t && ...
                   (valor_arm   >= THR_SIG)  && ...
                   (valor_arm_h >= THR_SIG1) && ...
                   is_local_max(ecg_h,pos_h)

                    % --------- ACEITO ----------
                    a = a+1;
                    r_idx(a)   = pos_h;
                    r_amp_h(a) = valor_arm_h;
                    accept = true;

                    SIG_LEV  = SIG_LEV  + EMA_COEF*(valor_arm   - SIG_LEV);
                    SIG_LEV1 = SIG_LEV1 + EMA_COEF*(valor_arm_h - SIG_LEV1);

                    % patch 1: piso
                    SIG_LEV  = max(SIG_LEV , SIG_FLOOR);
                    SIG_LEV1 = max(SIG_LEV1, SIG1_FLOOR);

                    if last_R>0
                        rr_last_for_thresh = pos_h - last_R;
                        rr_count = min(rr_count+1,8);
                        RR_buf   = circ_push(RR_buf,rr_last_for_thresh);
                    end
                    last_R = pos_h;

                    hyst_count = HYST_SAMPLES;
                    timeout_processed = false;

                else
                    % --------- REJEITADO ----------
                    if (valor_arm < THR_NOISE) && (valor_arm_h < THR_NOISE1)
                        % canal m
                        coef = EMA_NOISE_DOWN;
                        if valor_arm > NOISE_LEV, coef = EMA_NOISE_UP; end
                        cand = NOISE_LEV + coef*(valor_arm - NOISE_LEV);
                        cand = min(cand, NOISE_LEV + MAX_NOISE_STEP_FRAC*SIG_LEV);
                        NOISE_LEV = min(cand, NOISE_MAX_FRACTION*SIG_LEV);

                        % canal h
                        coef = EMA_NOISE_DOWN;
                        if valor_arm_h > NOISE_LEV1, coef = EMA_NOISE_UP; end
                        cand = NOISE_LEV1 + coef*(valor_arm_h - NOISE_LEV1);
                        cand = min(cand, NOISE_LEV1 + MAX_NOISE_STEP_FRAC*SIG_LEV1);
                        NOISE_LEV1 = min(cand, NOISE_MAX_FRACTION*SIG_LEV1);
                    end
                end

                % ------ Ajuste de escala por RR ------
                if rr_last_for_thresh>0 && rr_count>=8
                    rr_mean = mean_nonzero(RR_buf);
                    if rr_mean>0
                        if (rr_last_for_thresh < 0.92*rr_mean) || ...
                           (rr_last_for_thresh > 1.16*rr_mean)
                            thr_scale_m = max(thr_scale_m*0.5, MIN_THR_SCALE);
                            thr_scale_h = max(thr_scale_h*0.5, MIN_THR_SCALE);
                        end
                    end
                end

                % Histerese só se aceitou
                if accept && hyst_count>0
                    thr_scale_m = min(thr_scale_m*(1+HYST_GAIN),  MAX_THR_SCALE);
                    thr_scale_h = min(thr_scale_h*(1+HYST_GAIN),  MAX_THR_SCALE);
                end

                recalc_thresholds();
            end

            valor_arm   = 0;
            valor_arm_h = 0;
            state = IDLE;
    end

    % ----------------- Histerese -----------------
    if hyst_count>0
        hyst_count = hyst_count - 1;
    end

    % ----------------- Timeout / search-back -----------------
    if rr_count>=3
        RR_med = median_nonzero(RR_buf);
        RR_to  = round(RR_TIMEOUT_FACTOR * RR_med);
    else
        RR_to  = 0;
    end

    if (last_R>0) && (rr_count>=3) && (RR_to>0) && ((i-last_R)>=RR_to) && ~timeout_processed
        timeout_processed = true;
        missed_count = missed_count + 1;

        lo = max(last_R + RR_MIN, 1);
        hi = max(i - RR_MIN, lo);

        if hi>lo
            [best_v,best_i] = max_in_range(ecg_m,lo,hi);

            bs_from = max(best_i-DELAY+1,1);
            [best_h,best_h_i] = max_in_range(ecg_h,bs_from,best_i);

            if (best_v > SEARCH_REL*THR_SIG) && ...
               (best_h > SEARCH_REL*THR_SIG1) && ...
               is_local_max(ecg_h,best_h_i)

                % ---- RECUPERADO ----
                a = a+1;
                r_idx(a)   = best_h_i;
                r_amp_h(a) = best_h;

                SIG_LEV  = SIG_LEV  + EMA_COEF*(best_v - SIG_LEV);
                SIG_LEV1 = SIG_LEV1 + EMA_COEF*(best_h - SIG_LEV1);

                SIG_LEV  = max(SIG_LEV , SIG_FLOOR);
                SIG_LEV1 = max(SIG_LEV1, SIG1_FLOOR);

                rr_last = best_h_i - last_R;
                rr_count = min(rr_count+1,8);
                RR_buf   = circ_push(RR_buf, rr_last);

                last_R = best_h_i;
                recovered_count = recovered_count+1;

                hyst_count = round(0.03*fs);

                if rr_count>=8
                    rr_mean = mean_nonzero(RR_buf);
                    if rr_mean>0
                        if (rr_last < 0.92*rr_mean) || (rr_last > 1.16*rr_mean)
                            thr_scale_m = max(thr_scale_m*0.5, MIN_THR_SCALE);
                            thr_scale_h = max(thr_scale_h*0.5, MIN_THR_SCALE);
                        end
                    end
                end

                thr_scale_m = min(thr_scale_m*(1+HYST_GAIN/2), MAX_THR_SCALE);
                thr_scale_h = min(thr_scale_h*(1+HYST_GAIN/2), MAX_THR_SCALE);

                recalc_thresholds();
            else
                if best_v < THR_NOISE
                    coef = EMA_NOISE_DOWN;
                    if best_v > NOISE_LEV, coef = EMA_NOISE_UP; end
                    cand = NOISE_LEV + coef*(best_v - NOISE_LEV);
                    cand = min(cand, NOISE_LEV + MAX_NOISE_STEP_FRAC*SIG_LEV);
                    NOISE_LEV = min(cand, NOISE_MAX_FRACTION*SIG_LEV);
                end
                recalc_thresholds();
            end
        end
    end

    % ----------------- Fallback global -----------------
    if last_R>0
        if rr_count>=3 && RR_to>0
            limit = round(FALLBACK_RR_MULT * RR_to);
        else
            limit = round(2.0*fs);
        end

        if (i-last_R) > limit
            thr_scale_m = max(thr_scale_m*FALLBACK_GAIN, MIN_THR_SCALE);
            thr_scale_h = max(thr_scale_h*FALLBACK_GAIN, MIN_THR_SCALE);

            SIG_LEV  = SIG_LEV  * 0.995;
            SIG_LEV1 = SIG_LEV1 * 0.995;

            SIG_LEV  = max(SIG_LEV , SIG_FLOOR);
            SIG_LEV1 = max(SIG_LEV1, SIG1_FLOOR);

            recalc_thresholds();
        end
    end

    % ----------------- Clamp de NOISE -----------------
    if SIG_LEV>0 && NOISE_LEV > NOISE_HARD_FRACTION*SIG_LEV
        NOISE_LEV = NOISE_HARD_FRACTION*SIG_LEV;
    end
    if SIG_LEV1>0 && NOISE_LEV1 > NOISE_HARD_FRACTION*SIG_LEV1
        NOISE_LEV1 = NOISE_HARD_FRACTION*SIG_LEV1;
    end

    recalc_thresholds();

    % ----------------- Trace -----------------
    trace.SIG_LEV(i)   = SIG_LEV;
    trace.NOISE_LEV(i) = NOISE_LEV;
    trace.THR_SIG(i)   = THR_SIG;

    trace.SIG_LEV1(i)   = SIG_LEV1;
    trace.NOISE_LEV1(i) = NOISE_LEV1;
    trace.THR_SIG1(i)   = THR_SIG1;
end

r_idx = r_idx(1:a);

% =====================================================================
% Funções auxiliares
% =====================================================================
    function recalc_thresholds()
        % base
        base_m = NOISE_LEV  + ALPHA_THR_Q*(SIG_LEV  - NOISE_LEV);
        base_h = NOISE_LEV1 + ALPHA_THR_Q*(SIG_LEV1 - NOISE_LEV1);

        raw_m = thr_scale_m * base_m;
        raw_h = thr_scale_h * base_h;

        % ---- PATCH: suavização diferente para m e h ----
        % m (ecg_m): filtro mais pesado (0.8 / 0.2)
        if THR_SIG_prev == 0
            THR_SIG = raw_m;
        else
            THR_SIG = 0.8*THR_SIG_prev + 0.2*raw_m;
        end

        % h (ecg_h): filtro 0.5 / 0.5 (mais ágil)
        if THR_SIG1_prev == 0
            THR_SIG1 = raw_h;
        else
            THR_SIG1 = 0.5*THR_SIG1_prev + 0.5*raw_h;
        end

        % pisos
        if SIG_LEV>0
            THR_SIG = max(THR_SIG, MIN_THR_FRAC*SIG_LEV);
        end
        if SIG_LEV1>0
            THR_SIG1 = max(THR_SIG1, MIN_THR_FRAC*SIG_LEV1);
        end

        THR_NOISE  = 0.5*THR_SIG;
        THR_NOISE1 = 0.5*THR_SIG1;

        THR_SIG_prev  = THR_SIG;
        THR_SIG1_prev = THR_SIG1;
    end

    function [vmax,imax] = max_in_range(x,i1,i2)
        vmax = x(i1); imax = i1;
        for k=i1+1:i2
            if x(k)>vmax
                vmax = x(k); imax = k;
            end
        end
    end

    function ok = is_local_max(x,idx)
        L = length(x);
        if idx<=1 || idx>=L
            ok=false; return;
        end
        ok = (x(idx)>x(idx-1)) && (x(idx)>x(idx+1));
    end

    function m = median_nonzero(buf)
        nz = buf(buf>0);
        if isempty(nz), m=0; return; end
        nz = sort(nz);
        m  = nz(floor((length(nz)+1)/2));
    end

    function m = mean_nonzero(buf)
        nz = buf(buf>0);
        if isempty(nz), m=0; else m = mean(nz); end
    end

    function b = circ_push(buf,val)
        b = [buf(2:end); val];
    end

    function pass = pass_slope_and_width(m_sig,pos_m_loc,fs_loc,last_R_loc)
        W = max(round(0.075*fs_loc),1);
        pass = true;

        if pos_m_loc <= W || (last_R_loc>0 && last_R_loc<=W)
            return;
        end

        s1 = 0;
        for t = (pos_m_loc-W+1):pos_m_loc
            s1 = s1 + (m_sig(t) - m_sig(t-1));
        end
        s1 = s1/W;

        if last_R_loc>0
            s2 = 0;
            for t = (last_R_loc-W+1):last_R_loc
                s2 = s2 + (m_sig(t) - m_sig(t-1));
            end
            s2 = s2/W;
        else
            s2 = s1;
        end

        if abs(s1) <= 0.5*abs(s2)
            pass = false; return;
        end

        half = 0.5*m_sig(pos_m_loc);
        max_half_width = round(0.200*fs_loc);
        Lm = length(m_sig);

        le = pos_m_loc;
        for k=1:max_half_width
            if le<=1, break; end
            if m_sig(le) <= half, break; end
            le = le - 1;
        end

        ri = pos_m_loc;
        for k=1:max_half_width
            if ri>=Lm, break; end
            if m_sig(ri) <= half, break; end
            ri = ri + 1;
        end

        width = ri - le;
        if width < round(0.060*fs_loc) || width > round(0.200*fs_loc)
            pass = false;
        end
    end

end
