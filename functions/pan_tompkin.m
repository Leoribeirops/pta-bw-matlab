function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs,gr,winN,i0)
%% function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs,gr,winN,i0)
% Complete implementation of Pan-Tompkins algorithm (Sedghamiz)
%
% Inputs
%   ecg  : raw ecg vector signal (1D)
%   fs   : sampling frequency (Hz)
%   gr   : 1 plots / 0 no plots  (default = 1)
%   winN : (optional) window length in samples for visualization (default = 6000)
%   i0   : (optional) window start sample for visualization (default = 1)
%
% Outputs
%   qrs_amp_raw : amplitude of R waves amplitudes (bandpass)
%   qrs_i_raw   : index of R waves (bandpass)
%   delay       : number of samples which the signal is delayed due to filtering
%
% Ref:
% https://www.researchgate.net/publication/313673153_Matlab_Implementation_of_Pan_Tompkins_ECG_QRS_detector

if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
if nargin < 3
    gr = 1;
end
if nargin < 4 || isempty(winN)
    winN = 6000;
end
if nargin < 5 || isempty(i0)
    i0 = 1;
end

ecg = ecg(:); % vectorize

% ---- visualization window indices (based on raw ECG length) ----
Nsig = length(ecg);
i0 = max(1, min(i0, Nsig));
i1 = min(Nsig, i0 + winN - 1);

%% ======================= Initialize =============================== %
delay = 0;
skip = 0;                                                                  % becomes one when a T wave is detected
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0;
ax = zeros(1,6);

%% ============ Noise cancelation(Filtering)( 5-15 Hz) =============== %%
if fs == 200
% ------------------ remove the mean of Signal -----------------------%
  ecg = ecg - mean(ecg);
%% Low Pass (Butterworth replacement)
   Wn = 12*2/fs;
   N = 3;                                                                  % order of 3 less processing
   [a,b] = butter(N,Wn,'low');
   ecg_l = filtfilt(a,b,ecg);
   ecg_l = ecg_l/ max(abs(ecg_l));

   if gr
    figure;
    ax(1) = subplot(321);plot(ecg);axis tight;title('Raw signal');
    ax(2) = subplot(322);plot(ecg_l);axis tight;title('Low pass filtered');
   end

%% High Pass (Butterworth replacement)
   Wn = 5*2/fs;
   N = 3;
   [a,b] = butter(N,Wn,'high');
   ecg_h = filtfilt(a,b,ecg_l);
   ecg_h = ecg_h/ max(abs(ecg_h));

   if gr
    ax(3) = subplot(323);plot(ecg_h);axis tight;title('High Pass Filtered');
   end
else
%% bandpass for other fs
 f1=5;
 f2=15;
 Wn=[f1 f2]*2/fs;
 N = 3;
 [a,b] = butter(N,Wn);
 ecg_h = filtfilt(a,b,ecg);
 ecg_h = ecg_h/ max(abs(ecg_h));

 if gr
  ax(1) = subplot(3,2,[1 2]);plot(ecg);axis tight;title('Raw Signal');
  ax(3) = subplot(323);plot(ecg_h);axis tight;title('Band Pass Filtered');
 end
end

%% ==================== derivative filter ========================== %%
% H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
if fs ~= 200
 int_c = (5-1)/(fs*1/40);
 b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
 b = [1 2 0 -2 -1].*(1/8)*fs;
end
ecg_d = filtfilt(b,1,ecg_h);
ecg_d = ecg_d/max(ecg_d);

if gr
  ax(4)=subplot(324);plot(ecg_d); axis tight; title('Filtered with the derivative filter');
end

%% ========== Squaring nonlinearly enhance the dominant peaks ========== %%
ecg_s = ecg_d.^2;
if gr
  ax(5)=subplot(325); plot(ecg_s); axis tight; title('Squared');
end

%% ============  Moving average ================== %%
% Y(nt) = (1/N) sum_{k=0..N-1} x(n-k)
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + round(0.150*fs)/2;

if gr
  ax(6)=subplot(326);plot(ecg_m); axis tight;
  title('Averaged (MWI), Black=noise, Green=threshold, Red=sig level, Magenta=QRS');
end

%% -------- apply visualization window to the "classic 6-subplot" figure -----
if gr && any(ax)
    ax(~ax) = [];
    linkaxes(ax,'x'); zoom on;
    xlim([i0 i1]);
end

%% ===================== Fiducial Marks ============================== %%
[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*fs));

%% =================== Initialize Some Other Parameters =============== %%
LLp = length(pks);

qrs_c = zeros(1,LLp);           % amplitude of R (MWI)
qrs_i = zeros(1,LLp);           % index (MWI)
qrs_i_raw = zeros(1,LLp);       % index (bandpass)
qrs_amp_raw= zeros(1,LLp);      % amplitude (bandpass)

nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);

SIGL_buf  = zeros(1,LLp);
NOISL_buf = zeros(1,LLp);
SIGL_buf1 = zeros(1,LLp);
NOISL_buf1 = zeros(1,LLp);
THRS_buf1 = zeros(1,LLp);
THRS_buf  = zeros(1,LLp);

%% training phase (2 seconds)
THR_SIG   = max(ecg_m(1:2*fs))*1/3;
THR_NOISE = mean(ecg_m(1:2*fs))*1/2;
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;

THR_SIG1   = max(ecg_h(1:2*fs))*1/3;
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2;
SIG_LEV1 = THR_SIG1;
NOISE_LEV1 = THR_NOISE1;

%% ============ Thresholding and decision rule ============= %%
Beat_C = 0;
Beat_C1 = 0;
Noise_Count = 0;

for i = 1 : LLp

   %% locate corresponding peak in filtered signal (ecg_h)
    if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_h)
          [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
       else
          if i == 1
            [y_i,x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
          elseif locs(i)>= length(ecg_h)
            [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
          end
    end

  %% update heart rate estimate
    if Beat_C >= 9
        diffRR = diff(qrs_i(Beat_C-8:Beat_C));
        mean_RR = mean(diffRR);
        comp = qrs_i(Beat_C)-qrs_i(Beat_C-1);

        if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
            THR_SIG  = 0.5*(THR_SIG);
            THR_SIG1 = 0.5*(THR_SIG1);
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

    %% search-back
    if test_m
        if (locs(i) - qrs_i(max(1,Beat_C))) >= round(1.66*test_m) && Beat_C>0
            [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs)));
            locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1;

            if pks_temp > THR_NOISE
                Beat_C = Beat_C + 1;
                qrs_c(Beat_C) = pks_temp;
                qrs_i(Beat_C) = locs_temp;

                if locs_temp <= length(ecg_h)
                    [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):locs_temp));
                else
                    [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):end));
                end

                if y_i_t > THR_NOISE1
                    Beat_C1 = Beat_C1 + 1;
                    qrs_i_raw(Beat_C1) = locs_temp-round(0.150*fs)+ (x_i_t - 1);
                    qrs_amp_raw(Beat_C1) = y_i_t;
                    SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;
                end

                SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV;
            end
        end
    end

    %% classify peaks: QRS vs noise (MWI domain)
    if pks(i) >= THR_SIG

       % T-wave discrimination
       if Beat_C >= 3
          if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs)
              Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i))));
              Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C))));
              if abs(Slope1) <= abs(0.5*(Slope2))
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
                 qrs_i_raw(Beat_C1)= locs(i)-round(0.150*fs)+ (x_i - 1);
              end
              qrs_amp_raw(Beat_C1) =  y_i;
              SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;
          end
          SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV;
        end

    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
         NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
         NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;
    elseif pks(i) < THR_NOISE
        Noise_Count = Noise_Count + 1;
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV  = 0.125*pks(i) + 0.875*NOISE_LEV;
    end

    %% adjust thresholds (SNR)
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end

    %% buffers
    SIGL_buf(i)  = SIG_LEV;
    NOISL_buf(i) = NOISE_LEV;
    THRS_buf(i)  = THR_SIG;

    SIGL_buf1(i)  = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i)  = THR_SIG1;

    %% reset
    skip = 0;
    ser_back = 0;
end

%% ======================= Adjust Lengths ============================ %%
qrs_i_raw   = qrs_i_raw(1:Beat_C1);
qrs_amp_raw = qrs_amp_raw(1:Beat_C1);
qrs_c = qrs_c(1:Beat_C);
qrs_i = qrs_i(1:Beat_C);

%% ======================= Plottings (MWI figure) ===================== %%
if gr
  % janela aplicada (MWI domain)
  inW_mwi = (qrs_i >= i0) & (qrs_i <= i1);
  qrs_i_w = qrs_i(inW_mwi);
  qrs_c_w = qrs_c(inW_mwi);

  inL = (locs >= i0) & (locs <= i1);
  locs_w = locs(inL);

  hold on; scatter(qrs_i_w,qrs_c_w,'m');
  hold on; plot(locs_w,NOISL_buf(inL),'--k','LineWidth',2);
  hold on; plot(locs_w,SIGL_buf(inL),'--r','LineWidth',2);
  hold on; plot(locs_w,THRS_buf(inL),'--g','LineWidth',2);

  if any(ax)
    ax(~ax) = [];
    linkaxes(ax,'x');
    zoom on;
  end
  xlim([i0 i1]);
end

%% ================== overlay on the signals (3 subplots) ============= %%
if gr
   figure;

   % indices within window
   inW_raw = (qrs_i_raw >= i0) & (qrs_i_raw <= i1);
   qrs_i_raw_w   = qrs_i_raw(inW_raw);
   qrs_amp_raw_w = qrs_amp_raw(inW_raw);

   inW_mwi = (qrs_i >= i0) & (qrs_i <= i1);
   qrs_i_w = qrs_i(inW_mwi);
   qrs_c_w = qrs_c(inW_mwi);

   inL = (locs >= i0) & (locs <= i1);
   locs_w = locs(inL);

   az(1)=subplot(311);
   plot(ecg_h);
   title('QRS on Filtered Signal');
   axis tight;
   hold on; scatter(qrs_i_raw_w,qrs_amp_raw_w,'m');
   hold on; plot(locs_w,NOISL_buf1(inL),'LineWidth',2,'Linestyle','--','color','k');
   hold on; plot(locs_w,SIGL_buf1(inL),'LineWidth',2,'Linestyle','-.','color','r');
   hold on; plot(locs_w,THRS_buf1(inL),'LineWidth',2,'Linestyle','-.','color','g');
   xlim([i0 i1]);

   az(2)=subplot(312);
   plot(ecg_m);
   title('QRS on MVI signal and Noise level(black),Signal Level (red) and Adaptive Threshold(green)');
   axis tight;
   hold on; scatter(qrs_i_w,qrs_c_w,'m');
   hold on; plot(locs_w,NOISL_buf(inL),'LineWidth',2,'Linestyle','--','color','k');
   hold on; plot(locs_w,SIGL_buf(inL),'LineWidth',2,'Linestyle','-.','color','r');
   hold on; plot(locs_w,THRS_buf(inL),'LineWidth',2,'Linestyle','-.','color','g');
   xlim([i0 i1]);

   az(3)=subplot(313);
   ecg0 = ecg-mean(ecg);
   plot(ecg0);
   title('Pulse train of the found QRS on ECG signal');
   axis tight;

   % pulse train only inside window
   if ~isempty(qrs_i_raw_w)
       ylo = min(ecg0(i0:i1))/2;
       yhi = max(ecg0(i0:i1))/2;
       line(repmat(qrs_i_raw_w,[2 1]), repmat([ylo; yhi], size(qrs_i_raw_w)), ...
           'LineWidth',2.5,'LineStyle','-.','Color','r');
   end
   xlim([i0 i1]);

   linkaxes(az,'x');
   zoom on;
end

end
