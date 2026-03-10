function filt = panTompkinsFilterBW(ecg)
% panTompkinsFilter
% Executa as etapas de filtragem do pipeline Pan-Tompkins:
% bandpass -> derivada -> quadrado -> MWI
%
% Entrada:
%   ecg  : vetor do sinal ECG
%
% Saída:
%   filt : struct com os sinais intermediários

    %% Garantir vetor coluna
    %ecg = ecg(:);
    ecg = ecg(:)';
    %% Bandpass
    b_bw = [1/16, 0, -1/8, 0, 1/16];
    a_bw = [1, -3.176, 3.808, -2.076, 0.444];
    bw_out = filter(b_bw, a_bw, ecg);

    %% Derivative
    b_diff = [1/4, 1/8, 0, -1/8, -1/4];
    a_diff = 1;
    diff_out = filter(b_diff, a_diff, bw_out);

    %% Squaring
    square_out = diff_out.^2;

    %% Moving Window Integration
    N = 30;
    b_mwi = (1/N) * ones(1, N+1);
    a_mwi = 1;
    mwi_out = filter(b_mwi, a_mwi, square_out);

    %% Organizar saída
    filt = struct();
    filt.ecg = ecg;
    filt.bandpass = bw_out;
    filt.derivative = diff_out;
    filt.square = square_out;
    filt.mwi = mwi_out;

    % opcional: guardar coeficientes usados
    filt.params.b_bw = b_bw;
    filt.params.a_bw = a_bw;
    filt.params.b_diff = b_diff;
    filt.params.a_diff = a_diff;
    filt.params.N_mwi = N;
    filt.params.b_mwi = b_mwi;
    filt.params.a_mwi = a_mwi;
end