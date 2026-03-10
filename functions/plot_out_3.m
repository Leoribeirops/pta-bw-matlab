function janela = plot_out_3(fs, N, ecg_m, ecg_h, r_idx, trace, modo, j_ini, j_fim)
% plot_out_3 - Rotina de plotagem para ecg_m e ecg_h (PTA)
%
% modos:
%   0 -> plota apenas o sinal completo (full)
%   1 -> plota apenas a janela
%   2 -> plota janela + sinal completo
%
% Parâmetros opcionais:
%   j_ini, j_fim -> limites da janela (amostras). Se não fornecidos, usa default.

    %% ==============================
    %  1) Tratamento dos argumentos
    % ==============================
    if nargin < 7 || isempty(modo)
        modo = 1;  % default = janela
    end

    % Janela default, se não vier j_ini/j_fim
    if nargin >= 9 && ~isempty(j_ini) && ~isempty(j_fim)
        janela = [j_ini, min(j_fim, N)];
    else
        janela = [1000, min(8000, N)];  % ajuste livre
    end

    % Garante que início/fim estão dentro do intervalo
    j1 = max(janela(1), 1);
    j2 = min(janela(2), N);

    % Opcional: força início depois de 2*fs (como você fazia antes)
    j1 = max(j1, 2*fs);

    if j1 >= j2
        warning('Janela inválida. Ajustando para [2*fs, N].');
        j1 = max(2*fs, 1);
        j2 = N;
        janela = [j1, j2];
    else
        janela = [j1, j2];   % mantém janela real usada
    end

    idx_win  = j1:j2;  % índices da janela
    idx_full = 1:N;    % índices do sinal completo

    %% ==============================
    %  2) Tratamento ROBUSTO de r_idx
    % ==============================
    % Queremos no final: r_idx = [r1 r2 r3 ...] numérico, linha.

    if iscell(r_idx)
        % Remove células vazias
        r_idx = r_idx(~cellfun('isempty', r_idx));

        if isempty(r_idx)
            r_idx = [];
        elseif numel(r_idx) == 1 && isnumeric(r_idx{1})
            % Caso: um único cell com vetor numérico dentro
            r_idx = r_idx{1};
        else
            % Caso geral: várias células com escalares/vetores
            tmp = cellfun(@(x) x(:), r_idx, 'UniformOutput', false); % força colunas
            try
                r_idx = vertcat(tmp{:}); % empilha tudo
            catch ME
                warning('Falha ao concatenar r_idx (cell). Detalhe: %s', ME.message);
                % Debug opcional:
                % celldisp(r_idx);
                r_idx = [];  % para evitar crash; não plota R depois
            end
        end
    end

    % Garante que agora é numérico
    if ~isnumeric(r_idx)
        warning('r_idx não é numérico após tratamento. Ignorando R-peaks.');
        r_idx = [];
    end

    r_idx = double(r_idx(:)');      % vetor linha
    r_idx = r_idx(~isnan(r_idx));   % remove NaNs
    r_idx = r_idx(r_idx >= 1 & r_idx <= N); % clamp dentro do sinal

    % Índices válidos de R-peaks
    valid_idx_win  = r_idx(r_idx >= j1 & r_idx <= j2);
    valid_idx_full = r_idx;  % já filtrado para [1, N]

    %% ==============================
    %  3) MODO 1 ou 2 : Plotar JANELA
    % ==============================
    if modo == 1 || modo == 2
        figure('Name','Janela - ecg_m e ecg_h');
        tlo = tiledlayout(2,1);
        title(tlo, sprintf('Janela [%d, %d] amostras', j1, j2));

        % ---- Subplot 1: Integrador (ecg_m) ----
        nexttile;
        plot(idx_win, ecg_m(idx_win), 'LineWidth', 1); hold on; grid on; grid minor
        plot(idx_win, trace.NOISE_LEV(idx_win),'k--','LineWidth',1.0);
        plot(idx_win, trace.SIG_LEV(idx_win),  'r-.','LineWidth',1.0);
        plot(idx_win, trace.THR_SIG(idx_win),  'g-','LineWidth',1.0);
        if ~isempty(valid_idx_win)
            scatter(valid_idx_win, ecg_m(valid_idx_win), 25, 'm', 'filled');
        end
        legend('ecg_m','NOISE','SIG','THR','R','Location','best');
        ylabel('Amp');
        title('Integrador (ecg_m) - janela');

        % ---- Subplot 2: Bandpass (ecg_h) ----
        nexttile;
        plot(idx_win, ecg_h(idx_win), 'LineWidth', 1); hold on; grid on; grid minor
        plot(idx_win, trace.NOISE_LEV1(idx_win),'k--','LineWidth',1.0);
        plot(idx_win, trace.SIG_LEV1(idx_win),  'r-.','LineWidth',1.0);
        plot(idx_win, trace.THR_SIG1(idx_win),  'g-','LineWidth',1.0);
        if ~isempty(valid_idx_win)
            scatter(valid_idx_win, ecg_h(valid_idx_win), 25, 'm', 'filled');
        end
        legend('ecg_h','NOISE1','SIG1','THR1','R','Location','best');
        xlabel('Amostras'); ylabel('Amp');
        title('Bandpass (ecg_h) - janela');
    end

    %% ==============================
    %  4) MODO 0 ou 2 : Plotar SINAL COMPLETO
    % ==============================
    if modo == 0 || modo == 2
        % -------- Figura 1: ecg_h completo --------
        figure('Name','Sinal completo - ecg_h');
        plot(idx_full, ecg_h, 'LineWidth', 1); hold on; grid on; grid minor
        plot(idx_full, trace.NOISE_LEV1, 'k--', 'LineWidth', 1.0);
        plot(idx_full, trace.SIG_LEV1,   'r-.', 'LineWidth', 1.0);
        plot(idx_full, trace.THR_SIG1,   'g-',  'LineWidth', 1.0);
        if ~isempty(valid_idx_full)
            scatter(valid_idx_full, ecg_h(valid_idx_full), 15, 'm', 'filled');
        end
        legend('ecg_h','NOISE1','SIG1','THR1','R','Location','best');
        title('Sinal completo: ecg_h + limiares');
        xlabel('Amostras'); ylabel('Amplitude');

        % -------- Figura 2: ecg_m completo --------
        figure('Name','Sinal completo - ecg_m');
        plot(idx_full, ecg_m, 'LineWidth', 1); hold on; grid on; grid minor
        plot(idx_full, trace.NOISE_LEV, 'k--', 'LineWidth', 1.0);
        plot(idx_full, trace.SIG_LEV,   'r-.', 'LineWidth', 1.0);
        plot(idx_full, trace.THR_SIG,   'g-',  'LineWidth', 1.0);
        if ~isempty(valid_idx_full)
            scatter(valid_idx_full, ecg_m(valid_idx_full), 15, 'm', 'filled');
        end
        legend('ecg_m','NOISE','SIG','THR','R','Location','best');
        title('Sinal completo: ecg_m + limiares');
        xlabel('Amostras'); ylabel('Amplitude');
    end

end
