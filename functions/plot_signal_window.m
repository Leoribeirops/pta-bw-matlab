function plot_signal_window(x, plot_win, fig_id, fig_title)
%% plot_signal_window
% Plot genérico de sinal 1D com estilo padronizado
%
% Entradas:
%   x         : vetor do sinal (1D)
%   plot_win  : janela de plotagem em AMOSTRAS
%               0 ou []       -> plota o sinal inteiro
%               [ini fim]     -> intervalo em amostras
%   fig_id    : número da figura (ex: 1,2,3...)
%   fig_title : título do gráfico (string)
%
% Exemplo:
%   plot_signal_window(ecg_h, [5000 15000], 2, 'Bandpass')

%% ---------------- Validações ---------------- %%
if nargin < 4
    error('Uso: plot_signal_window(x, plot_win, fig_id, fig_title)');
end

if ~isvector(x)
    error('x deve ser um vetor 1D.');
end

x = x(:);
N = length(x);

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

%% ---------------- Plot ---------------- %%
figure(fig_id); clf;

if useWin
    idx = plot_win(1):plot_win(2);
else
    idx = 1:N;
end

plot(idx, x(idx), 'LineWidth', 1.0);
grid on; grid minor;

title(fig_title, 'FontSize', 12, 'FontWeight', 'bold');

xlabel('Samples [n]', 'Interpreter', 'latex');
ylabel('Amplitude',   'Interpreter', 'latex');

set(gcf, 'Color', 'w');
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'LineWidth', 0.8, ...
    'TickLabelInterpreter', 'latex');

end
