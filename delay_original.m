clear; clc; close all;

% load simulation environment data
load('sim_results\supporting_simdata.mat');
% load('sim_results\supporting_simdata2.mat');

% load simulation a-mode data
theta     = 0:-2.5:-20;
% theta     = -5;
n_element = 3; 
tx_idcs   = 1:n_element;

% all_data is a [M,N] cell matrix, 
% M -> angle, N-> emmitting transducer
all_data = {};
for i=1:length(theta)
    for j=tx_idcs
        str_matname = sprintf('signal_%s_%d.mat', num2str(theta(i), '%+03.f'), j);
        str_matpath = fullfile(pwd, 'sim_results', str_matname);
        load(str_matpath);
        all_data{i, j} = transducers;
    end
end
clear transducers str_matname str_matpath;

% [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20] degree;
select_idx = 5;

% FMC_singalall is a [n_Tx, n_Rx] matrix
FMC_singalall = cat(1, all_data{select_idx, tx_idcs});

% FMC_signalraw will be [n_sample, n_tx, n_rx];
FMC_signal = cat(1, FMC_singalall(:,:).signal_corr);
FMC_signal = reshape(FMC_signal', [], 3, 3);

%% Setup transducer configuration
% All of these configuration below known from the simulation. Later, if you
% want another simulation, don't forget to change this configuration

transducer_cfg.element.n     = 3;
transducer_cfg.element.width = 6 * 1e-3;   % [m]
transducer_cfg.element.kerf  = 2 * 1e-3;   % [m]
transducer_cfg.element.pitch = 8 * 1e-3;   % [m]
transducer_cfg.period        = kgrid.dt;   % [m]
transducer_cfg.freq          = 1/kgrid.dt; % [s]
transducer_cfg.t_array       = kgrid.t_array; % [[s]]
transducer_cfg.numsamples    = kgrid.Nt;

%%

% select us element
tx = 2;
rx = 3;

% discard near field disturbance
cutoff_s   = 2 * 1e-6; % [s]
cutoff_idx = round(cutoff_s/kgrid.dt);

% get the signal
s1 = FMC_signal(cutoff_idx:end,tx,tx);
s2 = FMC_signal(cutoff_idx:end,tx,rx);
t_arraycut = transducer_cfg.t_array(cutoff_idx:end);

% normalize the signal
s1_norm = s1 / max(s1);
s2_norm = s2 / max(s2);

% align the signal
maxlag_s   = 3 * 1e-6;
maxlag_idx = round(maxlag_s/kgrid.dt);

% let's consider first param is the ref
% [~, ~, D] = alignsignals(s1_norm, s2_norm, 'Method', 'xcorr', 'MaxLag', maxlag_idx);
d = finddelay(s1_norm, s2_norm, maxlag_idx);
% if d>0, y is delayed, let's bring them to x
if(d>0)
    s1_new = s1;
    s2_new = [s2(d+1:end); zeros(d, 1)];
% if d<0, x is delayed
elseif (d<0)
    D = abs(d);
    s1_new = s1;
    s2_new = [zeros(D,1); s2(1:end-D)];
else
    s1_new = s1;
    s2_new = s2;
end

% delay
delay_s = d*kgrid.dt;
fprintf('Delay est. : %.2f us\n', delay_s*1e6);

% display
fig1 = figure('Name', 'Test');
fig1.WindowState = "maximized";
tl1 = tiledlayout(fig1, 2, 1);

ax1  = nexttile(tl1);
plot(ax1, t_arraycut*1e6, s1, '-r', 'LineWidth', 1);
title(ax1, 'Before alignment', 'Interpreter', 'latex');
hold(ax1, 'on'); grid(ax1, 'on'); grid(ax1, 'minor');
plot(ax1, t_arraycut*1e6, s2, '-b', 'LineWidth', 1);
axis(ax1, 'tight');
fontsize(ax1, 14, "points");
set(gca,'TickLabelInterpreter','latex');

ax2  = nexttile(tl1);
plot(ax2, t_arraycut*1e6, s1_new, '-r', 'LineWidth', 1);
title(ax2, 'After alignment', 'Interpreter', 'latex');
hold(ax2, 'on'); grid(ax2, 'on'); grid(ax2, 'minor'); axis(ax2, 'tight');
plot(ax2, t_arraycut*1e6, s2_new, '-b', 'LineWidth', 1);
fontsize(ax2, 14, "points");
set(gca,'TickLabelInterpreter','latex');

xlabel(tl1, 'Time ($\mu$s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(tl1, 'Amplitude', 'Interpreter', 'latex', 'FontSize', 14);
hl = legend(sprintf('Tx: %d', tx), sprintf('Rx: %d', rx));
hl.Layout.Tile = 'East';






