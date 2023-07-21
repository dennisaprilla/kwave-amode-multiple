clear; clc; close all;

% load simulation environment data
load('sim_results\supporting_simdata.mat');

% load simulation a-mode data
theta     = 0:-2.5:-15;
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

% [0, 2.5, 5, 7.5, 10, 12.5, 15] degree;
select_idx = 5; 

% FMC_singalall is a [n_Tx, n_Rx] matrix
FMC_singalall = cat(1, all_data{select_idx, tx_idcs});

% FMC_signalraw will be [n_sample, n_tx, n_rx];
FMC_signal = cat(1, FMC_singalall(:,:).signal_raw);
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


%% 2) TEST 1, WRONG EQUATION

%{
v   = 1540;     % [m/s]
ndx = 1 * 8e-3; % [m] 
dt  = 1 / 2e6;  % [s]

% fogoh requires lower bound, or else there is square root of negative number in the equation
t_lim      = (2*ndx)/v;
% https://nl.mathworks.com/matlabcentral/answers/334202-how-can-i-round-numbers-to-multiple-of-any-number
t_limround = ceil(t_lim/dt)*dt;
% specifies domain for fogoh function (i.e. lower-bounded t_vector)
t_vector   = t_limround:dt:((40e-6)-dt);

% the equation is on my notebook
fogoh   = @(x) 2/v * sqrt( ((v*x)/2).^2 - ndx^2 );
t_prime = fogoh(t_vector);

fig1 = figure;
ax1 = axes(fig1);
plot(ax1, t_vector, t_prime, '-o', 'LineWidth', 1);
title(ax1, sprintf('Time transformation from receiving to transmitting probe (ndx = %.4f)', ndx), 'Interpreter','latex');
xlabel(ax1, 't in adjacent receiving probe ($\mu$s)', 'Interpreter','latex');
ylabel(ax1, 't in transmiting probe ($\mu$s)', 'Interpreter','latex');
axis(ax1, 'equal');
grid(ax1, 'on');
grid(ax1, 'minor');
ylim(ax1, [0, max(t_prime)]);

fontsize(ax1, 14, "points");
set(gca,'TickLabelInterpreter','latex');
%}

%% 2) TEST 2: CORRECT EQUATION

v   = 1540;     % [m/s]
ndx = 1 * 8e-3; % [m] 
dt  = 1 / 2e6;  % [s]

% fogoh requires lower bound, or else there is square root of negative number in the equation
t_lim      = ndx/v;
% https://nl.mathworks.com/matlabcentral/answers/334202-how-can-i-round-numbers-to-multiple-of-any-number
t_limround = ceil(t_lim/dt)*dt;
% specifies domain for fogoh function (i.e. lower-bounded t_vector)
t_vector   = t_limround:dt:((40e-6)-dt);

% t_vector = 0:dt:((40e-6)-dt);
% t_vector(1) = eps;

% the inverse is confirmed by matlab, so i am confident with the equation.
% syms x;
% hogof(x) = ( (v*x)/2 + sqrt( ((v*x)/2).^2 + ndx^2) ) /v;
% fogoh = finverse(hogof);
hogof = @(x) ( (v*x)/2 + sqrt( ((v*x)/2).^2 + ndx^2) ) /v;
fogoh = @(x) ( ((v*x).^2 - ndx^2) ./ (v^2*x) );

t_prime = fogoh(t_vector);

fig1a = figure;
ax1a = axes(fig1a);
plot(ax1a, t_vector, t_prime, '-o', 'LineWidth', 1);
title(ax1a, sprintf('Time transformation from receiving to transmitting probe (ndx = %.4f)', ndx), 'Interpreter','latex');
xlabel(ax1a, 't in adjacent receiving probe ($\mu$s)', 'Interpreter','latex');
ylabel(ax1a, 't in transmiting probe ($\mu$s)', 'Interpreter','latex');
axis(ax1a, 'equal');
grid(ax1a, 'on');
grid(ax1a, 'minor');
ylim(ax1a, [0, max(t_prime)]);

fontsize(ax1a, 14, "points");
set(gca,'TickLabelInterpreter','latex');

%% 3) TEST

% transducer_idx = -1 * ((1:transducer_cfg.element.n)' - (1:transducer_cfg.element.n));
transducer_idx = abs((1:transducer_cfg.element.n)' - (1:transducer_cfg.element.n));

v     = 1540;
ndx   = transducer_idx * transducer_cfg.element.pitch;
dt    = transducer_cfg.period;

%{
% fogoh requires lower bound, or else there is square root of negative number in the equation
t_lim      = (2*ndx)/v;
% https://nl.mathworks.com/matlabcentral/answers/334202-how-can-i-round-numbers-to-multiple-of-any-number
t_limround = ceil(t_lim/dt)*dt;

% the equation is on my notebook
fogoh   = @(x, ndx) 2/v * sqrt( ((v*x)/2).^2 - ndx^2 );

% specifies domain for fogoh function (i.e. lower-bounded t_vector)
t_vectors = {};
t_primes  = {};
for i=1:size(t_limround,1) % TXs
    for j=1:size(t_limround,2) % RXs
        t_vector       = t_limround(i,j):dt:(transducer_cfg.t_array(end));
        t_vectors{i,j} = t_vector;
        t_primes{i,j}  = fogoh(t_vector, ndx(i,j));
    end
end

offset     = length(t_primes{2,2}) - length(t_primes{2,3});
offset_idx = (offset+1):length(t_primes{2,2});

fig2 = figure('Name', 'Test');
ax2  = axes(fig2);
plot(ax2, transducer_cfg.t_array, FMC_signal(:,2,2), '-r', 'LineWidth', 1);
title(ax2, 'Ultrasound Signal', 'Interpreter', 'latex');
xlabel(ax2, 'Time ($\mu$s)', 'Interpreter', 'latex' );
ylabel(ax2, 'Amplitude', 'Interpreter', 'latex');
hold(ax2, 'on'); grid(ax2, 'on'); grid(ax2, 'minor');
plot(ax2, transducer_cfg.t_array, FMC_signal(:,2,3), '-b', 'LineWidth', 1);
plot(ax2, t_primes{2,3}, FMC_signal(offset_idx,2,3), '-m', 'LineWidth', 1);
legend('Transmitting transducer (2)', 'Receiving transducer (3)', 'Receiving transducer (3), transformed', 'Interpreter', 'latex');
axis(ax2, 'tight');
ylim(ax2, [-5 5])

fontsize(ax2, 14, "points");
set(gca,'TickLabelInterpreter','latex');
%}

% fogoh requires lower bound, or else there is square root of negative number in the equation
t_lim      = ndx/v;
% https://nl.mathworks.com/matlabcentral/answers/334202-how-can-i-round-numbers-to-multiple-of-any-number
t_limround = ceil(t_lim/dt)*dt;

% the equation is on my notebook
fogoh = @(x, ndx) ( ((v*x).^2 - ndx^2) ./ (v^2*x) );

% specifies domain for fogoh function (i.e. lower-bounded t_vector)
t_vectors = {};
t_primes  = {};
for i=1:size(t_limround,1) % TXs
    for j=1:size(t_limround,2) % RXs
        t_vector       = t_limround(i,j):dt:(transducer_cfg.t_array(end));
        t_vectors{i,j} = t_vector;
        t_primes{i,j}  = fogoh(t_vector, ndx(i,j));
    end
end


offset     = length(t_primes{2,2}) - length(t_primes{2,3});
offset_idx = (offset+1):length(t_primes{2,2});

fig2 = figure('Name', 'Test');
ax2  = axes(fig2);
plot(ax2, transducer_cfg.t_array, FMC_signal(:,2,2), '-r', 'LineWidth', 1);
title(ax2, 'Ultrasound Signal', 'Interpreter', 'latex');
xlabel(ax2, 'Time ($\mu$s)', 'Interpreter', 'latex' );
ylabel(ax2, 'Amplitude', 'Interpreter', 'latex');
hold(ax2, 'on'); grid(ax2, 'on'); grid(ax2, 'minor');
plot(ax2, transducer_cfg.t_array, FMC_signal(:,2,3), '-b', 'LineWidth', 1);
plot(ax2, t_primes{2,3}, FMC_signal(offset_idx,2,3), '-m', 'LineWidth', 1);
legend('Transmitting transducer (2)', 'Receiving transducer (3)', 'Receiving transducer (3), transformed', 'Interpreter', 'latex');
axis(ax2, 'tight');
ylim(ax2, [-5 5])

fontsize(ax2, 14, "points");
set(gca,'TickLabelInterpreter','latex');











