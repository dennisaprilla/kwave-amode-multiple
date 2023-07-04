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


%%

transducer_idx = -1 * ((1:transducer_cfg.element.n)' - (1:transducer_cfg.element.n));
x   = transducer_idx * transducer_cfg.element.pitch;
x_t = x ./ c0;

matrix_idcs = zeros(transducer_cfg.numsamples, transducer_cfg.element.n, transducer_cfg.element.n);

for i=1:transducer_cfg.element.n
    matrix_t           = ( sqrt((x_t(i,:).^2)' + (transducer_cfg.t_array.^2)) )'; % [s]
    matrix_idcs(:,:,i) = round(matrix_t/transducer_cfg.period)+1; % [array indices]
end


%%

linear_idcs_correction = 0:transducer_cfg.numsamples:transducer_cfg.numsamples*(transducer_cfg.element.n-1);

PNR_signal = zeros(transducer_cfg.numsamples, transducer_cfg.element.n);

for i=1:transducer_cfg.element.n
    IMG   = FMC_signal(:,:,i);
    Tx_Rx = matrix_idcs(:,:,i);
    Tx_Rx(Tx_Rx>transducer_cfg.numsamples) = transducer_cfg.numsamples;
    Tx_Rx = bsxfun(@plus, Tx_Rx, linear_idcs_correction);
    TMP   = IMG(Tx_Rx);
    
    PNR_signal(:,i) = sum(TMP, 2);
end

%%
figure;
yyaxis left;
plot(transducer_cfg.t_array, FMC_signal(:,2,2), '-r');
hold on;
plot(transducer_cfg.t_array, FMC_signal(:,2,1), '-g');
plot(transducer_cfg.t_array, FMC_signal(:,2,3), '-b');
ylim([-2 2])

yyaxis right; grid on; grid minor;
plot(transducer_cfg.t_array, PNR_signal(:,2), '-m');
ylim([-2 2])


