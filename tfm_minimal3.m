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
select_idx = 6; 

% FMC_singalall is a [n_Tx, n_Rx] matrix
FMC_singalall = cat(1, all_data{select_idx, tx_idcs});

% FMC_signalraw will be [n_sample, n_tx, n_rx];
FMC_signal = cat(1, FMC_singalall(:,:).signal_env);
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


%% 2) TEST

v   = 1540;     % [m/s]
ndx = 1 * 8e-3; % [m] 
dt  = 1 / 2e6;  % [s]

dtheta       = 2.5;
theta_vector = 0:dtheta:15;      % [deg]

dx       = 1 * 1e-3;
x_vector = 0:dx:20*1e-3; % [[m]]

% TDiff     = @(x, theta) x/v * ( ( 1/cos(deg2rad(theta)) - 1 ) );
TDiff2    = @(x, theta) ( ( 1./cos(deg2rad(2*theta)) - 1 ) )' * (x./v);

tdiff    = TDiff2(x_vector, theta_vector);     % [s]
ddiff    = TDiff2(x_vector, theta_vector) * v; % [m]

[X, Y] = meshgrid(x_vector, theta_vector);

fig1 = figure('Name', 'Surface', 'Position', [100 100 1800 600]);
tl1  = tiledlayout(fig1, 1,3);

ax1  = nexttile(tl1, 1);
surf(ax1, X* 1e3, Y, ddiff* 1e3,'FaceAlpha',0.4, 'HandleVisibility','off');
title('Depth and Angle');
xlabel(ax1, 'Depth (mm)');
ylabel(ax1, 'Plane (deg)');
zlabel(ax1, 'Depth Difference (mm)');
hold(ax1, 'on');

highlight_depth = 10 * 1e-3;
highlight_theta = 10;
highlight_depthidx = round(highlight_depth/dx)+1;
highlight_thetaidx = round(highlight_theta/dtheta)+1;

plot3(ax1, X(:, highlight_depthidx)*1e3, theta_vector, ddiff(:,highlight_depthidx)*1e3, '-or', 'LineWidth', 3);
plot3(ax1, x_vector*1e3, Y(highlight_thetaidx,:), ddiff(highlight_thetaidx,:)*1e3, '-ob', 'LineWidth', 3);

hL = legend(sprintf('%d mm Depth', x_vector(highlight_depthidx)*1e3), sprintf('%d deg Angle', theta_vector(highlight_thetaidx)));
hL.Layout.Tile = 'East';


ax2  = nexttile(tl1, 2);
surf(ax2, X* 1e3, Y, ddiff* 1e3,'FaceAlpha',0.4, 'HandleVisibility','off');
title(sprintf('Depth Fixed (%d mm)', x_vector(highlight_depthidx)*1e3));
xlabel(ax2, 'Depth (mm)');
ylabel(ax2, 'Plane (deg)');
zlabel(ax2, 'Depth Difference (mm)');
hold(ax2, 'on');
plot3(ax2, X(:, highlight_depthidx)*1e3, theta_vector, ddiff(:,highlight_depthidx)*1e3, '-or', 'LineWidth', 3);
view(ax2, [90, 0]);

ax3  = nexttile(tl1, 3);
surf(ax3, X* 1e3, Y, ddiff* 1e3,'FaceAlpha',0.4,'HandleVisibility','off');
title(sprintf('Angle Fixed (%d deg)', theta_vector(highlight_thetaidx)));
xlabel(ax3, 'Depth (mm)');
ylabel(ax3, 'Plane (deg)');
zlabel(ax3, 'Depth Difference (mm)');
hold(ax3, 'on');
plot3(ax3, x_vector*1e3, Y(highlight_thetaidx,:), ddiff(highlight_thetaidx,:)*1e3, '-ob', 'LineWidth', 3);
view(ax3, [0, 0]);


%%

v   = 1540;     % [m/s]
ndx = 1 * 8e-3; % [m] 
dt  = 1 / 2e6;  % [s]

theta_vector = theta(select_idx);  % [deg]
% angles       = 5;
% theta_vector = angles(select_idx);

% h(t) = t * 0.5*(1+1/cos(deg2rad(refl_angle)));
% h-1(t) = t * 2/(1+cos(deg2rad(refl_angle)));
correction_factor = 2/(1+1/cos(deg2rad(2*theta_vector)));

fig2 = figure('Name', 'Test');
%fig2.WindowState = 'maximized';
ax2  = axes(fig2);
plot(ax2, transducer_cfg.t_array, FMC_signal(:,2,2), '-r', 'LineWidth', 1);
title(ax2, 'Ultrasound Signal', 'Interpreter', 'latex');
xlabel(ax2, 'Time ($\mu$s)', 'Interpreter', 'latex' );
ylabel(ax2, 'Amplitude', 'Interpreter', 'latex');
hold(ax2, 'on'); grid(ax2, 'on'); grid(ax2, 'minor');
plot(ax2, transducer_cfg.t_array, FMC_signal(:,2,3), '-b', 'LineWidth', 1);
plot(ax2, transducer_cfg.t_array*correction_factor, FMC_signal(:,2,2), '--r', 'LineWidth', 1);
plot(ax2, transducer_cfg.t_array*correction_factor, FMC_signal(:,2,3), '--b', 'LineWidth', 1);
axis(ax2, 'tight');
% ylim(ax2, [-5 5]);
ylim(ax2, [0 1*1e5]);