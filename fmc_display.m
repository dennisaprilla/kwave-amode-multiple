clear; close all;

load('sim_results\supporting_simdata.mat');

theta     = 0:-2.5:-20;
n_element = 3; 
tx_idcs   = 1:n_element;

% all_data is a [M,N] cell matrix, M -> angle, N-> emmitting transducer
all_data = {};
for i=1:length(theta)
    for j=tx_idcs
        str_matname = sprintf('signal_%s_%d.mat', num2str(theta(i), '%+03.f'), j);
        str_matpath = fullfile(pwd, 'sim_results', str_matname);
        load(str_matpath);
        all_data{i, j} = transducers;
    end
end
clear transducers;

%%

select_idx   = 9;
select_theta = theta(select_idx);

% select_data is a [M,N] matrix, M-> Tx, N-> Rx
FMC  = cat(1, all_data{select_idx, tx_idcs});

ds     = USRaw_dvector(2);   % [mm]
ds_end = USRaw_dvector(end); % [mm]
dt     = USRaw_tvector(2);   % [µs]
dt_end = USRaw_tvector(end); % [µs]

s_resolution = 0.1; %[mm]
s_ticks_mm   = [0:s_resolution:ds_end, ds_end];
s_ticks_idx  = knnsearch(USRaw_dvector', s_ticks_mm', 'K', 1);

FMC_grouped = zeros(length(tx_idcs), length(tx_idcs), length(s_ticks_idx));
for i=tx_idcs % row, Tx
    for j=tx_idcs % col, Rx
        for k=1:length(s_ticks_idx)-1
            group = FMC(i,j).signal_env(s_ticks_idx(k):s_ticks_idx(k+1));
            FMC_grouped(i, j, k) = mean(group);
        end
    end
end

%%

[X, Y, Z] = meshgrid(tx_idcs, tx_idcs, 1:length(s_ticks_mm));
XYZ_FMC   = [X(:), Y(:), Z(:), FMC_grouped(:)];

offset     = 5;
offset_idx = (length(tx_idcs)^2 * offset) + 1;

scale      = 200;
color_var  = round( scale * XYZ_FMC(offset_idx:end,4) / max(XYZ_FMC(offset_idx:end,4))) + 1;

f1 = figure;
scatter3( XYZ_FMC(offset_idx:end,1), ...
               XYZ_FMC(offset_idx:end, 2),  ...
              -XYZ_FMC(offset_idx:end, 3) * s_resolution, ...
               color_var, ...
               XYZ_FMC(offset_idx:end,4), ...
              'filled', 'MarkerEdgeColor', 'black', 'MarkerFaceAlpha', .7, 'MarkerEdgeAlpha',.7);
xlabel('RX');
ylabel('TX');
zlabel('depth (mm)');
title(sprintf('FMC (Plane angle=%.1f°)', select_theta));
axis equal;
colorbar;
colormap('hot');
view(35,60);

%%

f2 = figure;
h_tabgroup = uitabgroup(f2);

for i=tx_idcs
    str_tabtitle = sprintf('Tx %d', i);
    tab          = uitab(h_tabgroup, 'Title', str_tabtitle);
    tab_axes     = axes('Parent', tab);
    hold(tab_axes, 'on');
    
    for j=tx_idcs
        sub_axes = subplot(length(tx_idcs), 1, j);
        hold(sub_axes, 'on');

        plot(sub_axes, USRaw_dvector, FMC(i,j).signal_corr);
        xlabel(sub_axes, sprintf('Depth (mm), 1-trip, v=%d m/s', c0));
        ylabel(sub_axes, 'Amplitude');
        grid(sub_axes, 'on'); 
        axis(sub_axes, 'tight'); 

        str_axtitle = sprintf('US element-%d', j);
        if(j==i)
            str_axtitle = strcat(str_axtitle, ' (Emits)');
        end
        title(sub_axes, str_axtitle);

        plot(sub_axes, USRaw_dvector, FMC(i,j).signal_env, '-g', 'LineWidth', 1);
        plot(sub_axes, FMC(i,j).signal_peak.loc_mm, FMC(i,j).signal_peak.amp, 'or', 'MarkerFaceColor', 'r');
        text(sub_axes, FMC(i,j).signal_peak.loc_mm, FMC(i,j).signal_peak.amp, ...
             strcat('$\,\leftarrow$', sprintf('%.2f', FMC(i,j).signal_peak.loc_mm)), 'Interpreter', 'latex');
        ylim([-2e5 2e5]);
    end
end





















