%%
% Before understanding this script, you should understand the general idea
% of how TFM work. Please check this paper: 
% Holmes, C., Drinkwater, B. and Wilcox, P., 2004. The post-processing of 
% ultrasonic array data using the total focusing method. Insight-Non-
% Destructive Testing and Condition Monitoring, 46(11), pp.677-680.
% https://doi.org/10.1784/insi.46.11.677.52285

% Or, you can go to one of this sites, to understand it moreintuitively:
% 1) https://www.olympus-ims.com/en/applications/using-the-total-focusing-method-to-improve-phased-array-ultrasonic-imaging/
% 2) https://www.eddyfi.com/en/technology/total-focusing-method-tfm
% 3) https://themysticwaves.com/total-focusing-method/

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
select_idx = 1; 

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
transducer_cfg.numsamples    = kgrid.Nt;

%% 2) Setup Transducer Coordinates
x = (0:transducer_cfg.element.n-1) * transducer_cfg.element.pitch;

if(mod(transducer_cfg.element.n, 2))
    mid_value = x( ceil(transducer_cfg.element.n/2) );
else
    mid_value = 0.5 * ( x(transducer_cfg.element.n/2) + x(transducer_cfg.element.n/2+1) );
end

x = x - mid_value;  % [m]
y = 5 * 1e-3;      % [m]
transducer_cfg.coordinates = [x; ones(size(x))*y]';

%% Setup Grid coordinates
% note: in kgrid, x_vec represents vertical axis, y-vec represents
% horizontal axis. here, i change them, i want x as horizontal axis, and z
% as vertical axis.

% Resolution of the pixel (grid)
pixel_res       = 0.1 * 1e-3;

% Define focus area
minmax_x        = [-0.010, 0.010];
minmax_z        = [0.005, -0.010]; % a-mode start from up to down

% Make a range from the specified res within the focus area
x_range         = (minmax_x(1):pixel_res:minmax_x(2))';
z_range         = (minmax_z(1)-pixel_res:-pixel_res:minmax_z(2))';

% Create the coordinate of the grid. If the focus area defined as an [M,N]
% enclosed area, grid_coordinate is a [MxN, 2] array (columns as x,y). The
% variable grid_coordinate will look like:
% grid_coordinate = [x1, y1]    -> 1
%                   [x2, y2]    -> 2
%                   [x3, y3]    -> 3
%                   ...         -> .
%                   [x12, y12]  -> 12
% each of the row of grid_coordinate represents a grid which indexed like:
% <example_grid> = [ 1 4 7 10 ]
%                  [ 2 5 8 11 ]
%                  [ 3 6 9 12 ]
[X, Z]          = meshgrid(x_range, z_range);
grid_coordinate = [X(:), Z(:)];

%% Create distance/time/indices matrix

% We need to calculate the distance of every grid_coordinate to every
% transducer coordinate.
matrix_s    = pdist2(grid_coordinate, transducer_cfg.coordinates); % [m]

% convert it to time
matrix_t    = matrix_s / c0; % [s]

% Convert it to A-mode indices. We can do this by dividing the time by the
% interval time of each sample, (dt, period, 1/f) of the transducer. 
% Dont forget to round it because the time can be "in-between" the index.
matrix_idcs = round(matrix_t/transducer_cfg.period)+1; % [array indices]

%% TFM

% Later, i will do linear indexing, a matrix but indexed as a linear array.
% In matlab, linear indexing a matrix is working like this
% [ 1  4  7 ]
% [ 2  5  8 ]
% [ 3  6  9 ]

% Our variable, matrix_idcs, is a matrix which the values will be used for
% indexing our MxN "image" (FMC_signalraw(:,i)), M -> n_sample, N -> n_rx. 
% This "image" is actually the A-mode signal, stacked together. 
% FMC_signalraw(:,i) means the "image" with i as the Tx index, and every 
% other, including i, as the Rx.
% 
% Each column of matrix_idcs corresponds to each column in FMC_signalraw(:,i). 
% If we want to index second column of FMC_signalraw(:,i) linearly, we need 
% to add the index accordingly, by the number of the row (that is, M, n_sample). 
% So, it would be like:
% [ 1 M+1 2M+1 ]
% [ 2 M+2 2M+2 ]
% [ 3 M+3 2M+3 ]
% ...
% [ M  2M   3M ]
linear_idcs_correction = 0:transducer_cfg.numsamples:transducer_cfg.numsamples*(transducer_cfg.element.n-1);

% Allocate a variable which will store our "image". This is linear image, a
% long 1D array, with the length of grid_coordinate (see section above). 
% If we have [M,N] grid, this "linear image" will be the length of MxN.
linear_image = zeros(size(matrix_idcs, 1), 1);

% This is matlab-efficient implementation of TFM, please check this link:
% https://nl.mathworks.com/matlabcentral/fileexchange/56971-total-focusing-method
for i=1:transducer_cfg.element.n
    
    % Each loop, represent each Tx (which the signal will caught by all Rx). 
    % TFM work by adding all the extracted amplitude for each combination 
    % of Tx and Rx. If we have 3-Tx and 3-Rx, there will be 9 combination. 
    % But we can do it sequentially for each Tx. So the resulting value 
    % will be stored (linear_image) and added as the loop progresses

    % Get the image (let's consider the first Tx)
    IMG     = FMC_signal(:,:,i);

    % Calculate delay. By adding each column matrix_idcs (let's consider
    % the first column, first Tx) to the whole matrix itself, you create a 
    % complete sound-travel that originated from first Tx. In the next
    % loop, you will create complete sound-travel that originated from 2nd
    % Tx. And so on until the last column, last Tx.
    Tx_Rx = bsxfun(@plus, matrix_idcs, matrix_idcs(:,i));

    % Clamp the delay (if the delay more than the ultrasound fast-time)
    Tx_Rx(Tx_Rx>transducer_cfg.numsamples) = transducer_cfg.numsamples;

    % Correct the linear indexing. This is the whole thing i described in
    % the beginning of this section.
    Tx_Rx = bsxfun(@plus, Tx_Rx, linear_idcs_correction);

    % Index the image to take the amplitude value from the "image"
    TMP   = IMG(Tx_Rx);

    % Add up all the amplitudes (and store it, and wait for the next loop
    % to add more amplitudes from the next Tx). Here, we should expect
    % un-synchronized amplitude, will cancel each other.
    linear_image = linear_image + sum(TMP, 2);
end

%% Display

image = reshape(linear_image, [length(z_range),length(x_range)]);

normalize = @(x) x./max(x(:));
image_processed1 = 20*log10(normalize(abs(image)));
image_processed2 = normalize(abs(image));

fig = figure();
ax = axes(fig);

imagesc(ax, x_range*1000, -z_range*1000, image_processed1);
title(ax, 'Reconstructed image by discrete-TFM');
xlabel(ax, 'X (mm)', 'Interpreter','latex');
ylabel(ax, 'Depth (mm)', 'Interpreter','latex');
hold(ax, 'on'); 
axis(ax, 'equal');
grid(ax, 'on');
grid(ax, 'minor');
plot(ax, transducer_cfg.coordinates(:,1)* 1000, -transducer_cfg.coordinates(:,2)* 1000, 'or');

c = colorbar;
c.Label.String = 'Power (dB)';
c.Label.Rotation = 270;
c.Label.VerticalAlignment = "bottom";
c.Label.Interpreter = 'latex';
set(c,'TickLabelInterpreter','latex')

fontsize(ax, 14, "points");
set(gca,'TickLabelInterpreter','latex');
















