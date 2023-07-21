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

% NOTE: The difference between tfm_clipped (this script) and tfm_original
% is that this script has some modification on the tfm algorithm, that is,
% each of the transdcuer is represented by multiple points.

% NOTE2: tfm_clipped2 (this script) is 99% the same to tfm_clipped, the
% only difference is only the display. The purpose of this script is to
% look at what happen if we play with the number of points which represent
% a single transducer.

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
transducer_cfg.numsamples    = kgrid.Nt;

%% 2) Setup Transducer Coordinates
% Here, i did something different to normal TFM (the tfm_original.mat). In
% our case, our transducers is too big and too far apart, and the resulting
% image from TFM is not good (a lot of circular arc effects).

% We model the transducer (big and far-apart) as multiple points instead of
% a single point. For example, like below:
% [1111111--2222222--33333333], instead of: [---1--------2--------3----].
% This make sure that the "circular effect" from the picture reduced.

% Resolution of the point gap within one transducer
transducer_res = 2 * 1e-3;

% Make x-points with the length of the width of single transducer, then
% replicate it to the number of the transducers. For ex:
% tmp_x = [ 0 1 2 3 4 ]
%         [ 0 1 2 3 4 ]
%         [ 0 1 2 3 4 ]
tmp_x       = repmat( 0:transducer_res:(transducer_cfg.element.width-transducer_res), transducer_cfg.element.n, 1)';

% add an x-offset (pitch) to the points according to each of transducers.
% With pitch=6, our matrix from comment above will be:
% x = [  0 1 2 3 4,  7 8 9 10 11,  14 15 16 17 18 ]
tmp_offset  = 0:transducer_cfg.element.pitch:(transducer_cfg.element.pitch*(transducer_cfg.element.n-1));
x           = reshape(tmp_x+tmp_offset, 1, []);

% Each x-points belongs to a transducer, this variable stores the index of
% the transducer for each point. For example: 
% x_idcs = [ 1 1 1 1 1,  2 2 2 2 2,  3 3 3 3 3]
tmp_x_idcs  = repmat( 1:transducer_cfg.element.n, size(tmp_x,1), 1);
x_idcs      = tmp_x_idcs(:);

% Add another offset, so that x=0 is in the middle. For example:
% [-9 -8 -7 -6 -5,  -2 -1 0 1 2,  5 6 7 8 9]
x = x - median(x); % [m]

% Set the y-position (based from the simulation)
y = 5 * 1e-3;      % [m]
% Encapsulate them as transducers coordinate
transducer_cfg.coordinates = [x; ones(size(x))*y]';

%% Setup Grid coordinates
% Note: in kgrid (from k-wave), x_vec represents vertical axis, y-vec 
% represents horizontal axis. Here, i change them, i want x as horizontal 
% axis, and z as vertical axis.

% Resolution of the pixel (grid)
pixel_res       = 0.01 * 1e-3;

% Define focus area
minmax_x        = [-0.010, 0.010];
minmax_z        = [0.0020, -0.008]; % a-mode start from up to down

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
matrix_idcs = round(matrix_t/transducer_cfg.period); % [array indices]

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
% (ORIGINAL TFM) Each column of matrix_idcs corresponds to each column in 
% FMC_signalraw(:,i). If we want to index second column of FMC_signalraw(:,i) 
% linearly, we need to add the index accordingly, by the number of the row 
% (that is, M, n_sample). So, it would be like:
% [ 1 M+1 2M+1 ]
% [ 2 M+2 2M+2 ]
% [ 3 M+3 2M+3 ]
% ...
% [ M  2M   3M ]

% But here, since we have multiple points for one transducer, we will have
% more column from matrix_idcs compare to our FMC_signalraw(:,i). For
% example: if we have 3 transducers, and each transducer represented as 5
% points, it means we have 3x5-columns for matrix_idcs and 3-columns for 
% FMC_signalraw(:,i). Each of the 5 is actually 1 transducer, so that means
% the indexing will be duplicated like below:
% [ 1 1 1 1 1  M+1 .. M+1  2M+1 .. 2M+1 ]  
% [ 2 2 2 2 2  M+2 .. M+2  2M+2 .. 2M+2 ]  
% [ 3 3 3 3 3  M+3 .. M+3  2M+3 .. 2M+3 ]  
% ...
% [ M M M M M  2M  ..  2M   3M  ..  3M  ]

% calculate the index correction, M = n_sample, each column will add another M
linear_idcs_correction = 0:transducer_cfg.numsamples:transducer_cfg.numsamples*(transducer_cfg.element.n-1);

% Since we have more rows for matrix_idcs, we duplicates the index correction 
linear_idcs_correction = repmat(linear_idcs_correction, size(tmp_x,1), 1);
linear_idcs_correction = linear_idcs_correction(:)';

% Allocate a variable which will store our "image". This is linear image, a
% long 1D array, with the length of grid_coordinate (see section above). 
% If we have [M,N] grid, this "linear image" will be the length of MxN.
linear_image = zeros(size(matrix_idcs, 1), 1);

% let's allocate another variable which will store also our "image". But
% this one will store each group of transducers.
linear_images = zeros(size(matrix_idcs, 1), transducer_cfg.element.n);
img_idx       = 1;

% This is matlab-efficient implementation of TFM, please check this link:
% https://nl.mathworks.com/matlabcentral/fileexchange/56971-total-focusing-method
% Since we are modelling the transducer as multiple points, the script is 
% modified to be suitable for our case. Please check the original script to
% understand what is going on, i put more comments in that script.
for i=1:size(transducer_cfg.coordinates,1)
    
    % Get the "image"
    IMG         = FMC_signal(:,:,x_idcs(i));
    IMG_h       = [hilbert(IMG(:,1)), hilbert(IMG(:,2)), hilbert(IMG(:,3))];

    % Calculate delay
    Tx_Rx = bsxfun(@plus, matrix_idcs, matrix_idcs(:,i));

    % Clamp the delay (if the delay more than the ultrasound fast-time)
    Tx_Rx(Tx_Rx>transducer_cfg.numsamples) = transducer_cfg.numsamples;

    % Correct the linear indexing
    Tx_Rx = bsxfun(@plus, Tx_Rx, linear_idcs_correction);

    % Index the image to take the amplitude value from the "image"
    TMP   = IMG(Tx_Rx);

    % Add up all the amplitudes (and store it, and wait for the next loop
    % to add more amplitudes from the next Tx). Here, we should expect
    % un-synchronized amplitude, will cancel each other.
    linear_image = linear_image + sum(TMP, 2);

    % test
    linear_images(:,img_idx) = linear_images(:,img_idx) + sum(TMP, 2);
    if(mod(i, size(tmp_x, 1))==0) img_idx = img_idx+1; end
end

%% Display 1

% reshape the linear image to actual dimension of image
image = reshape(linear_image, [length(z_range),length(x_range)]);

% post process the image
normalize = @(x) x./max(x(:));
image_processed1 = 20*log10(normalize(abs(image)));
image_processed2 = normalize(abs(image));

% display
fig1 = figure();
ax = axes(fig1);

imagesc(ax, x_range*1000, -z_range*1000, image_processed2);
title(ax, 'Reconstructed image by discrete-TFM', 'Interpreter','latex');
xlabel(ax, 'X (mm)', 'Interpreter','latex');
ylabel(ax, 'Depth (mm)', 'Interpreter','latex');
hold(ax, 'on'); 
axis(ax, 'equal');
grid(ax, 'on');
grid(ax, 'minor');
plot(ax, transducer_cfg.coordinates(:,1)* 1000, -transducer_cfg.coordinates(:,2)* 1000, 'or');

c = colorbar;
c.Label.String = 'Normalized Amplitude';
c.Label.Rotation = 270;
c.Label.VerticalAlignment = "bottom";
c.Label.Interpreter = 'latex';
set(c,'TickLabelInterpreter','latex')

fontsize(ax, 14, "points");
set(gca,'TickLabelInterpreter','latex');

%% Display 2

images = reshape(linear_images, [length(z_range),length(x_range),transducer_cfg.element.n]);

% post process the image
normalize = @(x) x./max(x(:));

% display
fig2 = figure();
tl2  = tiledlayout(fig2, 1, transducer_cfg.element.n);

for i=1:transducer_cfg.element.n
    % process the current image
    image_processed1 = abs(images(:,:,i));
    image_processed2 = normalize(abs(images(:,:,i)));

    ax = nexttile(tl2, i);
    imagesc(ax, x_range*1000, -z_range*1000, image_processed2);
    xlabel(ax, 'X (mm)', 'Interpreter','latex');
    ylabel(ax, 'Depth (mm)', 'Interpreter','latex');
    hold(ax, 'on'); 
    axis(ax, 'equal');
    axis(ax, 'tight');
    grid(ax, 'on');
    grid(ax, 'minor');
    plot(ax, transducer_cfg.coordinates(:,1)* 1000, -transducer_cfg.coordinates(:,2)* 1000, 'or');
    
    fontsize(ax, 14, "points");
    set(gca,'TickLabelInterpreter','latex');
end


title(tl2, 'Reconstructed image by discrete-TFM', 'Interpreter','latex');

c = colorbar;
c.Label.String = 'Normalized Amplitude';
c.Label.Rotation = 270;
c.Label.VerticalAlignment = "bottom";
c.Label.Interpreter = 'latex';
c.Layout.Tile = 'east';
set(c,'TickLabelInterpreter','latex');