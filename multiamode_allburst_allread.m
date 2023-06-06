%% =========================================================================
% SIMULATION
% =========================================================================

clearvars; close all;

% 1) create the computational grid ----------------------------------------
Nx = 512;      % number of grid points in the x (row) direction
Ny = Nx;       % number of grid points in the y (column) direction
dx = 0.05e-3;  % grid point spacing in the x direction [m]
dy = dx;       % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% 2) Create medium --------------------------------------------------------

% define the properties of the propagation medium
% medium.sound_speed = 1500;  % [m/s]
% medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
% medium.alpha_power = 1.25;

angle = 10;
medium_logic = makemedium_v2(angle, [Nx, Ny], dx, -5e-3);

c1 = 1560; % [m/s]
c2 = 3500; % [m/s]
d1 = 1049; % [kg/m^3]
d2 = 1908; % [kg/m^3]
a1 = 0.54; % [dB/(MHz^y cm]
a2 = 6.90; % [dB/(MHz^y cm]

speed1   = c1 * medium_logic;       
speed2   = c2 * (medium_logic==0);
density1 = d1 * medium_logic;      
density2 = d2 * (medium_logic==0); 
alpha1   = a1 * medium_logic;
alpha2   = a2 * (medium_logic==0); 

medium.sound_speed = speed1 + speed2;      
medium.density     = density1 + density2;      
medium.alpha_coeff = alpha1 + alpha2;
medium.alpha_power = 1.25;

% create the time array
kgrid.makeTime(medium.sound_speed,  [], 0.5e-6);

% 3) Create transducer and sensor -----------------------------------------

n_element      = 3;
idx_element    = 1:n_element; 
source_masks   = zeros([size(kgrid.x), n_element]);

% spesification of transducer
line_center_m  = [-10 0]*1e-3; % [m]
line_length_m  = 6e-3;         % [m]
line_kerf_m    = 2e-3;         % [m]

% specification in term of grids
line_center_grid     = [(Nx/2), Ny/2] + round(line_center_m/dx); % [grid]
line_length_grid     = round(line_length_m/dx);
line_halflength_grid = round((line_length_m/2)/dx);
line_kerf_grid       = round(line_kerf_m/dx);

% make indexes, for ex, 3 = [-1 0 1], 4 = [-1.5 -0.5 0.5 1.5]
if(mod(n_element,2))
    idxpos_element = idx_element - (ceil(n_element/2));
else
    idxpos_element = idx_element - (n_element/2+0.5);
end

for idx_current=idx_element
    % determining the start grid position of the transducer
    idxpos_current  = idxpos_element(idx_current);
    start           = idxpos_current * (line_length_grid + line_kerf_grid) + line_halflength_grid;
    line_start_grid = line_center_grid + [0 start];

    % storing all the mask
    source_masks(:,:, idx_current) = makeLine(Nx, Ny, line_start_grid, 0, line_length_grid);
end

% assign the mask for the transducer
source.p_mask = sum(source_masks, 3);

% create sensor
sensor.mask = source.p_mask;

% 4) Create waveform ------------------------------------------------------
sample_freq = 1/kgrid.dt; % [Hz]
signal_freq = 7.5e6;      % [Hz]
num_cycles  = 3;
source.p    = toneBurst(sample_freq, signal_freq, num_cycles, 'Plot', true);

% 5) Run the simulation ---------------------------------------------------
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor,'PlotLayout', true, 'PlotPML', false);








