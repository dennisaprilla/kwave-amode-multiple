%% sIMULATION PARAMETERS ==================================================
clearvars; 

simconfig = ini2struct('simconf.ini');

medium_angle    = str2double(simconfig.medium_anglemin):str2double(simconfig.medium_anglestep):str2double(simconfig.medium_anglemax); % [deg]
medium_depth1   = str2double(simconfig.medium_depth1);    % [m]
medium_depth2   = str2double(simconfig.medium_depth2);    % [m]
transducers_pos = str2double(simconfig.transducers_pos);  % [m]
transducers_elm = str2double(simconfig.transducers_elm);
display_fig     = logical(str2double(simconfig.display_fig));
record_mov      = logical(str2double(simconfig.record_mov));

for medium_angle_current = medium_angle
    %% PREPARING THE VIRTUAL SETUP AND PARAMETERS OF SIMULATION ===========
    
    % 1) Create the computational grid ------------------------------------
    
    Nx = 768;        % number of grid points in the x (row) direction
    Ny = 1024;        % number of grid points in the y (column) direction
    dx = 0.025e-3;   % grid point spacing in the x direction [m]
    dy = dx;          % grid point spacing in the y direction [m]
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    
    % 2) Create medium ----------------------------------------------------
    
    % define the properties of the propagation medium1
    angle1    = medium_angle_current; % [deg]
    position1 = medium_depth1;         % [m]
    medium1_logic = makemedium_v2(angle1, [Nx, Ny], dx, position1);
    % define the properties of the propagation medium2
    angle2    = angle1;               % [deg]
    position2 = medium_depth2;         % [m]
    medium2_logic = makemedium_v2(angle2, [Nx, Ny], dx, position2);
    % add the two medium together
    medium_logic = medium1_logic+medium2_logic;
    
    c0 = 1540; % soft   [m/s]
    c1 = 1440; % fat    [m/s]
    c2 = 1588; % muscle [m/s]
    c3 = 3514; % bone   [m/s]
    d1 = 911;  % fat    [kg/m^3]
    d2 = 1090; % muscle [kg/m^3]
    d3 = 1908; % bone   [kg/m^3]
    a1 = 0.48; % fat    [dB/(MHz^y cm]
    a2 = 1.09; % muscle [dB/(MHz^y cm]
    a3 = 7.38; % bone   [dB/(MHz^y cm]
    
    speed1   = c1 * (medium_logic==2);
    speed2   = c2 * (medium_logic==1);
    speed3   = c3 * (medium_logic==0);
    density1 = d1 * (medium_logic==2);
    density2 = d2 * (medium_logic==1);
    density3 = d3 * (medium_logic==0);
    alpha1   = a1 * (medium_logic==2);
    alpha2   = a2 * (medium_logic==1);
    alpha3   = a3 * (medium_logic==0);
    
    medium.sound_speed = speed1 + speed2 + speed3;      
    medium.density     = density1 + density2 + density3;      
    medium.alpha_coeff = alpha1 + alpha2 + alpha3;
    medium.alpha_power = 1.1;
    
    % create the time array
    % kgrid.makeTime(medium.sound_speed);
    kgrid.makeTime(medium.sound_speed,  [], 18e-6);
    
    % 3) Create transducer and sensor -------------------------------------
    
    n_element      = transducers_elm;
    idx_element    = 1:n_element; 
    
    % variables for storing source and sensor
    % for source, i use binary mask (i think there is no other option)
    source_masks   = zeros([size(kgrid.x), n_element]);
    % for sensor, i use opposing corner. The resulting signals from simulation
    % are easier to handle, it returns a structured variable.
    sensor_masks   = zeros(4, n_element);
    
    % spesification of transducer
    line_center_m  = -[transducers_pos 0]; % [m]
    line_length_m  = 6e-3;                 % [m]
    line_kerf_m    = 2e-3;                 % [m]
    
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
    
    % setting up the trasnducers
    for idx_current=idx_element
        % determining the start grid position of the transducer
        idxpos_current  = idxpos_element(idx_current);
        start           = idxpos_current * (line_length_grid + line_kerf_grid) + line_halflength_grid;
        line_start_grid = line_center_grid + [0 start];
    
        % storing all the binary mask
        source_masks(:,:, idx_current) = makeLine(Nx, Ny, line_start_grid, deg2rad(0), line_length_grid);
    
        % storing all the opposing corner for sensor
        sensor_masks(:, idx_current) = [line_start_grid(1); line_start_grid(2)-line_length_grid; line_start_grid(1); line_start_grid(2)];
    end
    
    % assign the mask for the transducer
    % source.p_mask = source_masks(:,:,idx_tx);
    
    % create sensor
    sensor.mask = sensor_masks;
    
    % 4) Create waveform --------------------------------------------------
    sample_freq = 1/kgrid.dt; % [Hz]
    signal_freq = 7.5e6;      % [Hz]
    num_cycles  = 3;
    source.p    = toneBurst(sample_freq, signal_freq, num_cycles);
    
    for idx_tx = idx_element
        %% RUN THE ULTRASOUND SIMULATION ==================================
        
        % i put the initialization of source (transducer) here, so i can loop it
        % troughout element
        source.p_mask = source_masks(:,:,idx_tx);
        
        % parameter for recording the simulation
        str_moviename = sprintf('sim_%s_%s', num2str(medium_angle_current, '%+03.f'), num2str(idx_tx));
        str_moviepath = fullfile(pwd, 'sim_results', str_moviename);
        movieargs     = {'FrameRate', 12};

        % run the simulation
        input_args    = {'RecordMovie', record_mov, 'MovieName', str_moviepath, 'MovieArgs', movieargs};
        sensor_data   = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        
        % reorganize the data
        sensor_matrix = zeros(n_element, kgrid.Nt);
        for i=1:n_element
            sensor_matrix(i,:) = sum(squeeze(sensor_data(i).p),1);
        end
        
        %% POST-PROCESSING THE SIGNAL =====================================
        
        % get a suitable scaling factor for the time axis
        [~, t_scale, t_prefix] = scaleSI(kgrid.t_array(end));
        
        % preparing neccessary variable
        USRaw_matrix     = sensor_matrix;
        USRaw_signalfreq = 1/kgrid.dt;  % [Hz]
        USRaw_samplefreq = kgrid.dt;    % [s]
        USRaw_tvector    = kgrid.t_array * t_scale; % [mus]
        USRaw_dvector    = (kgrid.t_array * c0) * 1e3 / 2; % [mm]
        
        USBurst_data       = source.p;
        USBurst_signalfreq = signal_freq; % [Hz]
        USBurst_samplefreq = sample_freq; % [Hz]
        USBurst_dt         = 1/USBurst_samplefreq; % [s]
        USBurst_tvector    = (0:length(USBurst_data) - 1) * (USBurst_dt* t_scale);
        
        % tgc parameter
        tgc_mastergain  = 20;  % [dB]
        tgc_dacslope    = 1;   % [dB/mus]
        tgc_dacdelay    = 2;   % [mus]
        tgc_maxgain     = 40;  % [dB]
        [tgc_x,tgc_y]   = my_tgc(tgc_mastergain, tgc_dacslope, tgc_dacdelay, tgc_maxgain, USRaw_samplefreq, USRaw_tvector);
        
        % offset for detecting peaks (there is a near field turbulence)
        offset = 200;
        
        signal_peak = struct('amp', 0, 'loc_mm', [], 'loc_t', 0);
        transducer  = struct('signal_raw', [], 'signal_tgc', [], 'signal_corr', [], 'signal_env', [], 'signal_peak', struct(signal_peak) );
        transducers = transducer;
        
        for i=1:n_element
            % amplify with TGC
            signal_tgc   = USRaw_matrix(i,:) .* 10.^(tgc_y/10);
            % correlate
            signal_corr  = my_corr(signal_tgc, USBurst_data);
            % envelope
            signal_env   = envelope(signal_corr, 200, 'analytic');
        
            % find peaks
            [peaks, locs] = findpeaks( signal_env(offset+1:end), ...
                                      'MinPeakHeight', 0.5, ...
                                      'MinPeakProminence', 0.5, ...
                                      'SortStr', 'descend');
        
            locs_idx = locs(1)+offset;
            locs_t   = locs_idx*(kgrid.dt*1e6);
            locs_mm  = (c0 * 1e3) * (locs_t*1e-6) / 2;
        
            my_peaks = peaks(1);
        
            % store all the data
            transducer.signal_raw  = USRaw_matrix(i,:);
            transducer.signal_tgc  = signal_tgc;
            transducer.signal_corr = signal_corr;
            transducer.signal_env  = signal_env;
            transducer.signal_peak.amp    = my_peaks;
            transducer.signal_peak.loc_mm = locs_mm;
            transducer.signal_peak.loc_t  = locs_t;
        
            transducers(i) = transducer;
        end
        
        %% Visualizaton
        
        if(display_fig)
            f1 = figure("Name", "Ultrasound Signal");
            for i=1:n_element

                subplot(n_element, 1, i);
                plot(USRaw_dvector, transducers(i).signal_corr);
                xlabel(sprintf('Depth (mm), 1-trip, v=%d m/s', c0));
                ylabel('Amplitude');

                str_title = sprintf('US element-%d', i);
                if(i==idx_tx)
                    str_title = strcat(str_title, ' (Emits)');
                end
                title(str_title);

                grid on; axis tight; hold on;
                plot(USRaw_dvector, transducers(i).signal_env, '-g', 'LineWidth', 1);
                plot(transducers(i).signal_peak.loc_mm, transducers(i).signal_peak.amp, 'or', 'MarkerFaceColor', 'r');
                text(transducers(i).signal_peak.loc_mm, transducers(i).signal_peak.amp, ...
                     strcat('$\,\leftarrow$', sprintf('%.2f', transducers(i).signal_peak.loc_mm)), 'Interpreter', 'latex');
                ylim([-2e5 2e5]);

            end

            % save figure
            str_figurename = sprintf('figUS_%s_%s', num2str(medium_angle_current, '%+03.f'), num2str(idx_tx));
            str_figurepath = fullfile(pwd, 'sim_results', str_figurename);
            saveas(f1, str_figurepath);

            f2 = figure("Name", "Virtual Medium");
            subplot(2,2,1);
            imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, source.p_mask); colorbar;
            title('Tx');
            subplot(2,2,2);
            imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, sum(source_masks, 3)); colorbar;
            title('Rx');
            subplot(2,2,3);
            imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, medium.sound_speed); colorbar;
            c=colorbar;
            c.Label.String = 'm/s';
            c.Label.Interpreter = 'latex';
            title('Speed of Sound');
            subplot(2,2,4);
            imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, medium.alpha_coeff); colorbar;
            c=colorbar;
            c.Label.String = 'dB/(MHz cm)';
            c.Label.Interpreter = 'latex';
            title('Attenuation');

            % save figure
            str_figurename = sprintf('figVMed_%s_%s', num2str(medium_angle_current, '%+03.f'), num2str(idx_tx));
            str_figurepath = fullfile(pwd, 'sim_results', str_figurename);
            saveas(f2, str_figurepath);
        end

        str_matname = sprintf('signal_%s_%s', num2str(medium_angle_current, '%+03.f'), num2str(idx_tx));
        str_matpath = fullfile(pwd, 'sim_results', str_matname);
        save(str_matpath, 'transducers');

        close all;

        % break;

    % end element within transducer
    end

    % break;
    
% end angle
end






