function signal_corr = my_corr(signal_main, signal_kernel)

    % correlate using matlab function
    [S_corr, ~]  = xcorr(signal_main, signal_kernel);
    S_corr       = S_corr';

    % start from halfway of barkercode entering the US signal
    sample_start = length(S_corr) - ...
                   size(signal_main,2) - ...
                   floor( 0.5 * length(signal_kernel) ) +1;
    sample_end   = length(S_corr) - ...
                   floor( 0.5 * length(signal_kernel) );

    % cut the signal
    signal_corr  = S_corr(sample_start:sample_end);
end

