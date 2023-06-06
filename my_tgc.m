function [tgc_x,tgc_y] = my_tgc(tgc_mastergain, tgc_dacslope, tgc_dacdelay, tgc_maxgain, f, t_vector)

    tgc_slopesample = tgc_dacslope * (f * 1e6); % [dB/sample]
    
    idx       = knnsearch(t_vector', tgc_dacdelay, "K", 1);
    tgc_x1    = t_vector(1:idx);
    tgc_x2    = t_vector(idx+1:end);
    tgc_y1    = zeros(1, length(tgc_x1));
    tgc_y2    = ( tgc_slopesample * (1:length(tgc_x2)) );
    
    tgc_x = [tgc_x1, tgc_x2];
    tgc_y = [tgc_y1, tgc_y2] + tgc_mastergain;
    tgc_y(tgc_y>tgc_maxgain) = tgc_maxgain;
end

