close all;
clear all;

if ~exist('gps_l1_ca_dll_pll_read_tracking_dump.m','file')
    addpath('./libs')
end


samplingFreq = 2600000;     %[Hz]
channels = 2;
first_channel = 0;

path = 'H:\working gnss-sdr\gnss-sdr\install\';  %% CHANGE THIS PATH

for N=1:1:channels
    tracking_log_path = [path 'tracking_ch_0' num2str(N+first_channel-1) '.dat']; %% CHANGE epl_tracking_ch_ BY YOUR dump_filename
    GNSS_tracking(N)= gps_l1_ca_dll_pll_read_tracking_dump(tracking_log_path);
end

% GNSS-SDR format conversion to MATLAB GPS receiver

for N=1:1:channels
        trackResults(N).status = 'T'; %fake track
        trackResults(N).codeFreq       = GNSS_tracking(N).code_freq_hz.';
        trackResults(N).carrFreq       = GNSS_tracking(N).carrier_doppler_hz.';
        trackResults(N).dllDiscr       = GNSS_tracking(N).code_error.';
        trackResults(N).dllDiscrFilt   = GNSS_tracking(N).code_nco.';
        trackResults(N).pllDiscr       = GNSS_tracking(N).carr_error.';
        trackResults(N).pllDiscrFilt   = GNSS_tracking(N).carr_nco.';

        trackResults(N).I_P = GNSS_tracking(N).prompt_I.';
        trackResults(N).Q_P = GNSS_tracking(N).prompt_Q.';

        trackResults(N).I_E = GNSS_tracking(N).E.';
        trackResults(N).I_L = GNSS_tracking(N).L.';
        trackResults(N).Q_E = zeros(1,length(GNSS_tracking(N).E));
        trackResults(N).Q_L = zeros(1,length(GNSS_tracking(N).E));
        trackResults(N).PRN = ones(1,length(GNSS_tracking(N).E));
        
        % Use original MATLAB tracking plot function
        settings.numberOfChannels = channels;
        settings.msToProcess = length(GNSS_tracking(N).E);
        plotTracking(N,trackResults,settings)
end
