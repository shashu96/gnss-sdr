% /*!
%  * \file gps_l1_ca_dll_kf_plot_sample.m
%  * \brief Read GNSS-SDR Tracking dump binary file using the provided
%  function and plot some internal variables
%  * \author Javier Arribas, 2011. jarribas(at)cttc.es
%  *         Shashanka Joisa, 2017. shjoisa1996(at)gmail.com
%  * -------------------------------------------------------------------------
%  *
%  * Copyright (C) 2010-2016  (see AUTHORS file for a list of contributors)
%  *
%  * GNSS-SDR is a software defined Global Navigation
%  *          Satellite Systems receiver
%  *
%  * This file is part of GNSS-SDR.
%  *
%  * GNSS-SDR is free software: you can redistribute it and/or modify
%  * it under the terms of the GNU General Public License as published by
%  * the Free Software Foundation, either version 3 of the License, or
%  * at your option) any later version.
%  *
%  * GNSS-SDR is distributed in the hope that it will be useful,
%  * but WITHOUT ANY WARRANTY; without even the implied warranty of
%  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  * GNU General Public License for more details.
%  *
%  * You should have received a copy of the GNU General Public License
%  * along with GNSS-SDR. If not, see <http://www.gnu.org/licenses/>.
%  *
%  * -------------------------------------------------------------------------
%  */ 
close all;
clear all;
samplingFreq       = 64e6/16;     %4 MHz[Hz]
channels=4;
path='H:\working-gnss-sdr\gnss-sdr\install\';
clear PRN_absolute_sample_start;
for N=1:1:channels
    tracking_log_path=[path 'tracking_ch' num2str(N-1) '.dat'];
    GNSS_tracking(N)= gps_l1_ca_dll_pll_read_tracking_dump(tracking_log_path);   
end

% GNSS-SDR format conversion to MATLAB GPS receiver

for N=1:1:channels
        trackResults(N).status='T'; %fake track
        trackResults(N).codeFreq=GNSS_tracking(N).code_freq_hz.';
        trackResults(N).carrFreq=GNSS_tracking(N).carrier_doppler_hz.';
        trackResults(N).dllDiscr       = GNSS_tracking(N).code_error.';
        trackResults(N).dllDiscrFilt   = GNSS_tracking(N).code_nco.';
        trackResults(N).pllDiscr       = GNSS_tracking(N).carr_error.';
        trackResults(N).pllDiscrFilt   = GNSS_tracking(N).carr_nco.';

        trackResults(N).I_P=GNSS_tracking(N).prompt_I.';
        trackResults(N).Q_P=GNSS_tracking(N).prompt_Q.';

        trackResults(N).I_E= GNSS_tracking(N).E.';
        trackResults(N).I_L = GNSS_tracking(N).L.';
        trackResults(N).Q_E = zeros(1,length(GNSS_tracking(N).E));
        trackResults(N).Q_L =zeros(1,length(GNSS_tracking(N).E));
        trackResults(N).PRN=N; %fake PRN
        
        % Use original MATLAB tracking plot function
        settings.numberOfChannels=channels;
        settings.msToProcess=length(GNSS_tracking(N).E);
        plotTracking(N,trackResults,settings)
end

% for N=1:1:channels
%  figure;
%  plot([GNSS_tracking(N).E,GNSS_tracking(N).P,GNSS_tracking(N).L],'-*');
%  title(['Early, Prompt, and Late correlator absolute value output for channel ' num2str(N)']);
%  figure;
%  plot(GNSS_tracking(N).prompt_I,GNSS_tracking(N).prompt_Q,'+');
%  title(['Navigation constellation plot for channel ' num2str(N)]);
%  figure;
%  
%  plot(GNSS_tracking(N).prompt_Q,'r');
%  hold on;
%  plot(GNSS_tracking(N).prompt_I);
%  title(['Navigation symbols I(red) Q(blue) for channel ' num2str(N)]);
% end


