/*!
 * \file gps_l1_ca_kf_tracking_cc.cc
 * \brief Implementation of Kalman filter tracking block
 * \author Carlos Aviles, 2010. carlos.avilesr(at)googlemail.com
 *         Javier Arribas, 2011. jarribas(at)cttc.es
 *         Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *         Carles Fernandez
 *         Shashanka Joisa
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2015  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * GNSS-SDR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNSS-SDR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNSS-SDR. If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "gps_l1_ca_kf_tracking_cc.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <gnuradio/io_signature.h>
#include <glog/logging.h>
#include <volk_gnsssdr/volk_gnsssdr.h>
#include "gps_sdr_signal_processing.h"
#include "tracking_transition.h"
#include "tracking_kf_discriminator.h"
#include "lock_detectors.h"
#include "GPS_L1_CA.h"
#include "control_message_factory.h"


/*!
 * \todo Include in definition header file
 */
#define CN0_ESTIMATION_SAMPLES 20
#define MINIMUM_VALID_CN0 25
#define MAXIMUM_LOCK_FAIL_COUNTER 50
#define CARRIER_LOCK_THRESHOLD 0.85
#define CN0_VALUE 35

using google::LogMessage;

gps_l1_ca_kf_tracking_cc_sptr
gps_l1_ca_kf_make_tracking_cc()
    long if_freq,
    long fs_in,
	unsigned int vector_length,
	bool dump,
	std::string dump_filename)
{
    return gps_l1_ca_kf_tracking_cc_sptr(new Gps_L1_Ca_Kf_Tracking_cc(if_freq,
    		fs_in, vector_length, dump, dump_filename))
}

Gps_L1_Ca_Kf_Tracking_cc::Gps_L1_Ca_Kf_Tracking_cc(
		long if_freq,
		long fs_in,
		unsigned int vector_length,
		bool dump,
		std::string dump_filename):
		gr::block("Gps_L1_Ca_Kf_Tracking_cc", gr::io_signature::make(1, 1, sizeof(gr_complex)),
		                gr::io_signature::make(1, 1, sizeof(Gnss_Synchro)))
{
	// Telemetry bit synchronization message port input
	this->message_port_register_in(pmt::mp("preamble_timestamp_s"));
	this->message_port_register_out(pmt::mp("events"));

	d_dump = dump;
	d_if_freq = if_freq;
	d_fs_in = fs_in;
	d_vector_length = vector_length;
	d_dump_filename = dump_filename;

	d_ts_in_sec = 1/d_fs_in;

	//Initialization
    x_new_old[3][1] = {0 , 50*d_ts_in_sec , 100*pow(d_ts_in_sec,2)}; //predicted state
    P_new_old[3][3] = {{1/12,0,0} , {0,1,0} , {0,0,1}}; //predicted error covariance


    // sample synchronization
    d_sample_counter = 0;
    d_acq_sample_stamp = 0;

	d_enable_tracking = false;
	d_pull_in = false;

    d_acquisition_gnss_synchro = 0;
    d_channel = 0;
    d_acq_code_phase_samples = 0.0;
    d_acq_carrier_doppler_hz = 0.0;
    d_carrier_doppler_hz = 0.0;
    d_acc_carrier_phase_rad = 0.0;
    d_code_phase_samples = 0.0;
    d_rem_code_phase_chips = 0.0;
    d_code_phase_step_chips = 0.0;
    d_carrier_phase_step_rad = 0.0;
}

void Gps_L1_Ca_Kf_Tracking_cc::start_tracking()
{
    /*
     *  correct the code phase according to the delay between acq and trk
     */

	//not sure about this -->start

    d_acq_code_phase_samples = d_acquisition_gnss_synchro->Acq_delay_samples;
    d_acq_carrier_doppler_hz = d_acquisition_gnss_synchro->Acq_doppler_hz;
    d_acq_sample_stamp = d_acquisition_gnss_synchro->Acq_samplestamp_samples;

    long int acq_trk_diff_samples;
    double acq_trk_diff_seconds;
    acq_trk_diff_samples = static_cast<long int>(d_sample_counter) - static_cast<long int>(d_acq_sample_stamp); //-d_vector_length;
    DLOG(INFO) << "Number of samples between Acquisition and Tracking =" << acq_trk_diff_samples;
    acq_trk_diff_seconds = static_cast<float>(acq_trk_diff_samples) / static_cast<float>(d_fs_in);
    // Doppler effect
    // Fd=(C/(C+Vr))*F
    double radial_velocity = (GPS_L1_FREQ_HZ + d_acq_carrier_doppler_hz) / GPS_L1_FREQ_HZ;
    // new chip and prn sequence periods based on acq Doppler
    double T_chip_mod_seconds;
    double T_prn_mod_seconds;
    double T_prn_mod_samples;
    d_code_freq_chips = radial_velocity * GPS_L1_CA_CODE_RATE_HZ;
    d_code_phase_step_chips = static_cast<double>(d_code_freq_chips) / static_cast<double>(d_fs_in);
    T_chip_mod_seconds = 1/d_code_freq_chips;
    T_prn_mod_seconds = T_chip_mod_seconds * GPS_L1_CA_CODE_LENGTH_CHIPS;
    T_prn_mod_samples = T_prn_mod_seconds * static_cast<double>(d_fs_in);

    d_current_prn_length_samples = round(T_prn_mod_samples);

    double T_prn_true_seconds = GPS_L1_CA_CODE_LENGTH_CHIPS / GPS_L1_CA_CODE_RATE_HZ;
    double T_prn_true_samples = T_prn_true_seconds * static_cast<double>(d_fs_in);
    double T_prn_diff_seconds = T_prn_true_seconds - T_prn_mod_seconds;
    double N_prn_diff = acq_trk_diff_seconds / T_prn_true_seconds;
    double corrected_acq_phase_samples, delay_correction_samples;
    corrected_acq_phase_samples = fmod((d_acq_code_phase_samples + T_prn_diff_seconds * N_prn_diff * static_cast<double>(d_fs_in)), T_prn_true_samples);
    if (corrected_acq_phase_samples < 0)
        {
            corrected_acq_phase_samples = T_prn_mod_samples + corrected_acq_phase_samples;
        }
    delay_correction_samples = d_acq_code_phase_samples - corrected_acq_phase_samples;

    d_acq_code_phase_samples = corrected_acq_phase_samples;

    d_carrier_doppler_hz = d_acq_carrier_doppler_hz;
    d_carrier_phase_step_rad = GPS_TWO_PI * d_carrier_doppler_hz / static_cast<double>(d_fs_in);

    // generate local reference ALWAYS starting at chip 1 (1 sample per chip)
    gps_l1_ca_code_gen_complex(d_ca_code, d_acquisition_gnss_synchro->PRN, 0);

    multicorrelator_cpu.set_local_code_and_taps(static_cast<int>(GPS_L1_CA_CODE_LENGTH_CHIPS), d_ca_code, d_local_code_shift_chips);

    d_carrier_lock_fail_counter = 0;
    d_rem_code_phase_samples = 0;
    d_rem_carr_phase_rad = 0.0;
    d_rem_code_phase_chips = 0.0;
    d_acc_carrier_phase_rad = 0.0;

    d_code_phase_samples = d_acq_code_phase_samples;

    std::string sys_ = &d_acquisition_gnss_synchro->System;
    sys = sys_.substr(0,1);

    //--> end

    // DEBUG OUTPUT
    std::cout << "Tracking start on channel " << d_channel << " for satellite " << Gnss_Satellite(systemName[sys], d_acquisition_gnss_synchro->PRN) << std::endl;
    LOG(INFO) << "Starting tracking of satellite " << Gnss_Satellite(systemName[sys], d_acquisition_gnss_synchro->PRN) << " on channel " << d_channel;

    // enable tracking
    d_pull_in = true;
    d_enable_tracking = true;

    LOG(INFO) << "PULL-IN Doppler [Hz]=" << d_carrier_doppler_hz
            << " Code Phase correction [samples]=" << delay_correction_samples
            << " PULL-IN Code Phase [samples]=" << d_acq_code_phase_samples;
}

int Gps_L1_Ca_Kf_Tracking_cc::general_work(int noutput_items __attribute__((unused)), gr_vector_int &ninput_items __attribute__((unused)),
        gr_vector_const_void_star &input_items, gr_vector_void_star &output_items)
{
	double wrap_sig;
	double sig_hz;
    double proc_cov_mat[3][3];
    double ele[3][3] = {{1/36,0,0},{0,1/4,0},{0,0,1}};
    double** est_out = 0;

	// Block input data and block output stream pointers
	const gr_complex* in = (gr_complex*) input_items[0]; //PRN start block alignment
	Gnss_Synchro **out = (Gnss_Synchro **) &output_items[0];

	// GNSS_SYNCHRO OBJECT to interchange data between tracking->telemetry_decoder
	Gnss_Synchro current_synchro_data = Gnss_Synchro();

	if (d_enable_tracking == true)
	    {
	        // Fill the acquisition data
	        current_synchro_data = *d_acquisition_gnss_synchro;
	        // Receiver signal alignment
	        if (d_pull_in == true)
	            {
                    int samples_offset;
                    double acq_trk_shif_correction_samples;
                    int acq_to_trk_delay_samples;
                    acq_to_trk_delay_samples = d_sample_counter - d_acq_sample_stamp;
                    acq_trk_shif_correction_samples = d_current_prn_length_samples - fmod(static_cast<float>(acq_to_trk_delay_samples), static_cast<float>(d_current_prn_length_samples));
                    samples_offset = round(d_acq_code_phase_samples + acq_trk_shif_correction_samples);
                    current_synchro_data.Tracking_timestamp_secs = (static_cast<double>(d_sample_counter) + static_cast<double>(d_rem_code_phase_samples)) / static_cast<double>(d_fs_in);
                    d_sample_counter = d_sample_counter + samples_offset; // count for the processed samples
                    d_pull_in = false;
                    // take into account the carrier cycles accumulated in the pull in signal alignment
                    d_acc_carrier_phase_rad -= d_carrier_phase_step_rad * samples_offset;
                    current_synchro_data.Carrier_phase_rads = d_acc_carrier_phase_rad;
                    current_synchro_data.Carrier_Doppler_hz = d_carrier_doppler_hz;
                    *out[0] = current_synchro_data;
                    consume_each(samples_offset); // shift input to perform alignment with local replica
                    return 1;

	            }

            // ########################### KALMAN FILTER #########################################
            //Phase input samples
	        wrap_sig = kf_two_quadrant_atan(in);
            sig_hz = wrap_sig/GPS_TWO_PI;

            //Phase noise variance
            cn0_lin_hz = pow(10,(CN0_VALUE/10));
            //d_ts_in_sec = 1/d_fs_in;
            phas_noise_var = (d_fs_in/(8*GPS_PI*GPS_PI*cn0_lin_hz))*(1+(d_fs_in/2*cn0_lin_hz));

            //proc_cov_mat = 1.0e-14*{{1/36,0,0},{0,1/4,0},{0,0,1}}; //process covariance matrix
            proc_cov_mat = cov_cal(ele); //process covariance matrix

            //Process begins here
            est_out = kf_impl_alg(sig_hz,phas_noise_var,proc_cov_mat,x_new_old,P_new_old);
	    }

    //assign the GNURadio block output data
	if(d_dump)
	    {
		    // MULTIPLEXED FILE RECORDING - Record results to file
		    try
		        {
		    	    // KF commands
		    	    d_dump_file.write(reinterpret_cast<char*>(&sig_hz), sizeof(double));
		    	    d_dump_file.write(reinterpret_cast<char*>(&phas_noise_var), sizeof(double));
		    	    d_dump_file.write(reinterpret_cast<char*>(&est_out), sizeof(double));
		        }
		    catch (const std::ifstream::failure &e)
		        {
		            LOG(WARNING) << "Exception writing trk dump file " << e.what();
		        }
	    }



	return 1;
}

void Gps_L1_Ca_Kf_Tracking_cc::set_channel(unsigned int channel)
{
    d_channel = channel;
    LOG(INFO) << "Tracking Channel set to " << d_channel;
    // ############# ENABLE DATA FILE LOG #################
    if (d_dump == true)
        {
            if (d_dump_file.is_open() == false)
                {
                    try
                    {
                        d_dump_filename.append(boost::lexical_cast<std::string>(d_channel));
                        d_dump_filename.append(".dat");
                        d_dump_file.exceptions (std::ifstream::failbit | std::ifstream::badbit);
                        d_dump_file.open(d_dump_filename.c_str(), std::ios::out | std::ios::binary);
                        LOG(INFO) << "Tracking dump enabled on channel " << d_channel << " Log file: " << d_dump_filename.c_str();
                    }
                    catch (const std::ifstream::failure &e)
                    {
                        LOG(WARNING) << "channel " << d_channel << " Exception opening trk dump file " << e.what();
                    }
                }
        }
}
