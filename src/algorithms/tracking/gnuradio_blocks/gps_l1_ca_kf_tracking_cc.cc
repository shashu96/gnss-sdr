/*!
 * \file gps_l1_ca_kf_tracking_cc.cc
 * \brief Implementation of Kalman filter tracking block
 * \author Carles Fernandez
 *         Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *         Carlos Aviles, 2010. carlos.avilesr(at)googlemail.com
 *         Javier Arribas, 2011. jarribas(at)cttc.es
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
{

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
    x_new_old = {0 , 50*d_ts_in_sec , 100*pow(d_ts_in_sec,2)}; //predicted state
    P_new_old = {{1/12,0,0} , {0,1,0} , {0,0,1}}; //predicted error covariance


	d_enable_tracking = false;

}

int Gps_L1_Ca_Kf_Tracking_cc::general_work(int noutput_items __attribute__((unused)), gr_vector_int &ninput_items __attribute__((unused)),
        gr_vector_const_void_star &input_items, gr_vector_void_star &output_items)
{
	double wrap_sig;
	double sig_hz;
    double proc_cov_mat;

	// Block input data and block output stream pointers
	//const gr_complex* in = (gr_complex*) input_items[0]; //PRN start block alignment
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

	            }

            // ########################### KALMAN FILTER #########################################
            //Phase input samples
	        wrap_sig = kf_two_quadrant_atan(input_signal);
            sig_hz = wrap_sig/GPS_TWO_PI;

            //Phase noise variance
            cn0_lin_hz = pow(10,(CN0_VALUE/10));
            //d_ts_in_sec = 1/d_fs_in;
            phas_noise_var = (d_fs_in/(8*GPS_PI*GPS_PI*cn0_lin_hz))*(1+(d_fs_in/2*cn0_lin_hz));

            proc_cov_mat = 1.0e-14*{{1/36,0,0},{0,1/4,0},{0,0,1}}; //process covariance matrix

            //Process begins here
            est_out = kf_impl_alg(sig_hz,phas_noise_var,proc_cov_mat,x_new_old,P_new_old);











}
