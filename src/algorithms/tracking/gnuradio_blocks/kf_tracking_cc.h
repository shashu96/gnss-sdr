/*!
 * \file gps_l1_ca_kf_tracking_cc.cc
 * \brief Interface of Kalman filter tracking block
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

#ifndef GNSS_SDR_GPS_L1_CA_KF_TRACKING_H
#define GNSS_SDR_GPS_L1_CA_KF_TRACKING_H

#include <fstream>
#include <map>
#include <string>
#include <gnuradio/block.h>
#include "gnss_synchro.h"

class Gps_L1_Ca_Kf_Tracking;

typedef boost_shared_ptr<Gps_L1_Ca_Kf_Tracking>
        gps_l1_ca_kf_tracking_cc_sptr;

gps_l1_ca_kf_tracking_cc_sptr
gps_l1_ca_kf_make_tracking_cc();

/*!
 *\brief This class implements Kalman Filter tracking block
 */
class Gps_L1_Ca_Kf_Tracking
{
public:
	~Gps_L1_Ca_Kf_Tracking(); //destructor
	void set_channel(unsigned int channel);
	void starttracking();
    int general_work (int noutput_items, gr_vector_int &ninput_items,
            gr_vector_const_void_star &input_items, gr_vector_void_star &output_items);

private:
	friend gps_l1_ca_kf_tracking_cc_sptr
	gps_l1_ca_kf_make_tracking_cc(long if_freq,
            long fs_in, unsigned
            int vector_length,
            bool dump,
            std::string dump_filename);

	Gps_L1_Ca_Kf_Tracking(
			long if_freq,
			long fs_in, unsigned
			int vector_length,
			bool dump,
		    std::string dump_filename); //consructor

	// tracking configuration vars
	unsigned int d_vector_length;
	bool d_dump;

	Gnss_Synchro* d_acquisition_gnss_synchro;
	unsigned int d_channel;

	long d_if_freq;
	long d_fs_in;
	long cn0_lin_hz;
	long d_ts_in_sec;

	// acquisition
	double d_acq_code_phase_samples;
	double d_acq_carrier_doppler_hz;

	double phas_noise_var;
	double x_new_old;
	double P_new_old;

	// tracking vars
	double d_code_freq_chips;
	double d_code_phase_step_chips;
	double d_carrier_doppler_hz;
	double d_carrier_phase_step_rad;
	double d_acc_carrier_phase_rad;
	double d_code_phase_samples;

	//PRN period in samples
	int d_current_prn_length_samples;

	//processing samples counters
	unsigned long int d_sample_counter;
	unsigned long int d_acq_sample_stamp;

	// control vars
	bool d_enable_tracking;
	bool d_pull_in;

	// file dump
	std::string d_dump_filename;
	std::ofstream d_dump_file;

	std::map<std::string, std::string> systemName;
	std::string sys;

};

#endif //GNSS_SDR_GPS_L1_CA_KF_TRACKING_H
