/*!
 * \file tracking_transition.h
 * \brief Interface of a library with tracking.
 * \authors <ul>
 *          <li> Carles Fernandez
 *          <li> Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *          </ul>
 *
 * Library with a set of code tracking and carrier tracking discriminators
 * that is used by the tracking algorithms.
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

#ifndef GNSS_SDR_TRACKING_TRANSITION_H_
#define GNSS_SDR_TRACKING_TRANSITION_H_

#include <gnuradio/gr_complex.h>
class Tracking_Transition
{
private:
	double error[][];
	double pred[][];
	double est[][];
	double stat_tran_mod[3][3];
	double obser_mod[3][3];
	double trans_obser_mod[3][3];

public:

    /*! \brief Kalman filter transition matrix
     *
     */
    void kf_transition_matrix();

    /*! \brief Initialization of predicted, estimated and error matrix
     *
     */
     void initiazation(long l);

     /*! \brief Kalman Filter algorithm implementation
      *
      */
     double kf_impl_alg(double signal, double R, double Q, double x_new_old, double P_new_old);
};
#endif
