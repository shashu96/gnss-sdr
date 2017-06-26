/*!
 * \file tracking_kf_discriminator.h
 * \brief Interface of a library with tracking.
 * \authors <ul>
 *          <li>
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

#ifndef GNSS_SDR_TRACKING_KF_DISCRIMINATOR_H_
#define GNSS_SDR_TRACKING_KF_DISCRIMINATOR_H_

#include <gnuradio/gr_complex.h>



/*! \brief Two quadrant arctan discriminator
 *
 */

double kf_two_quadrant_atan(gr_complex prompt_s1);

/*! \brief wrapping the signal within a range
 *
 */
double wrapping_filter(double wrap[][10], float range);
