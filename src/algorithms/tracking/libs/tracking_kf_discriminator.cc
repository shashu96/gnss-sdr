/*!
 * \file tracking_kf_discriminator.cc
 * \brief Implementation of a library used by the tracking algorithms.
 * \authors <ul>
 *          <li> Carles Fernandez
 *          <li> Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *          </ul>
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

#include "tracking_kf_discriminator.h"
#include <cmath>

//  All the outputs are in RADIANS
/*
 * two quadrant arctan
 */
double kf_two_quadrant_atan(gr_complex prompt_s1)
{
    if (prompt_s1.real() != 0.0)
        {
	    return atan(prompt_s1.imag() / prompt_s1.real());
	}
    else
	{
	    return 0;
	}
}

double wrapping_filter(double wrap[][10], float range)
{
    long a = sizeof(wrap);
    long b = sizeof(wrap[1][1]);
    long len = a/b; //length of the wrap signal

    int i,j=1;
    for(int i=1;i<=len;i++)
        {
	    if(wrap[i][j] > range)
	        while(wrap[i][j] > range)
		    wrap[i][j] =wrap[i][j] - 2*range;

	    else if(wrap[i][j] < -range)
	        while(wrap[i][j] < -range)
		    wrap[i][j] = wrap[i][j] + 2*range;
	}
	return wrap;

}


}


