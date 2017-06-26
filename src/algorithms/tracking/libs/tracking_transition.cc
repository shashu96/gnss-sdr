/*!
 * \file tracking_transition.cc
 * \brief Implementation of a library used by the tracking algorithms.
 * \authors <ul>
 *          <li>
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

#include "tracking_transition.h"
#include "tracking_kf_discriminator.h"
#include<cmath>

/*
 * Kalman filter transition matrix
 * F and H
 */
void kf_transition_matrix()
{
	stat_tran_mod[3][3] = {{1,1,0.5},{0,1,1},{0,0,1}};
	obser_mod[1][3] = {1,0,0}; //H=[1 0 0]
	trans_obser_mod[3][1] = {{1},{0},{0}};

}

/*
 * Initialization of predicted, estimated and error matrix
 */
void initiazation(long l)
{
	int i,j;
	for(i=1 ; i<=3 ; i++)
	{
		for(j=1 ; j<=l+1 ; j++)
		{
			pred[i][j] = 0;
		}
	}

	for(i=1 ; i<=3 ; i++)
	{
		for(j=1 ; j<=l ; j++)
		{
			est[i][j] = 0;
		}
	}

	for(i=1 ; i<=l ; i++)
		{
			    j=0;
				error[i][j] = 0;

		}

}

/*
 * Kalman Filter algorithm implementation
 */
double kf_impl_alg(double signal, double R, double Q, double x_new_old[3][3], double P_new_old[][3])
{
	long len = sizeof(signal);
	initiazlize(len);
	double kal_gain;

	int k,,p=1;
	for(k=1 ; k<=len ; k++)
		{
			//Measurement prediction
			error[k][1] = signal[k][1] - pred[1][k]; //error = y_k - y_k-1
            error[k][1] = wrapping_filter(error , 0.5);

            //Mearurement update
            //Kalman Gain calculation
            for(m=1;m<;m++)
            {
            	for(n=1;n<;n++)
            	{
            		//kal_gain = P_new_old*trantrans_obser_mod[3][1]
            		num[m][1] = num[m][1] + P_new_old[m][n]*trans_obser_mod[n][1]; //numerator of Kalman Gain
            		den_1[1][m] = den[1][m] + obser_mod[1][n]*P_new_old[n][m]; //Denominator part1 of Kalman Gain

            	}
            	den = den + den_1[1][m]*trans_obser_mod[m][1]; //Denominator part2 of Kalman Gain
            }
            den = den + R;

            for(m=0;m<;m++)
            {
            	kal_gain[m][1] = num[m][1]/den; //Kalman Gain = numerator/denominator
            }

            /*x_new_new =
            for(i=0;i<;i++)
            {

            }
            */




		}



}

