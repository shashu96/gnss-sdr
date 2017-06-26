/*!
 * \file tracking_transition.cc
 * \brief Implementation of a library used by the tracking algorithms.
 * \authors Carles Fernandez
 *          Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *          Shashanka Joisa
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
#include <cmath>

/*
 * Kalman filter transition matrix
 * F and H
 */
void kf_transition_matrix()
{
     stat_tran_mod[3][3] = {{1,1,0.5},{0,1,1},{0,0,1}}; //F = [1 1 1/2;0 1 1;0 0 1];
	 obser_mod[1][3] = {1,0,0}; //H=[1 0 0]
	 trans_obser_mod[3][1] = {{1},{0},{0}}; //H'
	 trans_stat_tran_mod[3][3] = {{1,0,0},{1,1,0},{0.5,1,1}}; //F'
	 eye[3][3]={{1,0,0,},{0,1,0},{0,0,1}}; //I identity matrix

}

/*
 * Initialization of predicted, estimated and error matrix
 */
void initiazation(long l)
{
	 int i;
	 int j;
	 for(i = 0; i < 3; i++)
	    {
		    for(j = 0; j < l + 1; j++)
		        {
			        pred[i][j] = 0;
		        }
	    }

	 for(i = 0; i < 3; i++)
	     {
		     for(j = 0; j < l; j++)
		         {
			         est[i][j] = 0; //estimation
		         }
	     }

	 for(i = 0; i < l; i++)
	     {
		     j = 0;
			 error[i][j] = 0;
         }

}

/*
 * Kalman Filter algorithm implementation
 */
double** kf_impl_alg(double signal, double R, double Q[3][3], double x_new_old[3][1], double P_new_old[][3])
{
    long len = sizeof(signal);
	initiazlize(len);
	double kal_gain[3][1]; //column matrix
	double x_new_new[3][1]; //column matrix
	double** est = 0;

	int k;
	int p=1;
	int i;
	int m;
	int n;
	int P_new_new[3][3];
	for(k = 1; k <= len; k++)
	    {
		    //Measurement prediction
			error[k][0] = signal[k][0] - pred[0][k]; //error = y_k - y_k-1

			//wrapping
			//check for passing of error whether it is array or single value
            error[k][0] = wrapping_filter(error[k][0] , 0.5);

            //Mearurement update
            //Kalman Gain calculation
            for(m = 0; m < 3; m++) //limit of m not yet known
                {
            	    for(n = 0; n < 3; n++) //limit of n not yet known
            	        {
            		        //kal_gain = P_new_old*trantrans_obser_mod[3][1]
            		        num[m][0] = num[m][0] + P_new_old[m][n]*trans_obser_mod[n][0]; //numerator of Kalman Gain
            		        den_1[0][m] = den_1[0][m] + obser_mod[0][n]*P_new_old[n][m]; //Denominator part1 of Kalman Gain

            	        }
            	        den = den + den_1[0][m]*trans_obser_mod[m][0]; //Denominator part2 of Kalman Gain
                 }
             den = den + R;

            for(m = 0; m < 3; m++)
                {
            	    kal_gain[m][0] = num[m][0]/den; //Kalman Gain = numerator/denominator
                }

            //x_new_new = x_new_old + K*error(k);
            for(i = 0; i < 3; i++)
                {
            	    x_new_new[i][0] = x_new_old[i][0] + kal_gain[i][0]*error[k][0];
                }

            //wrapping
            //check for passing of error whether it is array or single value
            x_new_new[0][0] = wrapping_filter(x_new_new[0][0],0.5);

            //Estimation of error covariance (P_new_new)
            for(m = 0; m < 3; m++)
                {
                    for(n = 0; n < 3; n++)
                        {
                            first[m][n] = kal_gain[m][0]*obser_mod[0][n]; //K*H
                        	second[m][n] = eye[m][n] + first[m][n];//I-(K*H)
                        }
                }

            //P_new_new = (eye(3) - K*H)*P_new_old;
            for(m = 0; m < 3; m++)
                {
                    for(n = 0; n < 3; n++)
                        {
                            for(i = 0; i < 3; i++)
                        	    {
                        	        P_new_new[m][n] += second[n][i] * P_new_old[i][n];
                        	    }

                        }
                }

            //Prediction
            //x_new_old = F*x_new_new;
            for(m = 0; m < 3; m++)
                {
                    for(n = 0; n < 3; n++)
                        {
                    	    x_new_old[m][0] += stat_tran_mod[m][n] * x_new_new[n][0];
                        }
                }

            //wrapping
            //check for passing of error whether it is array or single value
            x_new_old[0][0] = wrapping_filter(x_new_old[0][0],0.5);

            //predicted measurement
            //prediction(:,k+1) = H*x_new_old
            for(i = 0; i < 3; i++)
                {
                    sum += obser_mod[0][i] * x_new_old[i][0];
                }

            for(i = 0; i < 3; i++)
                {
                    pred[i][k+1] = sum;
                }

            //wrapping
            //check for passing of error whether it is array or single value
            pred[0][k+1] = wrapping_filter(pred[0][k+1],0.5);

            //Predicted error covariance
            //P_new_old = F*P_new_new*F.' + Q
            for(m = 0; m < 3; m++)
                {
                    for(n = 0; n < 3; n++)
            		    {
            		        for(i = 0; i < 3; i++)
            		            {
            		        	    P_new_old_fr[m][n] += stat_tran_mod[n][i] * P_new_new[i][n];
            		            }
            		    }
                 }

            for(m = 0; m < 3; m++)
                {
                    for(n = 0; n < 3; n++)
                        {
                            for(i = 0; i < 3; i++)
                                {
                        		    P_new_old[m][n] += P_new_old_fr[n][i] * trans_stat_tran_mod[i][n];
                        		}
                            P_new_old[m][n] += Q[m][n];
                        }
                }

            //Final estimation
            for(i = 0; i < 3; i++)
                {
            	    est[i][k] = x_new_new[i][0];
                }
		}
	return est;
}

