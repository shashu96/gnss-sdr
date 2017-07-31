/*!
 * \file tracking_2nd_PLL_filter.cc
 * \brief Implementation of a 2nd order PLL filter for tracking carrier loop.
 * \author Javier Arribas, 2011. jarribas(at)cttc.es
 *
 * Class that implements 2 order PLL filter for tracking carrier loop. The algorithm
 * is described in:
 * K.Borre, D.M.Akos, N.Bertelsen, P.Rinder, and S.~H.~Jensen, A Software-Defined
 * GPS and Galileo Receiver. A Single-Frequency Approach,
 * Birkhauser, 2007, Applied and Numerical Harmonic Analysis.
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

#include "tracking_2nd_PLL_filter.h"
#include "GPS_L1_CA.h"
#include <cmath>


void Tracking_2nd_PLL_filter::calculate_lopp_coef(float* tau1,float* tau2, float lbw, float zeta, float k)
{
    // Solve natural frequency
    float Wn;
    Wn = lbw*8*zeta / (4*zeta*zeta + 1);
    // solve for t1 & t2
    *tau1 = k / (Wn * Wn);
    *tau2 = (2.0 * zeta) / Wn;
}



void Tracking_2nd_PLL_filter::set_PLL_BW(float pll_bw_hz)
{
    //Calculate filter coefficient values
    d_pllnoisebandwidth = pll_bw_hz;
    calculate_lopp_coef(&d_tau1_carr, &d_tau2_carr, d_pllnoisebandwidth, d_plldampingratio, 0.25); // Calculate filter coefficient values
}

void Tracking_2nd_PLL_filter::initialize()
{
    // carrier/Costas loop parameters
    d_old_carr_nco   = 0.0;
    d_old_carr_error = 0.0;
    long l=3;
	est_out = 0;
	//ele[3][3] = {{1.0/36.0,0.0,0.0},{0.0,1.0/4.0,0.0},{0.0,0.0,1.0}};


    /*
     * Initialization of predicted, estimated and error matrix
     */
    int i;
    int j;
    for(i = 0; i < 3; i++)
        {
            for(j = 0; j < l + 1; j++)
                {
	                pred[i][j] = 0;
		        }
	    }

    est = new float*[3];
    for(i = 0; i < 3; i++)
        {
    	    est[i] = new int[1000];
    	    for(j = 0 ; j < 1000 ; j++)
    	        {
	                est[i][j] = 0; //estimation
    	        }
		}

    for(i = 0; i < l; i++)
	{
	    j = 0;
	    error[i][j] = 0;
    }

    /*
     * Kalman filter transition matrix
     * F and H
     */
    stat_tran_mod[3][3] = {{1,1,0.5},{0,1,1},{0,0,1}}; //F = [1 1 1/2;0 1 1;0 0 1];
    obser_mod[1][3] = {1,0,0}; //H=[1 0 0]
    trans_obser_mod[3][1] = {{1},{0},{0}}; //H'
    trans_stat_tran_mod[3][3] = {{1,0,0},{1,1,0},{0.5,1,1}}; //F'
    eye[3][3] = {{1,0,0,},{0,1,0},{0,0,1}}; //I identity matrix

}
/*
 * Wrapping within a bound
 */
double wrapping_filter(double wrap, float range)
{
	//long a = sizeof(wrap);
	//long b = sizeof(wrap[1][1]);
	//long len = a/b; //length of the wrap signal

	//int i,j=1;
	//for(int i=1;i<=len;i++)
	//{
		if(wrap > range)
			while(wrap > range)
				wrap = wrap - 2*range;

		else if(wrap < -range)
			while(wrap < -range)
				wrap = wrap + 2*range;
	//}
	return wrap;

}

/*
 * PLL second order FIR filter
 * Req Input in [Hz/Ti]
 * The output is in [Hz/s].
 */
float Tracking_2nd_PLL_filter::get_carrier_nco(float PLL_discriminator)
{
    float carr_nco;
    carr_nco = d_old_carr_nco + (d_tau2_carr/d_tau1_carr)*(PLL_discriminator - d_old_carr_error) + (PLL_discriminator + d_old_carr_error) * (d_pdi_carr/(2*d_tau1_carr));
    //carr_nco = d_old_carr_nco + (d_tau2_carr/d_tau1_carr)*(PLL_discriminator - d_old_carr_error) + PLL_discriminator * (d_pdi_carr/d_tau1_carr);
    d_old_carr_nco   = carr_nco;
    d_old_carr_error = PLL_discriminator;
    return carr_nco;
}

float Tracking_2nd_PLL_filter::get_carrier_kf_nco(float KF_discriminator, long d_fs_in)
{
    #define CN0_ESTIMATION_SAMPLES 20
	long cn0_lin_hz;
	double phas_noise_var;
	double proc_cov_mat[3][3];
	float carr_nco;
    long d_ts_in_sec = 1/d_fs_in;
	cn0_lin_hz = pow(10,(CN0_ESTIMATION_SAMPLES/10));
	//Initialization
    x_new_old[3][1] = {0 , 50*d_ts_in_sec , 100*pow(d_ts_in_sec,2)}; //predicted state
    P_new_old[3][3] = {{1/12,0,0} , {0,1,0} , {0,0,1}}; //predicted error covariance

	//Phase noise variance
	phas_noise_var = (d_fs_in/(8*GPS_PI*GPS_PI*cn0_lin_hz))*(1+(d_fs_in/2*cn0_lin_hz));
	proc_cov_mat = cov_cal(ele); //process covariance matrix

	//Process begins here
	est_out = kf_impl_alg(KF_discriminator,phas_noise_var,proc_cov_mat,x_new_old,P_new_old);
	carr_nco = est_out;
	return carr_nco;
}

/*
 * Kalman Filter algorithm implementation
 */
float** kf_impl_alg(double error_signal, double R, double Q[3][3], double x_new_old[3][1], double P_new_old[][3])
{
    long len = sizeof(error_signal);
    //initiazlize(len);
    double kal_gain[3][1]; //column matrix
    double x_new_new[3][1]; //column matrix
    //double** est = 0;

    int k;
    int p=1;
    int i;
    int m;
    int n;
    int P_new_new[3][3];
    for(k = 1; k <= len; k++)
        {

	        //Measurement prediction
            //error[k][0] = signal[k][0] - pred[0][k]; //error = y_k - y_k-1

	        //wrapping
	        //check for passing of error whether it is array or single value
            //error[k][0] = wrapping_filter(error[k][0] , 0.5);


            kal_gain[][] ={1;1;1};
            /*
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

            */

            //x_new_new = x_new_old + K*error(k);
            for(i = 0; i < 3; i++)
                {
            	    //x_new_new[i][0] = x_new_old[i][0] + kal_gain[i][0]*error[k][0];
            	x_new_new[i][0] = x_new_old[i][0] + (kal_gain[i][0]*error_signal);
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
    float carr = (float) est;
	return carr;
}

double cov_cal(double Qd[3][3])
{
    double Q[3][3];
    int i,j;

    for(i=0 ; i<3 ; i++)
        {
    	    for(j=0 ; j<3 ; j++)
		        {
    	    	    Q = 1.0e-14*Qd[i][j];
		        }
        }
    return Q;
}

Tracking_2nd_PLL_filter::Tracking_2nd_PLL_filter (float pdi_carr)
{
    //--- PLL variables --------------------------------------------------------
    d_pdi_carr = pdi_carr;// Summation interval for carrier
    //d_plldampingratio = 0.65;
    d_plldampingratio = 0.7;
}


Tracking_2nd_PLL_filter::Tracking_2nd_PLL_filter ()
{
    //--- PLL variables --------------------------------------------------------
    d_pdi_carr = 0.001;// Summation interval for carrier
    d_plldampingratio = 0.7;
}




Tracking_2nd_PLL_filter::~Tracking_2nd_PLL_filter ()
{}

void Tracking_2nd_PLL_filter::set_pdi(float pdi_carr)
{
    d_pdi_carr = pdi_carr; // Summation interval for code
}
