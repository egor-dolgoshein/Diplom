#include <stdio.h>
#include <string>
#include <cmath>
#include <cstring>
#include <iostream>
#include "parameters.h"

float Integrate_Trapeze(int finish, float h, float *function)
{
	float sum = 0.f;
	for (int i = 1; i < finish; i++)
	{
		sum += function[i];
	}
	return h * ((function[0] + function[finish]) * 0.5f + sum);
}

float K(float x)
{
	return 1.f - x * x;
}
float V(Params param, float x, float t, float MST, float r)
{
	// std::cout << (2.f * param.sigma * x + param.alp * cos(r * (t - ht__2)) + MST) * K(x) << '\n';
	return (2.f * param.sigma * x + param.alp * cos(r * (t - ht__2)) + MST) * K(x);
}

float A(Params param, float exp_delta_ht, float x, float t, float MST, float r)
{
	float x__h = x - hx__2;
	return exp_delta_ht * ht * (-K(x__h) - hx__2 * V(param, x__h, t, MST, r));
}
float C(Params param, float delta, float exp_delta_ht, float x, float t, float MST, float r)
{

	float x__h = x - hx__2;
	float x_h = x + hx__2;
	return exp_delta_ht * (hx_2 * (1.f + delta * ht) + ht * (K(x_h) + K(x__h) + hx__2 * V(param, x_h, t, MST, r) - hx__2 * V(param, x__h, t, MST, r)));
}
float B(Params param, float exp_delta_ht, float x, float t, float MST, float r)
{
	float x_h = x + hx__2;
	return exp_delta_ht * ht * (-K(x_h) + hx__2 * V(param, x_h, t, MST, r));
}
float F(float W_i)
{
	return hx_2 * W_i;
}

void Progonka(int N, float *A, float *B, float *C, float *F, float *result)
{
	float *alpha = new float[N + 1];
	float *beta = new float[N + 1];
	float inv_denom;

	alpha[1] = -C[0] / B[0];
	beta[1] = F[0] / B[0];

	for (int i = 1; i < N; i++)
	{
		inv_denom = 1.0f / (A[i] * alpha[i] + B[i]);
		alpha[i + 1] = -C[i] * inv_denom;
		beta[i + 1] = (F[i] - A[i] * beta[i]) * inv_denom;
	}

	result[N] = (F[N] - A[N] * beta[N]) / (B[N] + A[N] * alpha[N]);

	for (int i = N - 1; i >= 0; i--)
	{
		result[i] = alpha[i + 1] * result[i + 1] + beta[i + 1];
	}

	delete[] alpha;
	delete[] beta;
}

void find_last_W_and_M(Params params, float exp_delta_ht, float delta, int N_Tr, int N,
					   float *x, float *t, float *mst,
					   float *Wprev, float r, float *Wcurr, float *M, bool flag)
{
	float *A_values = new float[N + 1];
	float *B_values = new float[N + 1];
	float *C_values = new float[N + 1];
	float *F_values = new float[N + 1];
	float *W_on_x = new float[N + 1];
	for (int i = 0; i <= N; i++)
	{
		A_values[i] = B_values[i] = C_values[i] = F_values[i] = 0;
	}

	if (flag)
	{
		for (int i = 0; i <= N; i++)
		{
			W_on_x[i] = Wprev[i] * x[i];
		}
		M[0] = Integrate_Trapeze(N, hx, W_on_x);
	}

	for (int k = 1; k <= N_Tr; k++)
	{

		int i = 0;
		B_values[i] = C(params, delta, exp_delta_ht, x[i], t[k], mst[k], r);
		C_values[i] = B(params, exp_delta_ht, x[i], t[k], mst[k], r);
		F_values[i] = F(Wprev[i]);
		for (i = 1; i <= N - 1; i++)
		{
			A_values[i] = A(params, exp_delta_ht, x[i], t[k], mst[k], r);
			B_values[i] = C(params, delta, exp_delta_ht, x[i], t[k], mst[k], r);
			C_values[i] = B(params, exp_delta_ht, x[i], t[k], mst[k], r);
			F_values[i] = F(Wprev[i]);
		}

		i = N;
		A_values[i] = A(params, exp_delta_ht, x[i], t[k], mst[k], r);
		B_values[i] = C(params, delta, exp_delta_ht, x[i], t[k], mst[k], r);
		F_values[i] = F(Wprev[i]);
		Progonka(N, A_values, B_values, C_values, F_values, Wcurr);
		// normilizing koefficient now divided by 2*PI
		// originally it was without 2*PI
		float inv_norm = (1.f / (Integrate_Trapeze(N, hx, Wcurr)));
		// std::cout << "k = " << k << '\n';
		float checker = 0;
		for (int i = 0; i <= N; i++)
		{
			Wcurr[i] = inv_norm * Wcurr[i];
			checker += Wcurr[i];
			Wprev[i] = Wcurr[i];
			W_on_x[i] = Wcurr[i] * x[i];
		}
		// printf("%f\n", hx * checker);
		if (flag)
			M[k] = Integrate_Trapeze(N, hx, W_on_x);
		// printf("%f\n", M[k]);
	}
	delete[] A_values;
	delete[] B_values;
	delete[] C_values;
	delete[] F_values;
	delete[] W_on_x;
};

void Xi_calc(Params params)
{
	float XiL_3_PI_alp = (3.f * params.XiL) / (PI * params.alp);

	float delta = 3.f * params.sigma + params.alp;
	float exp_delta_ht = exp(-delta * ht);

	float delta_XiL = 3.f * params.sigma + params.alp + params.XiL; //?
	float exp_delta_ht_XiL = exp(-delta_XiL * ht);

	const int N = int(2.f / hx) - 1; //?
	float *x = new float[N + 1];

	for (int i = 0; i <= N; i++)
	{
		x[i] = -1.0f + (0.5f + i) * hx;
	}
	///??
	float d[deg_num];
	float r[deg_num]; // r = 2*tau*omega
	float r_2[deg_num];
	for (int i = 0; i < deg_num; i++)
	{
		d[i] = -2.0f + 0.1f * i;
		r[i] = 2 * pow(10.f, d[i]);
		r_2[i] = r[i] * 0.5f; // graphics abscissa is tau*omega
							  // r_2[i] = r[i];
							  // printf("%f\t%f\t%f\n", d[i], r[i], r_2[i]);
	}

	float Tr_max = (2.0f * PI / r[0]);
	int L_max = int(Tr_max / ht);
	float *t = new float[L_max + 1];

	t[0] = 0.0f;
	for (int m = 1; m <= L_max; m++)
	{
		t[m] = t[m - 1] + ht;
	}

	char buf[100];
	/*int res = snprintf(buf, sizeof(buf),
					   "REXi_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_hx%1.4f_ht%1.4f.txt",
					   params.alp, params.sigma, params.XiL, hx, ht);
	FILE *Re_Xi_file_interaction;
	Re_Xi_file_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "REXi_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, hx, ht);
	FILE *Re_Xi_file_without_interaction;
	Re_Xi_file_without_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "IMXi_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, hx, ht);
	FILE *Im_Xi_file_interaction;
	Im_Xi_file_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "IMXi_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, hx, ht);
	FILE *Im_Xi_file_without_interaction;
	Im_Xi_file_without_interaction = fopen(buf, "w");
	*/
	int res1 = snprintf(buf, sizeof(buf),
						"M_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_hx%1.4f_ht%1.4f.txt",
						params.alp, params.sigma, params.XiL, hx, params.ht);
	FILE *M_file_interaction;
	M_file_interaction = fopen(buf, "w");

	int res2 = snprintf(buf, sizeof(buf),
						"M_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_hx%1.4f_ht%1.4f.txt",
						params.alp, params.sigma, params.XiL, hx, params.ht);
	FILE *M_file_without_interaction;
	M_file_without_interaction = fopen(buf, "w");

	float *Wprev = new float[N + 1];
	float *Wprev_ = new float[N + 1];
	float *Wcurr = new float[N + 1];
	float *Wcurr_ = new float[N + 1];

	for (int s = 0; s < deg_num; s++)
	{
		float Tr = 2.f * PI / r[s];
		const int N_Tr = int(Tr / ht);
		float *M = new float[N_Tr + 1];
		float *mst = new float[N_Tr + 1];
		float *M_on_sin_rt = new float[N_Tr + 1];
		float *M_on_cos_rt = new float[N_Tr + 1];
		for (int i = 0; i <= N_Tr; i++)
		{
			mst[i] = 0.0f;
		}

		// find steady state value = W[t = t[N_Tr], x] starting from uniform value
		for (int i = 0; i <= N; i++)
		{
			Wprev[i] = 0.5f;
		}

		find_last_W_and_M(params, exp_delta_ht, delta, N_Tr, N, x, t, mst, Wprev, r[s], Wcurr, M, 0); // counting for 2*tau*omega

		// now Wcurr contains last value of density - we will use it as initial value
		for (int i = 0; i <= N; i++)
		{
			Wprev[i] = Wcurr[i];
		}

		find_last_W_and_M(params, exp_delta_ht, delta, N_Tr, N, x, t, mst, Wprev, r[s], Wcurr, M, 1);

		for (int i = 0; i <= N_Tr; i++)
		{
			M_on_sin_rt[i] = M[i] * sin(r[s] * t[i]);
			M_on_cos_rt[i] = M[i] * cos(r[s] * t[i]);
		}
		float im_chi = r[s] * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, ht, M_on_sin_rt);
		// fprintf(Im_Xi_file_without_interaction, "%f %f\n", r_2[s], im_chi);

		float re_chi = r[s] * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, ht, M_on_cos_rt);
		// fprintf(Re_Xi_file_without_interaction, "%f %f\n", r_2[s], re_chi); // abscissa is tau*omega
		for (int i = 0; i <= N_Tr; i++)
		{
			fprintf(M_file_without_interaction, "%f\n", M[i]);
		}
		//interactions on
		for (int k = 0; k <= N_Tr; k++) {
			mst[k] = params.XiL * M[k];
		}

		for (int i = 0; i <= N; i++)
		{
			Wprev_[i] = 0.5f;
		}

		find_last_W_and_M(params, exp_delta_ht_XiL, delta_XiL, N_Tr, N, x, t, mst, Wprev_, r[s], Wcurr_, M, 0);

		// now Wcurr contains last value of density - we will use it as initial value
		for (int i = 0; i <= N; i++)
		{
			Wprev_[i] = Wcurr_[i];
		}

		find_last_W_and_M(params, exp_delta_ht_XiL, delta_XiL, N_Tr, N, x, t, mst, Wprev_, r[s], Wcurr_, M, 1);

		for (int i = 0; i <= N_Tr; i++)
		{
			M_on_sin_rt[i] = M[i] * sin(r[s] * t[i]);
			M_on_cos_rt[i] = M[i] * cos(r[s] * t[i]);
		}

		im_chi = r[s] * XiL_3_PI_alp *  Integrate_Trapeze(N_Tr, ht, M_on_sin_rt);

		// fprintf(Im_Xi_file_interaction, "%f %f\n", r_2[s], im_chi);

		re_chi = r[s] * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, ht, M_on_cos_rt);
		// fprintf(Re_Xi_file_interaction, "%f %f\n", r_2[s], re_chi);
		for (int i = 0; i <= N_Tr; i++)
		{
			fprintf(M_file_interaction, "%f\n", M[i]);
		}

		delete[] M;
		delete[] mst;
		delete[] M_on_sin_rt;
		delete[] M_on_cos_rt;
	}

	printf("Calculation is finished. Files saved in current folder.\n");
	/*
		fclose(Re_Xi_file_interaction);
		fclose(Re_Xi_file_without_interaction);
		fclose(Im_Xi_file_interaction);
		fclose(Im_Xi_file_without_interaction);
	*/
	fclose(M_file_interaction);
	fclose(M_file_without_interaction);
	delete[] x;
	delete[] t;
	delete[] Wprev;
	delete[] Wprev_;
	delete[] Wcurr;
	delete[] Wcurr_;
}

int main(int argc, char *argv[])
{
    Params params;

    if (argc < 3)
        printf("Inditicate alpha, sigma, XiL\n");
	if (argc >= 3)
	{
		std::string vzaim_str = argv[3];

		sscanf(argv[1], "%f", &params.alp);
		sscanf(argv[2], "%f", &params.sigma);
		sscanf(argv[3], "%f", &params.XiL);
	}
	
	Xi_calc(params);

    return 0;
}
