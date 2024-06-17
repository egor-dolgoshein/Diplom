#include <stdio.h>
#include <string>
#include <cmath>
#include <cstring>
#include <vector>
#include "parameters.h"

#define VICHET 4
// #define DEG_NUM 21

// const float PI = std::atan(1.0f) * 4.f;

struct Params_delta
{
	float delta;
	float delta_ht;
	float delta_ht_hx_2;
	float delta_ht_hy_2;
	float exp_delta_ht;
	float exp_delta_ht_ht;
};

float Integrate_Trapeze_2(int finish1, float h1, int finish2, float h2, float **function)
{
	float sum = 0.f;
	float sum1 = 0.f;
	float sum2 = 0.f;
	sum = 0.25f * (function[0][0] + function[finish1][0] + function[0][finish2] + function[finish1][finish2]);
	for (int i = 1; i < finish1; i++)
	{
		sum1 += function[i][0] + function[i][finish2];
	}
	sum1 = 0.5f * sum1;
	for (int j = 1; j < finish2; j++)
	{
		sum2 += function[0][j] + function[finish1][j];
	}
	sum2 = 0.5f * sum2;
	sum += sum1 + sum2;
	for (int i = 1; i < finish1; i++)
	{
		for (int j = 1; j < finish2; j++)
		{
			sum += function[i][j];
		}
	}
	return h1 * h2 * sum;
}

float K1(float x)
{
	return 1.f - x * x;
}
float V1(Params param, float x, float y, float alpha_cos_tr, float mst1, float mst2, float mst3)
{
	// return (param.sigma2 * x + mst1) * K1(x) + (-x) * pow(K1(x), 0.5f) * (cos(y) * (alpha_cos_tr - mst3) - sin(y) * mst2); // another H orientation
	// return (param.sigma2 * x + mst1) * K1(x) + (-x) * pow(K1(x), 0.5f) * (sin(y) * (alpha_cos_tr - mst3) - cos(y) * mst2); // original case
	// return (param.sigma2 * x + mst1 + alpha_cos_tr * cos(param.angle)) * K1(x) + (-x) * pow(K1(x), 0.5f) * (sin(y) * (alpha_cos_tr * sin(param.angle) - mst3) - cos(y) * mst2); // with angle xi
	return K1(x) * (param.sigma2 * x - mst1 + alpha_cos_tr * param.cos_xi) -
		   x * pow(K1(x), 0.5f) * (sin(y) * (alpha_cos_tr * param.sin_xi * param.sin_psi - mst3) + (alpha_cos_tr * param.sin_xi * param.cos_psi - mst2) * cos(y));
}
float K2(float x)
{
	return 1 / (K1(x));
}
float V2(Params params, float x, float y, float alpha_cos_tr, float mst2, float mst3)
{
	// return -pow(K2(x), 0.5f) * ((alpha_cos_tr - mst3) * sin(y) + cos(y) * mst2); // another H orientation
	// return pow(K2(x), 0.5f) * ((alpha_cos_tr - mst3) * cos(y) + sin(y) * mst2); //* K2(x); original case
	return pow(K2(x), 0.5f) * ((alpha_cos_tr * params.sin_xi * params.sin_psi - mst3) * cos(y) +
							   sin(y) * (mst2 - alpha_cos_tr * params.sin_xi * params.cos_psi)); // with angles xi and psi
}

float Ax(Params param, Params_delta param_delta, float x, float y, float tr, float mst1, float mst2, float mst3)
{
	float x__h = x - hx__2;
	return param_delta.exp_delta_ht_ht * (-K1(x__h) - hx__2 * V1(param, x__h, y, tr, mst1, mst2, mst3));
}
float Bx(Params param, Params_delta param_delta, float x, float y, float tr, float mst1, float mst2, float mst3)
{
	float x__h = x - hx__2;
	float x_h = x + hx__2;
	return param_delta.exp_delta_ht * (param_delta.delta_ht_hx_2 + param.ht * (K1(x_h) + K1(x__h) + hx__2 * (V1(param, x_h, y, tr, mst1, mst2, mst3) - V1(param, x__h, y, tr, mst1, mst2, mst3))));
}
float Cx(Params param, Params_delta param_delta, float x, float y, float tr, float mst1, float mst2, float mst3)
{
	float x_h = x + hx__2;
	return param_delta.exp_delta_ht_ht * (-K1(x_h) + hx__2 * V1(param, x_h, y, tr, mst1, mst2, mst3));
}

float Bx_0(Params param, Params_delta param_delta, float x, float y, float tr, float mst1, float mst2, float mst3)
{
	float x_h = x + hx__2;
	return param_delta.exp_delta_ht * (param_delta.delta_ht_hx_2 + param.ht * (K1(x_h) + hx__2 * V1(param, x_h, y, tr, mst1, mst2, mst3)));
}
float Bx_N(Params param, Params_delta param_delta, float x, float y, float tr, float mst1, float mst2, float mst3)
{
	float x__h = x - hx__2;
	return param_delta.exp_delta_ht * (param_delta.delta_ht_hx_2 + param.ht * (K1(x__h) + hx__2 * V1(param, x__h, y, tr, mst1, mst2, mst3)));
}

float F(float W_prev, float W_curr, float W_next)
{
	return W_prev + W_curr + W_next;
}
float Ay(Params param, Params_delta param_delta, float x, float y, float tr, float mst2, float mst3)
{
	return param_delta.exp_delta_ht_ht * (-K2(x) - hy__2 * V2(param, x, y - hy__2, tr, mst2, mst3));
}
float By(Params param, Params_delta param_delta, float x, float y, float tr, float mst2, float mst3)
{
	return param_delta.exp_delta_ht * (param_delta.delta_ht_hy_2 + param.ht * (2 * K2(x) + hy__2 * (V2(param, x, y + hy__2, tr, mst2, mst3) - V2(param, x, y - hy__2, tr, mst2, mst3))));
}
float Cy(Params param, Params_delta param_delta, float x, float y, float tr, float mst2, float mst3)
{
	return param_delta.exp_delta_ht_ht * (-K2(x) + hy__2 * V2(param, x, y + hy__2, tr, mst2, mst3));
}

float By_0(Params param, Params_delta param_delta, float x, float y, float tr, float mst2, float mst3)
{
	return param_delta.exp_delta_ht * (param_delta.delta_ht_hy_2 + param.ht * (2 * K2(x) + hy__2 * V2(param, x, y + hy__2, tr, mst2, mst3)));
}
float By_MM(Params param, Params_delta param_delta, float x, float y, float tr, float mst2, float mst3)
{
	return param_delta.exp_delta_ht * (param_delta.delta_ht_hy_2 + param.ht * (2 * K2(x) + hy__2 * V2(param, x, y - hy__2, tr, mst2, mst3)));
}
void Wprev_equal_Wcurr(int N, int MM, float **Wprev, float **Wcurr)
{
	for (int j = 0; j <= MM; j++)
	{
		for (int i = 0; i <= N; i++)
		{
			Wprev[i][j] = Wcurr[i][j];
		}
	}
}
void W_normirovka(int N, int MM, float **W)
{

	float inv_norm = (1.f / (Integrate_Trapeze_2(N, hx, MM, hy, W)));
	float checker = 0;
	for (int j = 0; j <= MM; j++)
	{
		for (int i = 0; i <= N; i++)
		{
			W[i][j] = inv_norm * W[i][j];
			checker += hx * hy * W[i][j];
		}
	}
	// printf("%f\n", checker);
}

void M_calc(int k, int N, int MM, const std::vector<float> &x, const std::vector<float> &y, float **W, std::vector<float> &M, std::vector<float> &mst_calc1, std::vector<float> &mst_calc2, std::vector<float> &mst_calc3, Params params)
{
	// std::cout << "k = " << k << '\n';
	float W_on_sqrt;
	float **W_on_sin_phi = new float *[N + 1];

	W_on_sin_phi[0] = new float[(N + 1) * (MM + 1)];
	for (int i = 1; i < N + 1; i++)
	{
		W_on_sin_phi[i] = W_on_sin_phi[0] + i * (MM + 1);
	}
	float **W_on_cos_phi = new float *[N + 1];
	W_on_cos_phi[0] = new float[(N + 1) * (MM + 1)];
	for (int i = 1; i < N + 1; i++)
	{
		W_on_cos_phi[i] = W_on_cos_phi[0] + i * (MM + 1);
	}
	float **W_on_x = new float *[N + 1];
	W_on_x[0] = new float[(N + 1) * (MM + 1)];
	for (int i = 1; i < N + 1; i++)
	{
		W_on_x[i] = W_on_x[0] + i * (MM + 1);
	}

	float **W_Integrating_function = new float *[N + 1];
	W_Integrating_function[0] = new float[(N + 1) * (MM + 1)];

	for (int i = 1; i < N + 1; ++i)
	{
		W_Integrating_function[i] = W_Integrating_function[0] + i * (MM + 1);
	}

	for (int j = 0; j <= MM; j++)
	{
		for (int i = 0; i <= N; i++)
		{
			W_on_sqrt = W[i][j] * pow(1 - x[i] * x[i], 0.5f);
			W_on_sin_phi[i][j] = W_on_sqrt * sin(y[j]);
			W_on_x[i][j] = W[i][j] * x[i];
			W_on_cos_phi[i][j] = W_on_sqrt * cos(y[j]);
			W_Integrating_function[i][j] = W_on_sin_phi[i][j] * params.sin_xi * params.sin_psi +
										   W_on_cos_phi[i][j] * params.sin_xi * params.cos_psi + W_on_x[i][j] * params.cos_xi; // with angle xi
																															   // std::cout << "W[" << i << "][" << j << "] * x = " << W_on_x[i][j] << '\n';
		}
	}
	M[k] = Integrate_Trapeze_2(N, hx, MM, hy, W_Integrating_function); // on sin phi
	// printf("%f\n", M[k]);
	mst_calc1[k] = Integrate_Trapeze_2(N, hx, MM, hy, W_on_x);
	mst_calc2[k] = Integrate_Trapeze_2(N, hx, MM, hy, W_on_cos_phi); // cos
	mst_calc3[k] = Integrate_Trapeze_2(N, hx, MM, hy, W_on_sin_phi); // M[k]

	delete[] W_on_sin_phi[0];
	delete[] W_on_sin_phi;
	delete[] W_on_cos_phi[0];
	delete[] W_on_cos_phi;
	delete[] W_on_x[0];
	delete[] W_on_x;
	delete[] W_Integrating_function[0];
	delete[] W_Integrating_function;
}

float W_start(float alp, float x, float y, float r)
{
	return 1.f + alp * sin(y) * pow(K1(x), 0.5f) / (1 + pow(r / 2, 2));
}
void ProgonkaX(int N, int j, const std::vector<float> &A, const std::vector<float> &B, const std::vector<float> &C, const std::vector<float> &F, float **result)
{
	std::vector<float> alpha;
	alpha.reserve(N + 1);
	std::vector<float> beta;
	beta.reserve(N + 1);
	std::vector<float> gamma;
	gamma.reserve(N + 1);
	std::vector<float> w;
	w.reserve(N + 1);
	std::vector<float> v;
	v.reserve(N + 1);
	float inv_denom;
	float W0;

	alpha[1] = 0.f;
	beta[1] = 0.f;
	gamma[1] = 1.f;
	for (int i = 1; i < N; i++)
	{
		inv_denom = 1.0f / (A[i] * alpha[i] + B[i]);
		alpha[i + 1] = -C[i] * inv_denom;
		beta[i + 1] = (F[i] - A[i] * beta[i]) * inv_denom;
		gamma[i + 1] = (-A[i] * gamma[i]) * inv_denom;
	}
	w[N] = 0.f;
	v[N] = 1.f;

	for (int i = N - 1; i >= 0; i--)
	{
		w[i] = alpha[i + 1] * w[i + 1] + beta[i + 1];
		v[i] = alpha[i + 1] * v[i + 1] + gamma[i + 1];
	}
	W0 = (F[0] - A[0] * w[N - 1] - C[0] * w[1]) / (B[0] + A[0] * v[N - 1] + C[0] * v[1]);

	for (int i = 0; i < N + 1; i++)
	{
		result[i][j] = W0 * v[i] + w[i];
	}
	alpha.clear();
	beta.clear();
	gamma.clear();
	v.clear();
	w.clear();
}
void ProgonkaY(int N, int i, const std::vector<float> &A, const std::vector<float> &B, const std::vector<float> &C, const std::vector<float> &F, float **result)
{
	std::vector<float> alpha;
	alpha.reserve(N + 1);
	std::vector<float> beta;
	beta.reserve(N + 1);
	std::vector<float> gamma;
	gamma.reserve(N + 1);
	std::vector<float> w;
	w.reserve(N + 1);
	std::vector<float> v;
	v.reserve(N + 1);
	float inv_denom;
	float W0;

	alpha[1] = 0.f;
	beta[1] = 0.f;
	gamma[1] = 1.f;
	for (int j = 1; j < N; j++)
	{
		inv_denom = 1.0f / (A[j] * alpha[j] + B[j]);
		alpha[j + 1] = -C[j] * inv_denom;
		beta[j + 1] = (F[j] - A[j] * beta[j]) * inv_denom;
		gamma[j + 1] = (-A[j] * gamma[j]) * inv_denom;
	}
	w[N] = 0.f;
	v[N] = 1.f;
	for (int j = N - 1; j >= 0; j--)
	{
		w[j] = alpha[j + 1] * w[j + 1] + beta[j + 1];
		v[j] = alpha[j + 1] * v[j + 1] + gamma[j + 1];
	}
	W0 = (F[0] - A[0] * w[N - 1] - C[0] * w[1]) / (B[0] + A[0] * v[N - 1] + C[0] * v[1]);
	for (int j = 0; j < N + 1; j++)
	{
		result[i][j] = W0 * v[j] + w[j];
	}
	alpha.clear();
	beta.clear();
	gamma.clear();
	v.clear();
	w.clear();
}

void find_last_W_and_M(Params params, Params_delta params_delta, int N_Tr, int N, int MM,
					   const std::vector<float> &x, const std::vector<float> &y, const std::vector<float> &t, const std::vector<float> &mst1, const std::vector<float> &mst2, const std::vector<float> &mst3,
					   float **Wprev, float r, float **Wcurr, std::vector<float> &M, std::vector<float> &mst_calc1, std::vector<float> &mst_calc2, std::vector<float> &mst_calc3, bool flag)
{
	// bool FLAG_PAR_CASE = (params.angle == 0);
	std::vector<float> A_values_X;
	A_values_X.reserve(N + 1);
	std::vector<float> B_values_X;
	B_values_X.reserve(N + 1);
	std::vector<float> C_values_X;
	C_values_X.reserve(N + 1);
	std::vector<float> F_values_X;
	F_values_X.reserve(N + 1);
	std::vector<float> A_values_Y;
	A_values_Y.reserve(MM + 1);
	std::vector<float> B_values_Y;
	B_values_Y.reserve(MM + 1);
	std::vector<float> C_values_Y;
	C_values_Y.reserve(MM + 1);
	std::vector<float> F_values_Y;
	F_values_Y.reserve(MM + 1);
	float alpha_cos_t_sdvig_r = 0.f;
	float a, b, c;

	for (int i = 0; i < N + 1; i++)
	{
		A_values_X[i] = B_values_X[i] = C_values_X[i] = F_values_X[i] = 0.f;
	}
	for (int j = 0; j < MM + 1; j++)
	{
		A_values_Y[j] = B_values_Y[j] = C_values_Y[j] = F_values_Y[j] = 0.f;
	}
	if (flag)
		M_calc(0, N, MM, x, y, Wprev, M, mst_calc1, mst_calc2, mst_calc3, params);

	for (int k = 1; k <= N_Tr; k++)
	{
		alpha_cos_t_sdvig_r = params.alp * cos((t[k] - 0.5f * params.ht) * r);
		for (int j = 0; j <= MM; j++)
		{
			int i = 0;
			a = Ay(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			b = By(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			c = Cy(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			A_values_X[i] = hy_2 * Ax(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			B_values_X[i] = hy_2 * Bx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			C_values_X[i] = hy_2 * Cx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			if (j == 0)
				F_values_X[i] = -hx_2 * F(0.f, (By_0(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]) - VICHET) * Wprev[i][j], c * Wprev[i][j + 1]);
			if (j == MM)
				F_values_X[i] = -hx_2 * F(a * Wprev[i][j - 1], (By_MM(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]) - VICHET) * Wprev[i][j], 0.f);
			if ((j != 0) && (j != MM))
				F_values_X[i] = -hx_2 * F(a * Wprev[i][j - 1], (b - VICHET) * Wprev[i][j], c * Wprev[i][j + 1]);

			for (i = 1; i < N; i++)
			{
				a = Ay(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
				b = By(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
				c = Cy(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
				A_values_X[i] = hy_2 * Ax(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
				B_values_X[i] = hy_2 * Bx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
				C_values_X[i] = hy_2 * Cx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
				if (j == 0)
				{
					F_values_X[i] = -hx_2 * F(0.f, (By_0(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]) - VICHET) * Wprev[i][j], c * Wprev[i][j + 1]);
					continue;
				}
				if (j == MM)
				{
					F_values_X[i] = -hx_2 * F(a * Wprev[i][j - 1], (By_MM(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]) - VICHET) * Wprev[i][j], 0.f);
					continue;
				}
				F_values_X[i] = -hx_2 * F(a * Wprev[i][j - 1], (b - VICHET) * Wprev[i][j], c * Wprev[i][j + 1]);
			}
			ProgonkaX(N, j, A_values_X, B_values_X, C_values_X, F_values_X, Wcurr);
		}

		Wprev_equal_Wcurr(N, MM, Wprev, Wcurr);
		int i = 0;
		for (int j = 0; j < MM; j++)
		{
			b = Bx_0(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			c = Cx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			A_values_Y[j] = hx_2 * Ay(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			B_values_Y[j] = hx_2 * By(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			C_values_Y[j] = hx_2 * Cy(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			F_values_Y[j] = -hy_2 * F(0.f, (b - VICHET) * Wprev[i][j], c * Wprev[i + 1][j]);
		}
		ProgonkaY(MM, i, A_values_Y, B_values_Y, C_values_Y, F_values_Y, Wcurr);

		for (i = 1; i < N; i++)
		{
			for (int j = 0; j < MM; j++)
			{
				a = Ax(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
				b = Bx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
				c = Cx(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
				A_values_Y[j] = hx_2 * Ay(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
				B_values_Y[j] = hx_2 * By(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
				C_values_Y[j] = hx_2 * Cy(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
				F_values_Y[j] = -hy_2 * F(a * Wprev[i - 1][j], (b - VICHET) * Wprev[i][j], c * Wprev[i + 1][j]);
			}
			ProgonkaY(MM, i, A_values_Y, B_values_Y, C_values_Y, F_values_Y, Wcurr);
		}

		i = N;
		for (int j = 0; j < MM; j++)
		{
			a = Ax(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			b = Bx_N(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst1[k], mst2[k], mst3[k]);
			A_values_Y[j] = hx_2 * Ay(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			B_values_Y[j] = hx_2 * By(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			C_values_Y[j] = hx_2 * Cy(params, params_delta, x[i], y[j], alpha_cos_t_sdvig_r, mst2[k], mst3[k]);
			F_values_Y[j] = -hy_2 * F(a * Wprev[i - 1][j], (b - VICHET) * Wprev[i][j], 0.f);
		}
		ProgonkaY(MM, i, A_values_Y, B_values_Y, C_values_Y, F_values_Y, Wcurr);

		W_normirovka(N, MM, Wcurr);
		Wprev_equal_Wcurr(N, MM, Wprev, Wcurr);

		if (flag)
			M_calc(k, N, MM, x, y, Wprev, M, mst_calc1, mst_calc2, mst_calc3, params);

		A_values_X.clear();
		B_values_X.clear();
		C_values_X.clear();
		F_values_X.clear();
		A_values_Y.clear();
		B_values_Y.clear();
		C_values_Y.clear();
		F_values_Y.clear();
	}
}

void Xi_calc(Params params, Params_delta params_delta1, Params_delta params_delta2)
{
	float XiL_3_PI_alp = (3.f * params.XiL) / (PI * params.alp);

	const int N = int(2.f / hx) - 1; // numeration from zero
	const int MM = int(2.f * PI / hy) - 1;

	std::vector<float> x;
	x.reserve(N + 1);
	std::vector<float> y;
	y.reserve(MM + 1);
	for (int i = 0; i <= N; i++)
	{
		x[i] = -1.0f + (0.5f + i) * hx;
	}

	for (int j = 0; j <= MM; j++)
	{
		y[j] = (0.5f + j) * hy;
	}
	int deg_num = int(10 * std::log10(params.tau_omegaN / params.tau_omega0)) + 1;
	// int deg_num_first_half = int(10 * std::log10(params.tau_omegaN / params.tau_omega0)) + 1;
	// int deg_num_second_half = deg_num_first_half;
	// int deg_num = deg_num_first_half + deg_num_second_half;
	printf("deg_num = %d\n", deg_num);
	float *d = new float[deg_num];
	float *r = new float[deg_num]; // r = 2*tau*omega
	float *r_2 = new float[deg_num];
	r[0] = 2 * params.tau_omega0;
	r_2[0] = 0.5f * r[0];
	d[0] = std::log10(r[0]);
	// printf("%f\t%f\t%f\n", d[0], r[0], r_2[0]);
	for (int i = 1; i < deg_num; i++)
	{
		d[i] = d[0] + 0.1f * i;
		r[i] = pow(10.f, d[i]);
		r_2[i] = r[i] * 0.5f;
		// printf("%f\t%f\t%f\n", d[i], r[i], r_2[i]);
	}

	/*	float *d = new float[deg_num];
		float *r = new float[deg_num];
		float *r_2 = new float[deg_num];
		r_2[0] = params.tau_omega0;
		r[0] = 2.f * r_2[0];
		d[0] = std::log10(r[0]);
		// printf("%f\t%f\t%f\n", d[0], r[0], r_2[0]);
		for (int i = 1; i < deg_num_first_half; i++)
		{
			d[i] = d[0] + (std::log10(38.9 * 2) - d[0]) / (deg_num_first_half - 1) * i;
			// d[i] = d[0] + 0.087f * i;
			r[i] = pow(10.f, d[i]);
			r_2[i] = r[i] * 0.5f;
			printf("%f\t%f\t%f\n", d[i], r[i], r_2[i]);
		}
		std::cout << "\n";
		// d[0] = std::log10(38.9 * 2);
		double a = (std::log10(params.tau_omegaN * 2) - std::log10(38.9 * 2)) / (deg_num - deg_num_first_half);
		for (int i = deg_num_first_half; i < deg_num; i++)
		{
			d[i] = std::log10(38.9 * 2) + a * (i - deg_num_first_half + 1);
			r[i] = pow(10.f, d[i]);
			r_2[i] = r[i] * 0.5f;
		}
	*/
	/*for (int i = deg_num - 1; i > deg_num_first_half - 2; --i)
	{
		r[i] -= 30.562973;
		r_2[i] -= 15.281487;
	}*/

	// for (int i = 0; i < deg_num; ++i)
	//	printf("%f\t%f\t%f\n", d[i], r[i], r_2[i]);

	float Tr_max = (2.0f * PI / r[0]);
	int L_max = int(Tr_max / params.ht);

	std::vector<float> t;
	t.reserve(L_max + 1);
	t[0] = 0.0f;
	for (int m = 1; m <= L_max; m++)
	{
		t[m] = t[m - 1] + params.ht;
	}

	char buf[100];
	/*int res = snprintf(buf, sizeof(buf),
					   "REXi_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, xi / PI, psi / PI, hx, params.ht);
	FILE *Re_Xi_file_interaction;
	Re_Xi_file_interaction = fopen(buf, "w");

	int res = snprintf(buf, sizeof(buf),
					   "REXi_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
					   params.alp, params.sigma, params.XiL, xi / PI, psi / PI, hx, params.ht);
	FILE *Re_Xi_file_without_interaction;
	Re_Xi_file_without_interaction = fopen(buf, "w");

		res = snprintf(buf, sizeof(buf),
					   "IMXi_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, xi / PI, psi / PI, hx, params.ht);
		FILE *Im_Xi_file_interaction;
		Im_Xi_file_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "IMXi_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, xi / PI, psi / PI, hx, params.ht);
	FILE *Im_Xi_file_without_interaction;
	Im_Xi_file_without_interaction = fopen(buf, "w");
	*/
	int res1 = snprintf(buf, sizeof(buf),
						"M_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
						params.alp, params.sigma, params.XiL, params.xi / PI, params.psi / PI, hx, params.ht);
	FILE *M_file_interaction;
	M_file_interaction = fopen(buf, "w");

	int res2 = snprintf(buf, sizeof(buf),
						"M_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
						params.alp, params.sigma, params.XiL, params.xi / PI, params.psi / PI, hx, params.ht);
	FILE *M_file_without_interaction;
	M_file_without_interaction = fopen(buf, "w");

	float **Wprev = new float *[N + 1];
	Wprev[0] = new float[(N + 1) * (MM + 1)];
	float **Wprev_ = new float *[N + 1];
	Wprev_[0] = new float[(N + 1) * (MM + 1)];
	float **Wcurr = new float *[N + 1];
	Wcurr[0] = new float[(N + 1) * (MM + 1)];
	float **Wcurr_ = new float *[N + 1];
	Wcurr_[0] = new float[(N + 1) * (MM + 1)];

	for (int i = 1; i < N + 1; i++)
	{
		Wprev[i] = Wprev[0] + i * (MM + 1);
		Wprev_[i] = Wprev_[0] + i * (MM + 1);
		Wcurr[i] = Wcurr[0] + i * (MM + 1);
		Wcurr_[i] = Wcurr_[0] + i * (MM + 1);
	}

	for (int s = 0; s < deg_num; s++)
	{
		float Tr = 2.f * PI / r[s];
		const int N_Tr = int(Tr / params.ht);
		std::vector<float> M;
		M.reserve(N_Tr + 1); 
		std::vector<float> mst_calc1;
		mst_calc1.reserve(N_Tr + 1);
		std::vector<float> mst_calc2;
		mst_calc2.reserve(N_Tr + 1);
		std::vector<float> mst_calc3;
		mst_calc3.reserve(N_Tr + 1);
		std::vector<float> mst1;
		mst1.reserve(N_Tr + 1); 
		std::vector<float> mst2;
		mst2.reserve(N_Tr + 1); 
		std::vector<float> mst3;
		mst3.reserve(N_Tr + 1); 
		std::vector<float> M_on_sin_rt;
		M_on_sin_rt.reserve(N_Tr + 1); 
		std::vector<float> M_on_cos_rt;
		M_on_cos_rt.reserve(N_Tr + 1); 
		
		for (int k = 0; k <= N_Tr; k++)
		{
			mst1[k] = 0.0f;
			mst_calc1[k] = 0.0f;
			mst2[k] = 0.0f;
			mst_calc2[k] = 0.0f;
			mst3[k] = 0.0f;
			mst_calc3[k] = 0.0f;
			M[k] = 0.f;
			M_on_sin_rt[k] = 0.f;
			M_on_cos_rt[k] = 0.f;
		}

		// find steady state value = W[t = t[N_Tr], x] starting from uniform value
		for (int j = 0; j <= MM; j++)
		{
			for (int i = 0; i <= N; i++)
			{
				Wprev[i][j] = 1. / 4. / PI * W_start(params.alp, x[i], y[j], r[s]);
			}
		}
		W_normirovka(N, MM, Wprev);

		find_last_W_and_M(params, params_delta1, N_Tr, N, MM, x, y, t, mst1, mst2, mst3, Wprev, r[s], Wcurr, M, mst_calc1, mst_calc2, mst_calc3, 0); // counting for tau*omega
		// now Wcurr contains last value of density - we will use it as initial value
		Wprev_equal_Wcurr(N, MM, Wprev, Wcurr);
		find_last_W_and_M(params, params_delta1, N_Tr, N, MM, x, y, t, mst1, mst2, mst3, Wprev, r[s], Wcurr, M, mst_calc1, mst_calc2, mst_calc3, 1);
		Wprev_equal_Wcurr(N, MM, Wprev, Wcurr);
		
		float rt = 0.f;
		for (int k = 0; k <= N_Tr; k++)
		{
			rt = r[s] * t[k];
			M_on_sin_rt[k] = M[k] * sin(rt);
			M_on_cos_rt[k] = M[k] * cos(rt);		
		}
		float im_chi = r[s] * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, params.ht, M_on_sin_rt);
		float re_chi = r[s] * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, params.ht, M_on_cos_rt);

		// fprintf(M_file_without_interaction, "%f\n", r_2[s]);
		for (int i = 0; i <= N_Tr; i++)
		{
			fprintf(M_file_without_interaction, "%f\n", M[i]);
		}
		// fprintf(Re_Xi_file_without_interaction, "%f %f\n", r_2[s], re_chi); // abscissa is tau*omega/2
		// fprintf(Im_Xi_file_without_interaction, "%f %f\n", r_2[s], im_chi);
		printf("%f\n", r_2[s]);

		//interactions on
		float Xil_na_4_Pi = params.XiL / (4.f * PI);
		float coeff_mst2_mst_3 = -params.XiL;
		for (int k = 0; k <= N_Tr; k++)
		{
			mst1[k] = Xil_na_4_Pi * mst_calc1[k];
			mst2[k] = coeff_mst2_mst_3 * mst_calc2[k];
			mst3[k] = coeff_mst2_mst_3* mst_calc3[k];
		}

		for (int j = 0; j <= MM; j++)
		{
			for (int i = 0; i <= N; i++)
			{
				Wprev_[i][j] = W_start(params.alp, x[i], y[j], r[s]);
			}
		}
		W_normirovka(N, MM, Wprev_);

		find_last_W_and_M(params, params_delta2, N_Tr, N, MM, x, y, t, mst1, mst2, mst3, Wprev_, r[s], Wcurr_, M, mst_calc1, mst_calc2, mst_calc3, 0);

		// now Wcurr contains last value of density - we will use it as initial value
		Wprev_equal_Wcurr(N, MM, Wprev_, Wcurr_);
		find_last_W_and_M(params, params_delta2, N_Tr, N, MM, x, y, t, mst1, mst2, mst3, Wprev_, r[s], Wcurr_, M, mst_calc1, mst_calc2, mst_calc3, 1);

		for (int k = 0; k <= N_Tr; k++)
		{
			rt = r[s] * t[k];
			M_on_sin_rt[k] = M[k] * sin(rt);
			M_on_cos_rt[k] = M[k] * cos(rt);
		}
		for (int i = 0; i <= N_Tr; i++)
		{
			fprintf(M_file_interaction, "%f\n", M[i]);
		}
		im_chi = r[s] * XiL_3_PI_alp *  Integrate_Trapeze(N_Tr, params.ht, M_on_sin_rt);
		// fprintf(Im_Xi_file_interaction, "%f %f\n", r_2[s], im_chi);

		re_chi = r[s] * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, params.ht, M_on_cos_rt);
		// fprintf(Re_Xi_file_interaction, "%f %f\n", r_2[s], re_chi);
		printf("%f\n", r_2[s]);
	}
	// fclose(Re_Xi_file_without_interaction);
	// fclose(Im_Xi_file_without_interaction);
	//  fclose(Re_Xi_file_interaction);
	//  fclose(Im_Xi_file_interaction);
	fclose(M_file_interaction);
	fclose(M_file_without_interaction);
	delete[] Wprev[0];
	delete[] Wprev;
	delete[] Wprev_[0];
	delete[] Wprev_;
	delete[] Wcurr[0];
	delete[] Wcurr;
	delete[] Wcurr_[0];
	delete[] Wcurr_;
	delete[] d;
	delete[] r;
	delete[] r_2;
	//printf("Calculation is finished. Files saved in current folder.\n");
}

int main(int argc, char *argv[])
{
	Params params;
	Params_delta params_delta1;
	Params_delta params_delta2;
	float xi;
	float psi;
	if (argc < 8)
		printf("Inditicate alpha, sigma, XiL\n");
	if (argc >= 8)
	{
		//std::string vzaim_str = argv[3];
		sscanf(argv[1], "%f", &params.alp);
		sscanf(argv[2], "%f", &params.sigma);
		sscanf(argv[3], "%f", &params.XiL);
		sscanf(argv[4], "%f", &params.ht);
		sscanf(argv[5], "%f", &params.tau_omega0);
		sscanf(argv[6], "%f", &params.tau_omegaN);
		sscanf(argv[7], "%f", &xi); // в градусах
		sscanf(argv[8], "%f", &psi);
	}
	xi = xi * PI / 180;
	psi = psi * PI / 180;
	params.xi = xi;
	params.psi = psi;
	params.sin_xi = sin(xi);
	params.cos_xi = cos(xi);
	params.sin_psi = sin(psi);
	params.cos_psi = cos(psi);
	params.sigma2 = params.sigma * 2;
	
	params_delta1.delta = 6.f * params.sigma + 2 * params.alp;
	params_delta1.delta_ht = 1.f + params_delta1.delta * params.ht;
	params_delta1.delta_ht_hx_2 = params_delta1.delta_ht * hx_2;
	params_delta1.delta_ht_hy_2 = params_delta1.delta_ht * hy_2;
	params_delta1.exp_delta_ht = exp(-params_delta1.delta * params.ht / 2.f);
	params_delta1.exp_delta_ht_ht = params.ht * exp(-params_delta1.delta * params.ht / 2.f);	

	params_delta2.delta = 6.f * params.sigma + 2 * params.alp + 2 * params.XiL;
	params_delta2.delta_ht = 1.f + params_delta2.delta * params.ht;
	params_delta2.delta_ht_hx_2 = params_delta2.delta_ht * hx_2;
	params_delta2.delta_ht_hy_2 = params_delta2.delta_ht * hy_2;
	params_delta2.exp_delta_ht = exp(-params_delta2.delta * params.ht / 2.f);
	params_delta2.exp_delta_ht_ht = params.ht * exp(-params_delta2.delta * params.ht / 2.f);
	Xi_calc(params, params_delta1, params_delta2);
	return 0;
}

