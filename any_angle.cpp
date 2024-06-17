#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include "parameters.h"

std::string const PATH = "D:\\Documents\\Egor\\study\\Diplom\\"; // your working folder path

void calc(Params param, float r, std::vector<float> &t, float &im_chi, float &im_chi_int, float &re_chi, float &re_chi_int, int i,
		  std::ifstream &M_par_int, std::ifstream &M_par_w_int, std::ifstream &M_perp_int, std::ifstream &M_perp_w_int)
{

	float XiL_3_PI_alp = (3.f * param.XiL) / (PI * param.alp);
	float Tr = 2.f * PI / r;
	const int N_Tr = int(Tr / param.ht);
	// const int N = int(2.f / hx) - 1;
	// const int chunk = 150;
	// for (int l = 0; l < int(N_Tr / chunk) + 1; ++l)
	{
		// std::vector<float> M_ineraction(chunk + 1, 0);
		// std::vector<float> M_without_ineraction(chunk + 1, 0);
		float M_ineraction[N_Tr + 1];
		float M_without_ineraction[N_Tr + 1];
		float s = 0;

		for (int j = 0; j <= N_Tr; ++j)
		{
			M_perp_int >> s;
			float number = s; // std::stof(s);
			M_ineraction[j] = number * param.sin_xi * param.sin_xi;
		}
		s = 0;
		for (int j = 0; j <= N_Tr; ++j)
		{
			M_par_int >> s;
			float number = s; // std::stof(s);
			M_ineraction[j] += number * param.cos_xi * param.cos_xi;
		}
		s = 0;
		for (int j = 0; j <= N_Tr; ++j)
		{
			M_perp_w_int >> s;
			float number = s; // std::stof(s);
			M_without_ineraction[j] = number * param.sin_xi * param.sin_xi;
		}

		s = 0;
		for (int j = 0; j <= N_Tr; ++j)
		{
			M_par_w_int >> s;
			float number = s; // std::stof(s);
			M_without_ineraction[j] += number * param.cos_xi * param.cos_xi;
		}

		float M_on_sin_rt[N_Tr + 1];
		float M_on_cos_rt[N_Tr + 1];
		// std::vector<float> M_on_sin_rt(N_Tr + 1, 0);
		//	M_on_sin_rt_interaction.reserve(N_Tr + 1);
		// std::vector<float> M_on_cos_rt(N_Tr + 1, 0);
		//	M_on_cos_rt_interaction.reserve(N_Tr + 1);

		float rt = 0.f;
		for (int k = 0; k <= N_Tr; k++)
		{
			rt = r * t[k];
			M_on_sin_rt[k] = M_ineraction[k] * sin(rt);
			M_on_cos_rt[k] = M_ineraction[k] * cos(rt);
		}
		im_chi_int = r * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, param.ht, M_on_sin_rt);
		re_chi_int = r * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, param.ht, M_on_cos_rt);

		for (int k = 0; k <= N_Tr; k++)
		{
			rt = r * t[k];
			M_on_sin_rt[k] = M_without_ineraction[k] * sin(rt);
			M_on_cos_rt[k] = M_without_ineraction[k] * cos(rt);
		}
		im_chi = r * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, param.ht, M_on_sin_rt);
		re_chi = r * XiL_3_PI_alp * Integrate_Trapeze(N_Tr, param.ht, M_on_cos_rt);
	}
}

void init_and_set(Params params)
{
	char buf[100];
	int res = snprintf(buf, sizeof(buf),
					   "REXi_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
					   params.alp, params.sigma, params.XiL, params.xi / PI, params.psi / PI, hx, params.ht);
	FILE *Re_Xi_file_interaction;
	Re_Xi_file_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "REXi_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, params.xi / PI, params.psi / PI, hx, params.ht);
	FILE *Re_Xi_file_without_interaction;
	Re_Xi_file_without_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "IMXi_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, params.xi / PI, params.psi / PI, hx, params.ht);
	FILE *Im_Xi_file_interaction;
	Im_Xi_file_interaction = fopen(buf, "w");

	res = snprintf(buf, sizeof(buf),
				   "IMXi_without_interaction_alpha%1.1f_sigma%1.1f_chil%1.1f_xi%1.3fPI_psi%1.3fPI_hx%1.4f_ht%1.4f.txt",
				   params.alp, params.sigma, params.XiL, params.xi / PI, params.psi / PI, hx, params.ht);
	FILE *Im_Xi_file_without_interaction;
	Im_Xi_file_without_interaction = fopen(buf, "w");

	// res = snprintf(buf, sizeof(buf), "Relaxation_times.txt");
	// FILE *Relaxation_times;
	// Relaxation_times = fopen(buf, "w");

	float *d = new float[deg_num];
	float *r = new float[deg_num];
	float *r_2 = new float[deg_num];
	r[0] = 2 * params.tau_omega0;
	r_2[0] = 0.5f * r[0];
	d[0] = std::log10(r[0]);
	for (int i = 1; i < deg_num; i++)
	{
		d[i] = d[0] + 0.1f * i;
		r[i] = pow(10.f, d[i]);
		r_2[i] = r[i] * 0.5f;
		// printf("%f\t%f\t%f\n", d[i], r[i], r_2[i]);
	}

	float Tr_max = (2.0f * PI / r[0]);
	int L_max = int(Tr_max / params.ht);

	std::vector<float> t;
	t.reserve(L_max + 1);
	t[0] = 0.0f;
	for (int m = 1; m <= L_max; m++)
	{
		t[m] = t[m - 1] + params.ht;
	}
	std::ifstream M_perp_int;
	M_perp_int.open(PATH + "M_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.500PI_hx0.0500_ht0.0100.txt");
	std::ifstream M_perp_w_int;
	M_perp_w_int.open(PATH + "M_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.500PI_hx0.0500_ht0.0100.txt");
	std::ifstream M_par_int;
	M_par_int.open(PATH + "M_interaction_alpha0.0_sigma3.0_chil1.0_hx0.0500_ht0.0000.txt");
	std::ifstream M_par_w_int;
	M_par_w_int.open(PATH + "M_without_interaction_alpha0.0_sigma3.0_chil1.0_hx0.0500_ht0.0000.txt");
	std::ofstream Relaxation;
	Relaxation.open(PATH + "Relaxation_times_int_sigma_3.txt", std::fstream::app);
	float im_chi_int;
	float im_chi;
	float re_chi_int;
	float re_chi;
	float im_max = -1;
	float relax_time;
	for (int i = 0; i < deg_num; ++i)
	{

		calc(params, r[i], t, im_chi, im_chi_int, re_chi, re_chi_int, i, M_par_int, M_par_w_int, M_perp_int, M_perp_w_int);
		printf("%f\n", r_2[i]);
		fprintf(Re_Xi_file_without_interaction, "%f %f\n", r_2[i], re_chi); // abscissa is tau*omega/2
		if (im_max < im_chi_int)
		{
			im_max = im_chi_int;
			relax_time = r_2[i];
		}
		fprintf(Im_Xi_file_without_interaction, "%f %f\n", r_2[i], im_chi);
		fprintf(Re_Xi_file_interaction, "%f %f\n", r_2[i], re_chi_int); // abscissa is tau*omega/2
		fprintf(Im_Xi_file_interaction, "%f %f\n", r_2[i], im_chi_int);
	}
	//	Relaxation << params.xi * 180 / PI << ' ' << relax_time << '\n'; // uncomment when want to count relaxation times
	//  fprintf(Relaxation_times, "%f %f\n", params.xi * 180 / PI, r_2[std::distance(im_chi_values.begin(), std::find(im_chi_values.begin(), im_chi_values.end(), im_max))]);
	M_perp_int.close();
	M_perp_w_int.close();
	M_par_int.close();
	M_par_w_int.close();
	Relaxation.close();
	fclose(Re_Xi_file_without_interaction);
	fclose(Im_Xi_file_without_interaction);
	fclose(Re_Xi_file_interaction);
	fclose(Im_Xi_file_interaction);

	// fclose(Relaxation_times);
	delete[] d;
	delete[] r;
	delete[] r_2;
}

int main(int argc, char *argv[])
{
	Params params;
	if (argc < 8)
		printf("Inditicate alpha, sigma, XiL\n");
	if (argc >= 8)
	{
		// std::string vzaim_str = argv[3];
		sscanf(argv[1], "%f", &params.alp);
		sscanf(argv[2], "%f", &params.sigma);
		sscanf(argv[3], "%f", &params.XiL);
		sscanf(argv[4], "%f", &params.ht);
		sscanf(argv[5], "%f", &params.tau_omega0);
		sscanf(argv[6], "%f", &params.tau_omegaN);
		sscanf(argv[7], "%f", &params.xi); // в градусах
		sscanf(argv[8], "%f", &params.psi);
	}
	params.xi = params.xi * PI / 180;
	params.psi = params.psi * PI / 180;
	params.sin_xi = sin(params.xi);
	params.cos_xi = cos(params.xi);
	params.sin_psi = sin(params.psi);
	params.cos_psi = cos(params.psi);
	params.sigma2 = params.sigma * 2;
	// for (auto deg : {0, 30, 45, 60, 90})
	{
		//	params.xi = deg * PI / 180;
		//	params.sin_xi = sin(params.xi);
		//	params.cos_xi = cos(params.xi);
		init_and_set(params);
	}
	return 0;
}