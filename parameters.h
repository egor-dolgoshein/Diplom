#include <cmath>
#include <vector>

const float hx = 0.05f;
const float ht = 0.01f;
const int deg_num = 41;

const float hy = 0.01f;
const float hy_2 = hy * hy;
const float hy__2 = hy * 0.5f;

const float PI = std::atan(1.0f) * 4.f;
const float hx_2 = hx * hx;
const float hx__2 = hx * 0.5f;
const float ht__2 = ht * 0.5f;

struct Params
{
	float sigma;
	float alp;
	float XiL;

	float xi;
	float psi;
	float ht;		  // step of time
	float tau_omega0; // start frequency
	float tau_omegaN; // start frequency
	float sigma2;
	float sin_xi;
	float sin_psi;
	float cos_xi;
	float cos_psi;
};
