#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

int main()
{
    float tau_omega0 = 0.01;
    int deg_num = 41;
    float *d = new float[deg_num];
    float *r = new float[deg_num];
    float *r_2 = new float[deg_num];
    r[0] = 2 * tau_omega0;
    r_2[0] = 0.5f * r[0];
    d[0] = std::log10(r[0]);
    for (int i = 1; i < deg_num; i++)
    {
        d[i] = d[0] + 0.1f * i;
        r[i] = pow(10.f, d[i]);
        r_2[i] = r[i] * 0.5f;
        // printf("%f\t%f\t%f\n", d[i], r[i], r_2[i]);
    }
    std::vector<float> im_chi_values = {1, 2, 3, 4, 5, 456, 234, 523, 45, 1234, 234, 5, 1234, 123};
    float im_max = -1;
    for (auto v : im_chi_values)
        if (im_max < v)
            im_max = v;
    std::cout << '\n'
              << std::distance(im_chi_values.begin(), std::find(im_chi_values.begin(), im_chi_values.end(), im_max));
}