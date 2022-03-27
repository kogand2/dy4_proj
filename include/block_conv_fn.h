/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_BLOCKCONVFN_H
#define DY4_BLOCKCONVFN_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>


// declaration of a function prototypes
std::vector<float> fmDemod(std::vector<float>, std::vector<float>, std::vector<float> &);
std::vector<float> fmDemodArctan(std::vector<float>, std::vector<float>, float &);

void low_pass_coeff(float, float, int, std::vector<float> &);

void band_pass_coeff(float, float, float, int, std::vector<float> &);

void all_pass_coeff(std::vector<float> &,std::vector<float> &, std::vector<float> &);

void rrc_coeff(float, int, std::vector<float> &)

void fmPll(std::vector<float> &, std::vector<float> &, std::vector<float> &, float, float, float, float, float);

void mixer(std::vector<float> &, std::vector<float> &, std::vector<float> &);

void sq_non_linearity(std::vector<float> &)

void state_block_conv(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, std::vector<float> &);

void ds_block_conv(std::vector<float> &, const std::vector<float> , const std::vector<float> , std::vector<float> &, int, std::vector<float> &);

void rs_block_conv(std::vector<float> &, const std::vector<float> , const std::vector<float> , std::vector<float> &, int, int, std::vector<float> &);

#endif // DY4_BLOCKCONVFN_H
