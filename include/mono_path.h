/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_MONO_H
#define DY4_MONO_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>

// declaration of a function prototypes
void readStdinBlockData(unsigned int, std::vector<float> &);

void downsample(int, std::vector<float>, std::vector<float> &);

std::vector<float> mono_path(std::vector<float> , std::vector<float> , std::vector<float> &, int);

#endif // DY4_MONO_H
