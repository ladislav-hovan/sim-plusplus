#pragma once

struct ParamsLJ
{
	double epsilon = 119.8;  // Lennard-Jones energy parameter, in K
	double r_m = 1.12246f * 0.3405;  // Lennard-Jones distance at minimum, in nm
	double cutoff = 2.5f * 0.3405;  // Cutoff distance of potential function, in nm
};