#pragma once

namespace Constants
{
	inline constexpr double boltzmann = 1.38064852e-23;  // J/K
	inline constexpr double mu = 1.660539040e-27;  // kg
	inline constexpr double KToNatural = boltzmann / (1e6 * mu);
}
