/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

Modified 2026: Replaced Intel MKL VSL with PCG random number generator for portability.

Algorithm references:
  - PCG32: Melissa O'Neill, "PCG: A Family of Simple Fast Space-Efficient
    Statistically Good Algorithms for Random Number Generation" (2014)
    https://www.pcg-random.org/
  - Box-Muller transform: G. E. P. Box and Mervin E. Muller, "A Note on the
    Generation of Random Normal Deviates" (1958), Annals of Mathematical
    Statistics 29(2):610-611
  - Float conversion: Standard technique using upper bits of uint32
  - Edge case handling: Informed by GSL (GNU Scientific Library) gauss.c
    which guards against log(0) and division by zero. Our implementation
    clamps uniform values away from 0 to ensure log() domain validity.
    See: https://github.com/ampl/gsl/blob/master/randist/gauss.c
*/


#include "include/compute_random.h"
#include <math.h>

/*
 * Convert PCG uint32 to float in (0, 1) exclusive of endpoints.
 * Uses upper 24 bits of randomness (matches float mantissa precision).
 *
 * Edge case handling: Clamps output to [FLT_EPSILON, 1-FLT_EPSILON] to
 * prevent log(0) in Box-Muller transform. This follows the same principle
 * as GSL's gsl_rng_uniform_pos() which returns values in (0,1].
 */
static inline VAR_COMPUTE pcg32_to_float(pcg32_random_t* rng) {
    uint32_t bits = pcg32_random_r(rng);
    VAR_COMPUTE val = (bits >> 8) * (1.0f / 16777216.0f);
    // Clamp to (0, 1) to ensure log() domain validity in Box-Muller
    if (val < FLT_EPSILON) val = FLT_EPSILON;
    if (val > (1.0f - FLT_EPSILON)) val = 1.0f - FLT_EPSILON;
    return val;
}

#if VAR_COMPUTE_PRECISION != 1
/*
 * Double precision version using two 32-bit values for 53-bit mantissa.
 */
static inline VAR_COMPUTE pcg32_to_double(pcg32_random_t* rng) {
    uint64_t bits = ((uint64_t)pcg32_random_r(rng) << 32) | pcg32_random_r(rng);
    VAR_COMPUTE val = (bits >> 11) * (1.0 / 9007199254740992.0);  // 2^53
    if (val < DBL_EPSILON) val = DBL_EPSILON;
    if (val > (1.0 - DBL_EPSILON)) val = 1.0 - DBL_EPSILON;
    return val;
}
#endif

/*
 * Box-Muller transform for generating standard normal random variates.
 *
 * Algorithm: Given u1, u2 uniform in (0,1):
 *   z0 = sqrt(-2 * ln(u1)) * cos(2 * pi * u2)
 *   z1 = sqrt(-2 * ln(u1)) * sin(2 * pi * u2)
 *
 * Reference: Box & Muller (1958), "A Note on the Generation of Random
 * Normal Deviates", Annals of Mathematical Statistics 29(2):610-611
 *
 * Edge cases: u1 is clamped away from 0 by pcg32_to_float() to prevent
 * log(0) = -infinity. This follows GSL's approach of ensuring the log
 * argument is strictly positive.
 */
static VAR_COMPUTE box_muller_single(pcg32_random_t* rng) {
    VAR_COMPUTE u1, u2, z0;

    #if VAR_COMPUTE_PRECISION==1
        u1 = pcg32_to_float(rng);
        u2 = pcg32_to_float(rng);
        z0 = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (VAR_COMPUTE)M_PI * u2);
    #else
        u1 = pcg32_to_double(rng);
        u2 = pcg32_to_double(rng);
        z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    #endif

    return z0;
}


void rand_uniform(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd){

  __assume_aligned(v_rnd, 64);

  #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    #if VAR_COMPUTE_PRECISION==1
      v_rnd[i] = pcg32_to_float(stream);
    #else
      v_rnd[i] = pcg32_to_double(stream);
    #endif
  }

  return;
}


VAR_COMPUTE single_rand_uniform(VSLStreamStatePtr stream){

  #if VAR_COMPUTE_PRECISION==1
    return pcg32_to_float(stream);
  #else
    return pcg32_to_double(stream);
  #endif
}


void rand_normal(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_mu, VAR_COMPUTE *v_sigma){

  __assume_aligned(v_rnd, 64);
  __assume_aligned(v_mu, 64);
  __assume_aligned(v_sigma, 64);

  // Generate standard normal values using Box-Muller
  for (int i = 0; i < VLENGTH; i += 2) {
    VAR_COMPUTE u1, u2, z0, z1;

    #if VAR_COMPUTE_PRECISION==1
      u1 = pcg32_to_float(stream);
      u2 = pcg32_to_float(stream);
      VAR_COMPUTE r = sqrtf(-2.0f * logf(u1));
      VAR_COMPUTE theta = 2.0f * (VAR_COMPUTE)M_PI * u2;
      z0 = r * cosf(theta);
      z1 = r * sinf(theta);
    #else
      u1 = pcg32_to_double(stream);
      u2 = pcg32_to_double(stream);
      VAR_COMPUTE r = sqrt(-2.0 * log(u1));
      VAR_COMPUTE theta = 2.0 * M_PI * u2;
      z0 = r * cos(theta);
      z1 = r * sin(theta);
    #endif

    v_rnd[i] = z0;
    if (i + 1 < VLENGTH) {
      v_rnd[i + 1] = z1;
    }
  }

  // Scale by sigma and shift by mu
  #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_rnd[i] = v_sigma[i] * v_rnd[i] + v_mu[i];
  }

  return;
}


void rand_normal_zero(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_sigma){

  __assume_aligned(v_rnd, 64);
  __assume_aligned(v_sigma, 64);

  // Generate standard normal values using Box-Muller
  for (int i = 0; i < VLENGTH; i += 2) {
    VAR_COMPUTE u1, u2, z0, z1;

    #if VAR_COMPUTE_PRECISION==1
      u1 = pcg32_to_float(stream);
      u2 = pcg32_to_float(stream);
      VAR_COMPUTE r = sqrtf(-2.0f * logf(u1));
      VAR_COMPUTE theta = 2.0f * (VAR_COMPUTE)M_PI * u2;
      z0 = r * cosf(theta);
      z1 = r * sinf(theta);
    #else
      u1 = pcg32_to_double(stream);
      u2 = pcg32_to_double(stream);
      VAR_COMPUTE r = sqrt(-2.0 * log(u1));
      VAR_COMPUTE theta = 2.0 * M_PI * u2;
      z0 = r * cos(theta);
      z1 = r * sin(theta);
    #endif

    v_rnd[i] = z0;
    if (i + 1 < VLENGTH) {
      v_rnd[i + 1] = z1;
    }
  }

  // Scale by sigma (mu = 0)
  #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_rnd[i] = v_sigma[i] * v_rnd[i];
  }

  return;
}


VAR_COMPUTE single_rand_normal(VSLStreamStatePtr stream, VAR_COMPUTE mu, VAR_COMPUTE sigma){

  VAR_COMPUTE rnd = box_muller_single(stream);
  rnd = sigma * rnd + mu;

  return rnd;
}
