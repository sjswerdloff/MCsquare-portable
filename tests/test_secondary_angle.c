/*
 * Regression test for issue #16: the emission angle of a nuclear-inelastic secondary proton is
 * sampled from the angular table interpolated at the secondary's table energy (MeV).
 *
 * The code used to interpolate the table at the energy after rescaling to eV, extrapolating about
 * 1e6 bracket-widths past the table. This synthetic table gives angular weight only to the
 * 30-60 degree points, whose sampling bins span 25-65 degrees, so every emitted angle must fall in
 * [25, 65]. Row 2 is smaller than row 1, which made the old cumulative negative, so the old code
 * sent every sample past the table and emitted it at 165-180 degrees.
 *
 * Built by `make test_secondary_angle`.
 */
#include "include/define.h"
#include "include/struct.h"
#include "include/compute_nuclear_interaction.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

enum { ROWS = 3, NSAMPLES = 2000 };

static VAR_DATA energy_list[ROWS] = {0.0f, 20.0f, 40.0f};                 // secondary energies (MeV)
static VAR_DATA d_cross[ROWS] = {0.0f, 0.0f, 1.0f};                       // secondary energy always in row 1..2
static VAR_DATA dd[ROWS * 13];                                           // angular table, 13 angles per row
static VAR_DATA incident[2] = {100.0f, 110.0f};                           // incident energies (MeV)

int main(void) {
  memset(dd, 0, sizeof dd);
  for (int k = 3; k <= 6; k++) {   // 30, 40, 50, 60 degrees only
    dd[13 * 1 + k] = 4.0f;         // row 1 (20 MeV)
    dd[13 * 2 + k] = 1.0f;         // row 2 (40 MeV)
  }
  DATA_Nuclear_Inelastic blocks[2];
  memset(blocks, 0, sizeof blocks);
  for (int b = 0; b < 2; b++) {
    blocks[b].Proton_Mult = 1.0f;
    blocks[b].P_Nbr_data = ROWS;
    blocks[b].P_Energy_list = energy_list;
    blocks[b].P_D_Cross_section = d_cross;
    blocks[b].P_DD_Cross_section = dd;
  }
  Materials material;
  memset(&material, 0, sizeof material);
  material.Nbr_Inelastic_Energy = 2;
  material.Inelastic_Energy_List = incident;
  material.Nuclear_Inelastic = blocks;

  DATA_config config;
  memset(&config, 0, sizeof config);
  config.Ecut_Pro = 0.0f;
  config.Simulate_Secondary_Protons = 1;

  static Hadron hadron;
  memset(&hadron, 0, sizeof hadron);
  hadron.v_T[0] = 105.0f * UMeV;
  hadron.v_M[0] = 1.0f;
  hadron.v_u[0] = 0.6f; hadron.v_v[0] = 0.0f; hadron.v_w[0] = 0.8f;  // not along an axis

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 12345u, 0u);

  int outside = 0, backward = 0;
  double min_deg = 180, max_deg = 0;
  for (int n = 0; n < NSAMPLES; n++) {
    Hadron_buffer secondary[1];
    memset(secondary, 0, sizeof secondary);
    int nsec = 0;
    Compute_Nuclear_Inelastic_proton(0, &hadron, secondary, &nsec, &material, 0, &rng, &config);
    if (nsec != 1) { printf("FAIL: no secondary emitted (sample %d)\n", n); return 1; }
    double dot = secondary[0].u * hadron.v_u[0] + secondary[0].v * hadron.v_v[0] + secondary[0].w * hadron.v_w[0];
    double deg = acos(fmax(-1.0, fmin(1.0, dot))) * 180.0 / M_PI;
    if (deg < min_deg) min_deg = deg;
    if (deg > max_deg) max_deg = deg;
    if (deg < 25.0 - 0.01 || deg > 65.0 + 0.01) outside++;
    if (deg > 165.0) backward++;
  }
  printf("emitted angles %.1f to %.1f deg; %d of %d outside [25, 65]; %d beyond 165\n", min_deg, max_deg, outside, NSAMPLES, backward);
  if (outside != 0) { printf("FAIL: angles not drawn from the table\n"); return 1; }
  printf("PASS: secondary emission angle follows the angular table\n");
  return 0;
}
