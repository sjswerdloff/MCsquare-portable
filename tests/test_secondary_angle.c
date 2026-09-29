/*
 * Regression test for issue #16: the emission angle of a nuclear-inelastic secondary (proton,
 * deuteron, alpha) is sampled from its own angular table interpolated at the secondary's table
 * energy (MeV).
 *
 * The code used to interpolate the table at the energy after rescaling to eV, extrapolating about
 * 1e6 bracket-widths past the table. Each species here gets angular weight on different ICRU
 * angles, so its emitted angles must fall inside those points' sampling bins:
 *   proton    30-60 deg   -> [25, 65]
 *   deuteron  90-110 deg  -> [80, 120]
 *   alpha     10 deg      -> [5, 15]
 * The supports do not overlap, so reading another species' table also fails. Row 2 is smaller than
 * row 1, which made the old cumulative negative, so the old code sent every sample past the table
 * and emitted it at 165-180 degrees.
 *
 * The last checks (issue #24): an out-of-range angle index must abort, not be clamped to an angle. For
 * each species a forked child is given a table whose cumulative decreases, which drives the index to 13
 * even with the fix; the test passes only if the child dies of SIGABRT with the range guard's own
 * "FATAL <function>: angle_index 13" line on stderr.
 *
 * Built by `make test_angle`.
 */
#include "include/define.h"
#include "include/struct.h"
#include "include/compute_nuclear_interaction.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

enum { ROWS = 3, NSAMPLES = 2000 };

typedef VAR_COMPUTE (*inelastic_fn)(int, Hadron *, Hadron_buffer *, int *, Materials *, int, VSLStreamStatePtr, DATA_config *);

static VAR_DATA energy_list[ROWS] = {0.0f, 20.0f, 40.0f};  // secondary energies (MeV)
static VAR_DATA d_cross[ROWS] = {0.0f, 0.0f, 1.0f};        // secondary energy always in row 1..2
static VAR_DATA dd_p[ROWS * 13], dd_d[ROWS * 13], dd_a[ROWS * 13];
static VAR_DATA incident[2] = {100.0f, 110.0f};            // incident energies (MeV)

// Give angular weight to ICRU angle indices first..last in rows 1 and 2, row 2 smaller.
static void fill(VAR_DATA *dd, int first, int last) {
  memset(dd, 0, ROWS * 13 * sizeof(VAR_DATA));
  for (int k = first; k <= last; k++) {
    dd[13 * 1 + k] = 4.0f;
    dd[13 * 2 + k] = 1.0f;
  }
}

static int check(const char *name, inelastic_fn fn, int expected_type, double lo, double hi,
                 Materials *material, DATA_config *config) {
  static Hadron hadron;
  memset(&hadron, 0, sizeof hadron);
  hadron.v_T[0] = 105.0f * UMeV;
  hadron.v_M[0] = 1.0f;
  hadron.v_u[0] = 0.6f; hadron.v_v[0] = 0.0f; hadron.v_w[0] = 0.8f;  // not along an axis

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 12345u, 0u);

  int outside = 0;
  double min_deg = 180, max_deg = 0;
  for (int n = 0; n < NSAMPLES; n++) {
    Hadron_buffer secondary[1];
    memset(secondary, 0, sizeof secondary);
    int nsec = 0;
    fn(0, &hadron, secondary, &nsec, material, 0, &rng, config);
    if (nsec != 1) { printf("FAIL %s: no secondary emitted (sample %d)\n", name, n); return 1; }
    if ((int)secondary[0].type != expected_type) { printf("FAIL %s: wrong secondary type\n", name); return 1; }
    double dot = secondary[0].u * hadron.v_u[0] + secondary[0].v * hadron.v_v[0] + secondary[0].w * hadron.v_w[0];
    double deg = acos(fmax(-1.0, fmin(1.0, dot))) * 180.0 / M_PI;
    if (deg < min_deg) min_deg = deg;
    if (deg > max_deg) max_deg = deg;
    if (deg < lo - 0.01 || deg > hi + 0.01) outside++;
  }
  printf("%-8s emitted angles %.1f to %.1f deg; %d of %d outside [%.0f, %.0f]\n",
         name, min_deg, max_deg, outside, NSAMPLES, lo, hi);
  if (outside != 0) { printf("FAIL %s: angles not drawn from its own table\n", name); return 1; }
  return 0;
}

// Run fn on a broken angular table in a forked child. Passes only if the child dies of SIGABRT AND its
// stderr carries the range guard's own diagnostic, "FATAL <function>: angle_index 13", so an abort from
// anywhere else does not count.
static int expect_abort(const char *name, inelastic_fn fn, int species, DATA_Nuclear_Inelastic *blocks,
                        VAR_DATA *dd_bad, Materials *material, DATA_config *config) {
  int fds[2];
  if (pipe(fds) != 0) { perror("pipe"); return 1; }
  fflush(stdout);
  pid_t pid = fork();
  if (pid == 0) {
    close(fds[0]);
    dup2(fds[1], STDERR_FILENO);
    for (int b = 0; b < 2; b++) {
      if (species == 0) blocks[b].P_DD_Cross_section = dd_bad;
      if (species == 1) blocks[b].D_DD_Cross_section = dd_bad;
      if (species == 2) blocks[b].A_DD_Cross_section = dd_bad;
    }
    Hadron_buffer out[1];
    int nsec;
    static Hadron h;
    memset(&h, 0, sizeof h);
    h.v_T[0] = 105.0f * UMeV; h.v_M[0] = 1.0f; h.v_w[0] = 1.0f;
    pcg32_random_t r;
    pcg32_srandom_r(&r, 1u, 0u);
    for (int n = 0; n < NSAMPLES; n++) { nsec = 0; fn(0, &h, out, &nsec, material, 0, &r, config); }
    _exit(0);  // reached only if no sample ever produced an out-of-range index
  }
  close(fds[1]);
  char err[4096] = {0};
  size_t len = 0;
  ssize_t got;
  while (len < sizeof err - 1 && (got = read(fds[0], err + len, sizeof err - 1 - len)) > 0) len += (size_t)got;
  close(fds[0]);
  int status = 0;
  waitpid(pid, &status, 0);

  const char *fname = species == 0 ? "Compute_Nuclear_Inelastic_proton"
                    : species == 1 ? "Compute_Nuclear_Inelastic_deuteron" : "Compute_Nuclear_Inelastic_alpha";
  char want[128];
  snprintf(want, sizeof want, "FATAL %s: angle_index 13", fname);
  int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
  if (aborted && strstr(err, want)) {
    printf("broken table (%s): aborts at the range guard (SIGABRT, \"%s\")\n", name, want);
    return 0;
  }
  printf("FAIL broken table (%s): expected SIGABRT with \"%s\", got %s %d, stderr: %s\n", name, want,
         WIFSIGNALED(status) ? "signal" : "exit", WIFSIGNALED(status) ? WTERMSIG(status) : WEXITSTATUS(status),
         err[0] ? err : "(empty)");
  return 1;
}

int main(void) {
  fill(dd_p, 3, 6);  // 30, 40, 50, 60 deg
  fill(dd_d, 8, 9);  // 90, 110 deg
  fill(dd_a, 1, 1);  // 10 deg

  DATA_Nuclear_Inelastic blocks[2];
  memset(blocks, 0, sizeof blocks);
  for (int b = 0; b < 2; b++) {
    blocks[b].Proton_Mult = 1.0f;
    blocks[b].P_Nbr_data = ROWS;
    blocks[b].P_Energy_list = energy_list;
    blocks[b].P_D_Cross_section = d_cross;
    blocks[b].P_DD_Cross_section = dd_p;
    blocks[b].Deuteron_Mult = 1.0f;
    blocks[b].D_Nbr_data = ROWS;
    blocks[b].D_Energy_list = energy_list;
    blocks[b].D_D_Cross_section = d_cross;
    blocks[b].D_DD_Cross_section = dd_d;
    blocks[b].Alpha_Mult = 1.0f;
    blocks[b].A_Nbr_data = ROWS;
    blocks[b].A_Energy_list = energy_list;
    blocks[b].A_D_Cross_section = d_cross;
    blocks[b].A_DD_Cross_section = dd_a;
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
  config.Simulate_Secondary_Deuterons = 1;
  config.Simulate_Secondary_Alphas = 1;

  // A broken table: negative weights, so the cumulative decreases and angle_index reaches 13.
  static VAR_DATA dd_bad[ROWS * 13];
  memset(dd_bad, 0, sizeof dd_bad);
  for (int k = 3; k <= 12; k++) { dd_bad[13 * 1 + k] = -4.0f; dd_bad[13 * 2 + k] = -1.0f; }

  int failures = 0;
  failures += check("proton", Compute_Nuclear_Inelastic_proton, Proton, 25.0, 65.0, &material, &config);
  failures += check("deuteron", Compute_Nuclear_Inelastic_deuteron, Deuteron, 80.0, 120.0, &material, &config);
  failures += check("alpha", Compute_Nuclear_Inelastic_alpha, Alpha, 5.0, 15.0, &material, &config);
  failures += expect_abort("proton", Compute_Nuclear_Inelastic_proton, 0, blocks, dd_bad, &material, &config);
  failures += expect_abort("deuteron", Compute_Nuclear_Inelastic_deuteron, 1, blocks, dd_bad, &material, &config);
  failures += expect_abort("alpha", Compute_Nuclear_Inelastic_alpha, 2, blocks, dd_bad, &material, &config);
  if (failures) return 1;
  printf("PASS: secondary emission angles follow each species' angular table; out-of-range index aborts\n");
  return 0;
}
