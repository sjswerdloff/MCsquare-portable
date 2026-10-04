/*
 * Regression tests for Transport_to_CT (src/compute_beam_model.c), issue #8.
 *
 * A primary generated outside the CT is translated along its direction to the first CT face it
 * enters, losing energy in air on the way. A primary whose ray never enters the CT must be left
 * where it is (outside), with its energy untouched, so that Generate_PBS_particle counts it as
 * "generated outside the geometry" and drops it. Before the fix, that case read Translation[3]
 * out of bounds and could copy an uninitialised position.
 *
 * Built by `make test_transport` with UBSan, so an out-of-bounds read fails the run
 * deterministically rather than depending on stack contents.
 */
#include "include/define.h"
#include "include/struct.h"
#include "include/compute_beam_model.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); failures++; } } while (0)

static VAR_DATA CT[3] = {100.0f, 100.0f, 100.0f};
static DATA_config config;

static Hadron_buffer proton(float x, float y, float z, float u, float v, float w) {
  Hadron_buffer h;
  memset(&h, 0, sizeof(h));
  h.x = x; h.y = y; h.z = z;
  h.u = u; h.v = v; h.w = w;
  h.T = 100.0f * UMeV;
  h.M = 1; h.charge = 1; h.mass = 1;
  h.type = Proton;
  return h;
}

static void expect_unchanged(Hadron_buffer before, Hadron_buffer after, const char *what) {
  char msg[160];
  snprintf(msg, sizeof msg, "%s: position unchanged", what);
  CHECK(after.x == before.x && after.y == before.y && after.z == before.z, msg);
  snprintf(msg, sizeof msg, "%s: energy unchanged", what);
  CHECK(after.T == before.T, msg);
  snprintf(msg, sizeof msg, "%s: type unchanged", what);
  CHECK(after.type == before.type, msg);
}

/* At x = -10 moving along +y: x never changes, so the CT is never entered. The x and z
 * translations divide by a zero direction component and are +inf, so the loop tries those faces,
 * gets NaN positions, and finds no entry. */
static void test_ray_parallel_to_faces_misses(void) {
  Hadron_buffer h = proton(-10.0f, 50.0f, 50.0f, 0.0f, 1.0f, 0.0f);
  Hadron_buffer before = h;
  Transport_to_CT(&h, CT, &config);
  expect_unchanged(before, h, "parallel miss");
}

/* Outside on every axis and moving away on every axis: all three translations are negative, so
 * no face is ever tried. */
static void test_ray_moving_away_misses(void) {
  float c = (float)(1.0 / sqrt(3.0));
  Hadron_buffer h = proton(-10.0f, -10.0f, -10.0f, -c, -c, -c);
  Hadron_buffer before = h;
  Transport_to_CT(&h, CT, &config);
  expect_unchanged(before, h, "moving-away miss");
}

static void test_inside_is_untouched(void) {
  Hadron_buffer h = proton(50.0f, 50.0f, 50.0f, 0.0f, 0.0f, 1.0f);
  Hadron_buffer before = h;
  Transport_to_CT(&h, CT, &config);
  expect_unchanged(before, h, "already inside");
}

/* Enters through the z = 0 face: lands just inside, and loses energy in proportion to the air gap. */
static void test_ray_hits_z_face(void) {
  Hadron_buffer near_h = proton(50.0f, 50.0f, -10.0f, 0.0f, 0.0f, 1.0f);
  Hadron_buffer far_h = proton(50.0f, 50.0f, -20.0f, 0.0f, 0.0f, 1.0f);
  Transport_to_CT(&near_h, CT, &config);
  Transport_to_CT(&far_h, CT, &config);

  CHECK(near_h.x == 50.0f && near_h.y == 50.0f, "hit: lateral position kept");
  CHECK(near_h.z > 0.0f && near_h.z < 0.01f, "hit: placed just inside the z = 0 face");
  CHECK(near_h.type == Proton, "hit: still a proton");

  double loss_near = 100.0 * UMeV - near_h.T;
  double loss_far = 100.0 * UMeV - far_h.T;
  CHECK(loss_near > 0.0, "hit: energy lost in air");
  CHECK(fabs(loss_far / loss_near - 2.0) < 0.01, "hit: air loss proportional to path length");
}

int main(void) {
  memset(&config, 0, sizeof(config));
  config.Ecut_Pro = 0.5f;

  test_ray_parallel_to_faces_misses();
  test_ray_moving_away_misses();
  test_inside_is_untouched();
  test_ray_hits_z_face();

  if (failures) {
    printf("%d check(s) failed\n", failures);
    return 1;
  }
  printf("Transport_to_CT: all checks passed\n");
  return 0;
}
