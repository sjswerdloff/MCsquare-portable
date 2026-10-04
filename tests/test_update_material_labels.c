/*
 * Regression tests for Update_material_labels (src/compute_geometry.c), issue #20.
 *
 * In SPR transport every lane whose voxel material differs from its current label must get the
 * new label, and the stopping power of every lane is then recomputed from the labels. The
 * converted loop reused the lane index in its inner loop, so only the first changed lane was
 * relabelled.
 *
 * Built by `make test_material_labels` with UBSan.
 */
#include "include/define.h"
#include "include/struct.h"
#include "include/compute_geometry.h"
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); failures++; } } while (0)

enum { NVOX = VLENGTH, NDATA = 4 };

static unsigned short voxel_material[NVOX];
static VAR_DATA sp_water[NDATA] = {1.0f, 1.1f, 1.2f, 1.3f};
static VAR_DATA sp_bone[NDATA] = {2.0f, 2.1f, 2.2f, 2.3f};

static void setup(DATA_CT *ct, Materials material[2], int *v_index, int *v_data_index, int *v_label, VAR_COMPUTE *v_sp) {
  memset(ct, 0, sizeof(*ct));
  memset(material, 0, 2 * sizeof(Materials));
  ct->material = voxel_material;
  material[0].Stop_Pow = sp_water;
  material[1].Stop_Pow = sp_bone;
  for (int i = 0; i < VLENGTH; i++) {
    voxel_material[i] = 0;
    v_index[i] = i;
    v_data_index[i] = i % NDATA;
    v_label[i] = 0;
    v_sp[i] = -1.0f;  // sentinel: shows whether the stopping power was recomputed
  }
}

int main(void) {
  DATA_CT ct;
  Materials material[2];
  int v_index[VLENGTH], v_data_index[VLENGTH], v_label[VLENGTH];
  VAR_COMPUTE v_sp[VLENGTH];

  // Two lanes enter bone in the same pass, and neither is the last lane.
  setup(&ct, material, v_index, v_data_index, v_label, v_sp);
  voxel_material[2] = 1;
  voxel_material[5] = 1;
  int changed = Update_material_labels(&ct, material, v_index, v_data_index, v_label, v_sp);
  CHECK(changed == 2, "both changed lanes counted");
  CHECK(v_label[2] == 1, "first changed lane relabelled");
  CHECK(v_label[5] == 1, "second changed lane relabelled");
  int others_kept = 1;
  for (int i = 0; i < VLENGTH; i++) if (i != 2 && i != 5 && v_label[i] != 0) others_kept = 0;
  CHECK(others_kept, "unchanged lanes keep their label");
  CHECK(v_sp[5] == sp_bone[5 % NDATA], "second changed lane gets bone stopping power");
  CHECK(v_sp[0] == sp_water[0], "every lane's stopping power recomputed after a change");

  // No lane changes material: labels and stopping power are left alone.
  setup(&ct, material, v_index, v_data_index, v_label, v_sp);
  changed = Update_material_labels(&ct, material, v_index, v_data_index, v_label, v_sp);
  CHECK(changed == 0, "no change counted");
  CHECK(v_sp[3] == -1.0f, "stopping power untouched when nothing changed");

  // The last lane alone changes: the boundary of the lane loop.
  setup(&ct, material, v_index, v_data_index, v_label, v_sp);
  voxel_material[VLENGTH - 1] = 1;
  changed = Update_material_labels(&ct, material, v_index, v_data_index, v_label, v_sp);
  CHECK(changed == 1 && v_label[VLENGTH - 1] == 1, "last lane relabelled");

  if (failures == 0) printf("PASS: Update_material_labels\n");
  return failures == 0 ? 0 : 1;
}
