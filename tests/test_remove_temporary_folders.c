/*
 * Regression tests for Remove_temporary_folders (src/data_sparse.c), issue #12.
 *
 * The per-beamlet temporary folders are named from the config's Output_Directory, so the path
 * must never reach a shell. Before the fix it was quoted into `rm -r "..."` and run by system():
 * a directory name containing a double quote and $(...) ran the command and left the folder.
 *
 * Built by `make test_remove_tmp`.
 */
#include "include/define.h"
#include "include/data_sparse.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static int failures = 0;

#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); failures++; } } while (0)

static int exists(const char *p) { struct stat st; return lstat(p, &st) == 0; }

static void write_file(const char *p) {
  FILE *f = fopen(p, "w");
  if (f == NULL) { perror(p); exit(2); }
  fputs("x\n", f);
  fclose(f);
}

static void make_dir(const char *p) {
  if (mkdir(p, 0700) != 0) { perror(p); exit(2); }
}

int main(void) {
  char base[] = "/tmp/mcsq_rmtmp_XXXXXX";
  if (mkdtemp(base) == NULL) { perror("mkdtemp"); return 2; }
  if (chdir(base) != 0) { perror("chdir"); return 2; }  // a relative `touch PWNED` would land here

  // Kept outside the temporary folders; reached from inside them only through symlinks.
  char keep_file[PATH_SIZE], keep_dir[PATH_SIZE], keep_dir_file[PATH_SIZE];
  snprintf(keep_file, sizeof keep_file, "%s/keep.txt", base);
  snprintf(keep_dir, sizeof keep_dir, "%s/keepdir", base);
  snprintf(keep_dir_file, sizeof keep_dir_file, "%s/keepdir/inside.txt", base);
  write_file(keep_file);
  make_dir(keep_dir);
  write_file(keep_dir_file);

  // InputPath as the caller builds it ("<Output_Directory>tmp/Beamlet_"), with a hostile name.
  char input[PATH_SIZE];
  snprintf(input, sizeof input, "%s/Beamlet_\"$(touch PWNED)\"_", base);

  char d1[PATH_SIZE], d2[PATH_SIZE], p[PATH_SIZE];
  snprintf(d1, sizeof d1, "%s1", input);
  snprintf(d2, sizeof d2, "%s2", input);
  make_dir(d1);
  make_dir(d2);
  snprintf(p, sizeof p, "%s/Sparse_Dose.bin", d1); write_file(p);
  snprintf(p, sizeof p, "%s/nested", d1); make_dir(p);
  snprintf(p, sizeof p, "%s/nested/deep.txt", d1); write_file(p);
  snprintf(p, sizeof p, "%s/link_to_file", d1);
  if (symlink(keep_file, p) != 0) { perror("symlink"); return 2; }
  snprintf(p, sizeof p, "%s/link_to_dir", d2);
  if (symlink(keep_dir, p) != 0) { perror("symlink"); return 2; }

  // Three folders requested, two exist: the missing one must not stop the others.
  Remove_temporary_folders(input, 3);

  CHECK(!exists(d1), "folder 1 with nested content removed");
  CHECK(!exists(d2), "folder 2 removed");
  CHECK(!exists("PWNED"), "no command from the path was executed");
  CHECK(exists(keep_file), "symlinked file outside the folders kept");
  CHECK(exists(keep_dir_file), "symlinked directory's contents outside the folders kept");

  // Clean up what the test created (best effort).
  remove(keep_dir_file); rmdir(keep_dir); remove(keep_file); remove("PWNED");
  if (exists(d1) || exists(d2)) printf("note: leftovers under %s\n", base);
  else rmdir(base);

  if (failures == 0) printf("PASS: Remove_temporary_folders\n");
  return failures == 0 ? 0 : 1;
}
