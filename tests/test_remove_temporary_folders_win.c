/*
 * Windows regression test for Remove_temporary_folders (src/data_sparse.c), issues #12 and #28.
 *
 * Everything happens under a root the test creates itself: argv[1], which must not exist yet.
 *   <root>\tmp_1\a\b\f.txt      nested contents              -> removed, tmp_1 removed
 *   <root>\tmp_1\j              junction to <root>\keep      -> the junction is removed, NOT followed:
 *   <root>\keep\sentinel.txt    external retained target        keep\ and sentinel.txt must survive
 *   <root>\tmp_2                missing                      -> skipped silently, as ENOENT on POSIX
 *   <root>\tmp_3\g.txt          a second folder              -> removed
 * Then a prefix containing '<' (an invalid Windows name): GetFileAttributesA fails with an error other
 * than "not found", which must be reported, not skipped. The workflow checks stdout for that warning.
 *
 * Built by `make test_remove_tmp_win` (MSYS2 UCRT64 gcc).
 */
#include "include/define.h"
#include "include/data_sparse.h"
#include <stdio.h>
#include <string.h>
// As in data_sparse.c: full windows.h pulls in winioctl.h, whose MEDIA_TYPE enumerator `Unknown`
// collides with the one in struct.h.
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

static int failures = 0;
static void expect(int ok, const char *what) {
  printf("%s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) failures++;
}
static int exists(const char *p) { return GetFileAttributesA(p) != INVALID_FILE_ATTRIBUTES; }
static void write_file(const char *p) { FILE *f = fopen(p, "w"); if (f) { fputs("x\n", f); fclose(f); } }

int main(int argc, char **argv) {
  if (argc != 2) { printf("usage: %s <new test root>\n", argv[0]); return 2; }
  const char *root = argv[1];
  char p[1024], q[1024];
  if (exists(root)) { printf("test root %s already exists: refusing to use it\n", root); return 2; }
  if (!CreateDirectoryA(root, NULL)) { printf("cannot create %s\n", root); return 2; }
  // Proof for the workflow that THIS invocation created the root, so only then may it clean up.
  snprintf(p, sizeof p, "%s\\.created-by-test_remove_temporary_folders_win", root); write_file(p);

  snprintf(p, sizeof p, "%s\\keep", root); CreateDirectoryA(p, NULL);
  snprintf(p, sizeof p, "%s\\keep\\sentinel.txt", root); write_file(p);
  snprintf(p, sizeof p, "%s\\tmp_1", root); CreateDirectoryA(p, NULL);
  snprintf(p, sizeof p, "%s\\tmp_1\\a", root); CreateDirectoryA(p, NULL);
  snprintf(p, sizeof p, "%s\\tmp_1\\a\\b", root); CreateDirectoryA(p, NULL);
  snprintf(p, sizeof p, "%s\\tmp_1\\a\\b\\f.txt", root); write_file(p);
  snprintf(p, sizeof p, "%s\\tmp_3", root); CreateDirectoryA(p, NULL);
  snprintf(p, sizeof p, "%s\\tmp_3\\g.txt", root); write_file(p);
  // Junction tmp_1\j -> keep (a directory junction needs no special privilege).
  snprintf(q, sizeof q, "cmd /c mklink /J \"%s\\tmp_1\\j\" \"%s\\keep\" >NUL", root, root);
  int rc = system(q);
  snprintf(p, sizeof p, "%s\\tmp_1\\j", root);
  DWORD ja = GetFileAttributesA(p);
  expect(rc == 0 && ja != INVALID_FILE_ATTRIBUTES && (ja & FILE_ATTRIBUTE_REPARSE_POINT), "fixture: tmp_1\\j is a junction");
  snprintf(p, sizeof p, "%s\\tmp_1\\j\\sentinel.txt", root);
  expect(exists(p), "fixture: the sentinel is reachable through the junction");

  snprintf(p, sizeof p, "%s\\tmp_", root);
  Remove_temporary_folders(p, 3);

  snprintf(p, sizeof p, "%s\\tmp_1", root); expect(!exists(p), "tmp_1 and its nested contents removed");
  snprintf(p, sizeof p, "%s\\tmp_3", root); expect(!exists(p), "tmp_3 removed");
  snprintf(p, sizeof p, "%s\\keep", root); expect(exists(p), "junction target keep\\ still exists");
  snprintf(p, sizeof p, "%s\\keep\\sentinel.txt", root); expect(exists(p), "keep\\sentinel.txt still exists (junction not followed)");

  // Non-missing error: '<' is invalid in a Windows name. Remove_temporary_folders must print a warning.
  printf("--- invalid-name case, expect a 'unable to remove temporary folder' warning below:\n");
  fflush(stdout);
  snprintf(p, sizeof p, "%s\\bad<", root);
  Remove_temporary_folders(p, 1);
  fflush(stdout);

  printf(failures ? "FAILED: %d check(s)\n" : "PASS: Windows temporary-folder removal\n", failures);
  return failures ? 1 : 0;
}
