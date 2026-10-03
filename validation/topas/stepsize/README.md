# TOPAS step-size and nuclear-elastic diagnostics

Diagnostic runs for §5.4 (Table 6b). Not pre-specified and not part of the frozen halo addendum. Table 6 does not use them.

All runs use the TOPAS setup in `validation/topas/` at commit `6c4aaff81994`. Its `validation/topas` is the same as the halo-addendum freeze `5aedf3ac49d7`. Endpoints are from `validation/pencil_endpoints.py`.

Each prediction was committed before its runs, in clement-7074f29f's office repository:

| test | prediction | results |
|---|---|---|
| 200 MeV: 0.5 and 0.1 mm step limits | `1b1fbe5` 2026-10-03 20:28 | `542bd57` 20:49 |
| 100 and 150 MeV: 0.1 mm; 200 MeV: 0.05 mm | `946bc6e` 20:58 | `e41298b` 21:35 |
| 100 and 150 MeV: 0.05 mm | `3a21014` 2026-10-04 01:31 | `0336aaa` 01:46 |
| 200 MeV: `g4h-elastic_HP` removed | `10c261b` 01:36 | `a5a083d` 02:22 |

## Files

- `PREDICTION*.md`, `RESULTS*.md`: the prediction and result of each test.
- `run_*.sh`: the drivers. They write into the TOPAS output store on the Mac Studio (`/Volumes/T7 Shield/SMBWritable/topas/stepsize/6c4aaff81994/`). The outputs are read-only (mode 444) and are not in this repository.
- `endpoints.sh`, `followup2_endpoints.jsonl`: endpoint records for follow-up 2.
- `elastic_*_endpoints.jsonl`, `compare_elastic.py`, `elastic_compare.txt`: the elastic-off test (control = the 8 halo-addendum 200 MeV runs).
- `doseall_rings.sh`, `doseall_rings.jsonl`: `Dose` and `DoseAll` per ring, used for the scorer-filter check.
- `*.log`: run logs (start time, seed, exit code, wall time).
