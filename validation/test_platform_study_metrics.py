"""Synthetic tests for the R80/R20 crossing rule, run before any study data (design, "Distal edge").

Run: python validation/test_platform_study_metrics.py   (exit 0 = all pass)
"""

import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from platform_study_metrics import distal_level  # noqa: E402

z = (np.arange(150) + 0.5) * 2.0  # voxel-centre depths 1, 3, ..., 299 mm
failures = 0


def check(ok: bool, what: str) -> None:
    global failures
    print(("ok   " if ok else "FAIL ") + what)
    failures += 0 if ok else 1


# 1. Piecewise linear: flat 1.0 to 100 mm, then linearly to 0 at 200 mm. Interpolation is exact there.
lin = np.clip(1 - (z - 100) / 100, 0, 1)
r80, n80 = distal_level(z, lin, 0.8)
r20, n20 = distal_level(z, lin, 0.2)
check(r80 is not None and abs(r80 - 120) < 1e-9 and n80 == 1, f"linear: R80 = 120 mm exactly, one crossing (got {r80}, {n80})")
check(r20 is not None and abs(r20 - 180) < 1e-9 and n20 == 1, f"linear: R20 = 180 mm exactly, one crossing (got {r20}, {n20})")

# 2. Smooth fall-off 1 - Phi((z - mu)/sigma) after a rising build-up: R80 = mu - 0.8416 sigma, R20 = mu + 0.8416 sigma.
mu, sigma = 240.0, 6.0
phi = 0.5 * (1 + np.vectorize(math.erf)((z - mu) / (sigma * math.sqrt(2))))
smooth = (0.5 + 0.5 * np.minimum(z / 200, 1)) * (1 - phi)
zq = 0.8416212335729143
peak = float(np.max(smooth))
# The known depths are those where smooth equals frac * its own max; solve on a fine grid rather than assume.
zf = np.linspace(200, 299, 200001)
sf = (0.5 + 0.5 * np.minimum(zf / 200, 1)) * (1 - 0.5 * (1 + np.vectorize(math.erf)((zf - mu) / (sigma * math.sqrt(2)))))
true80 = float(zf[np.argmax(sf < 0.8 * peak)])
true20 = float(zf[np.argmax(sf < 0.2 * peak)])
r80, n80 = distal_level(z, smooth, 0.8)
r20, n20 = distal_level(z, smooth, 0.2)
check(abs(true80 - (mu - zq * sigma)) < 0.05 and abs(true20 - (mu + zq * sigma)) < 0.05, "smooth: fine-grid truth matches mu -/+ 0.8416 sigma")
check(r80 is not None and abs(r80 - true80) < 0.1 and n80 == 1, f"smooth: R80 within 0.1 mm of {true80:.3f} (got {r80}, {n80})")
check(r20 is not None and abs(r20 - true20) < 0.1 and n20 == 1, f"smooth: R20 within 0.1 mm of {true20:.3f} (got {r20}, {n20})")

# 3. Re-crossing: a bump lifts the curve back above 80% after the first crossing. First crossing kept, count 2.
bumped = lin.copy()
i = int(np.argmax(z > 130))
bumped[i] = 0.9
r80b, n80b = distal_level(z, bumped, 0.8)
check(r80b is not None and abs(r80b - 120) < 1e-9 and n80b == 2, f"re-crossing: first R80 (120) kept and flagged, 2 crossings (got {r80b}, {n80b})")

# 4. Absent: the curve never falls below 30%, so R20 is undefined; R80 still found.
floor = np.maximum(lin, 0.3)
r20a, n20a = distal_level(z, floor, 0.2)
r80a, _ = distal_level(z, floor, 0.8)
check(r20a is None and n20a == 0, f"absent: R20 undefined (got {r20a}, {n20a})")
check(r80a is not None and abs(r80a - 120) < 1e-9, f"absent: R80 still 120 (got {r80a})")

# 5. A proximal fall must not count: crossings are searched only distal to the maximum. The curve starts at 0.95,
#    dips to 0.5 around 30 mm (a downward crossing of 80% BEFORE the peak), then rises to the 1.0 plateau.
rise = np.where(z < 20, 0.95, np.where(z < 40, 0.5, np.where(z < 60, 0.5 + (z - 40) / 40, lin)))
rise = np.minimum(rise, np.where(z < 20, 0.95, 1.0))
r80r, n80r = distal_level(z, rise, 0.8)
check(r80r is not None and abs(r80r - 120) < 1e-9 and n80r == 1, f"proximal rise ignored (got {r80r}, {n80r})")

print("PASS" if not failures else f"FAILED: {failures}")
sys.exit(1 if failures else 0)
