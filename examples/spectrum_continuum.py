"""External continuum correction on three unequal chips; no FERRE required."""

import numpy as np
from aspcap import Spec1D, continuum

lengths = (120, 100, 90)
wave = [15000 + 200 * i + np.arange(n) for i, n in enumerate(lengths)]
model_flux = [1 - 0.2 * np.exp(-(((np.arange(n) - n / 2) / 3) ** 2)) for n in lengths]
reference = Spec1D(model_flux, wave=wave)
observed = Spec1D(
    [f * (2 + i) for i, f in enumerate(model_flux)],
    wave=wave,
    err=[np.full(n, 0.02) for n in lengths],
)
estimate = continuum.fit(observed, method="model_ratio", model=reference, window=11)
normalized = estimate.apply(observed)
for i, n in enumerate(lengths):
    np.testing.assert_allclose(normalized.flux[:n, i], model_flux[i])
print(normalized)
print("Per-chip scales:", [estimate.continuum[0, i] for i in range(3)])
print("Original fluxes preserved:", not observed.normalized)
