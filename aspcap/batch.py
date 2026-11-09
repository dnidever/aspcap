"""Ordered batch results with independent per-spectrum failures."""

from dataclasses import dataclass


@dataclass
class BatchItem:
    """One input's identity, returned stage value, or exception message."""

    index: int
    star_id: str
    value: object = None
    error: str = None

    @property
    def success(self):
        if self.error is not None:
            return False
        if isinstance(self.value, dict):
            return all(
                bool(getattr(value, "success", True)) for value in self.value.values()
            )
        return bool(getattr(self.value, "success", True))


def run_many(spectra, function, *, star_ids=None, parameters=None):
    """Apply a stage serially, retaining input order and isolating failures.

    parameters is an optional sequence of one complete parameter mapping per
    input. This explicit sequence avoids ambiguity between shared and per-star
    starts. The callback receives (spectrum, star_id, parameters_or_None).
    """
    spectra = list(spectra)
    if star_ids is None:
        star_ids = [
            getattr(s, "metadata", {}).get("star_id")
            or (getattr(s, "filename", None) or f"star-{i}")
            for i, s in enumerate(spectra)
        ]
    else:
        star_ids = list(star_ids)
    parameters = [None] * len(spectra) if parameters is None else list(parameters)
    if len(star_ids) != len(spectra) or len(parameters) != len(spectra):
        raise ValueError("star_ids and parameters must match the number of spectra")
    result = []
    for index, (spectrum, star_id, parameter) in enumerate(
        zip(spectra, star_ids, parameters)
    ):
        item = BatchItem(index, str(star_id))
        try:
            from .spectrum import Spec1D

            if not isinstance(spectrum, Spec1D):
                raise TypeError("Batch inputs must be Spec1D objects")
            item.value = function(spectrum, item.star_id, parameter)
        except (ValueError, TypeError, RuntimeError, OSError, KeyError) as error:
            item.error = str(error)
        result.append(item)
    return result
