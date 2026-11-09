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


def _run_batch(
    spectra, function, *, star_ids=None, parameters=None, raise_on_error=False
):
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
            if raise_on_error:
                raise
            item.error = str(error)
        result.append(item)
    return result


def batch_spectra(function):
    """Dispatch list/tuple spectra through an otherwise single-spectrum method."""
    from collections.abc import Mapping
    from functools import wraps
    from inspect import signature

    sig = signature(function)

    @wraps(function)
    def dispatch(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)
        spectrum = bound.arguments["spectrum"]
        star_ids = bound.arguments.pop("star_ids", None)
        if not isinstance(spectrum, (list, tuple)):
            if star_ids is not None:
                raise ValueError("star_ids applies only to a list of spectra")
            return function(*bound.args, **bound.kwargs)
        count = len(spectrum)
        per_input = {}
        for name in ("initial", "parameters"):
            value = bound.arguments.get(name)
            if (
                isinstance(value, (list, tuple))
                and value
                and isinstance(value[0], Mapping)
            ):
                if len(value) != count or not all(
                    isinstance(v, Mapping) for v in value
                ):
                    raise ValueError(
                        f"{name} must have one parameter mapping per spectrum"
                    )
                per_input[name] = value
        if star_ids is None and bound.arguments.get("star_id") is not None:
            star_ids = [bound.arguments["star_id"]] * count

        # Use an explicit index parameter, including for invalid inputs, so one
        # rejected object cannot shift subsequent per-spectrum starting values.
        def indexed(spec, identity, index):
            options = dict(bound.arguments, spectrum=spec)
            for name, values in per_input.items():
                options[name] = values[index]
            if "star_id" in sig.parameters:
                options["star_id"] = identity
            return function(**options)

        return _run_batch(
            spectrum,
            indexed,
            star_ids=star_ids,
            parameters=range(count),
            raise_on_error=bound.arguments.get("raise_on_error", False),
        )

    dispatch.__doc__ = (function.__doc__ or "") + """

        A list or tuple of Spec1D inputs returns an ordered list of BatchItem
        objects; a single Spec1D retains the normal return type. Batch failures
        are isolated unless raise_on_error=True is available and requested.
        initial/parameters may be one shared mapping or one mapping per input.
        star_ids optionally supplies one identity per batch input. Other options
        are shared across the batch. Processing is serial; inputs are copied by
        the same rules as single-spectrum calls.
        """
    return dispatch
