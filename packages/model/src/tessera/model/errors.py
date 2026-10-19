"""Error types shared across packages."""


class InputError(ValueError):
    """A problem in what the user asked for (unknown part, missing material, bad face...).

    Clients show its message as-is; anything else escaping ``tessera.core`` is a bug.
    """
