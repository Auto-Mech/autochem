"""Library of automol errors."""


class FailedGeometryGenerationError(RuntimeError):
    """exception for when we fail to generate a correct geometry."""


class FailedInchiGenerationError(RuntimeError):
    """exception for when we fail to generate a correct inchi."""
