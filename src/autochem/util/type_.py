"""Type utilities."""

import abc
from collections.abc import Callable, Mapping
from typing import Annotated, ClassVar, Self, TypeAlias, TypeVar

import numpy as np
import pint
import pydantic
from numpy.typing import ArrayLike, NDArray
from pydantic_core import core_schema

from . import form
from .form import Formula

# Type variables
T = TypeVar("T")
S = TypeVar("S", bound="SubclassTyped")


# Type aliases
Number: TypeAlias = float | int
Scalers: TypeAlias = dict[str, Callable[[ArrayLike, ArrayLike], NDArray[np.float64]]]


# Abstract base classes
class Frozen(pydantic.BaseModel, abc.ABC):
    """Abstract base class for frozen models.

    Enforces faux-immutability to prevent model fields from changing, preventing data
    from being corrupted by inconsistent mutations.
    """

    model_config = pydantic.ConfigDict(frozen=True)


class SubclassTyped(pydantic.BaseModel, abc.ABC):
    """Abtract base class for models with typed subclasses.

    Allows subclasses to be tagged with a "type" field to distinguish them from each
    other.

    Usage:
    ```
    from typing import ClassVar

    class A(SubclassTyped):
        x: int

    class B(A):
        type: ClassVar[str] = "b"

    class C(A):
        type: ClassVar[str] = "c"

    A.from_data({"x": 1, "type": "b"})

    # output: B(x=5)
    ```
    """

    type_: ClassVar[str] = ""

    @pydantic.computed_field
    def type(self) -> str:
        """Subclass type."""
        return self.__class__.type_

    @classmethod
    def model_validate(
        cls,
        obj: object,
        *,
        strict: bool | None = None,
        from_attributes: bool | None = None,
        context: object | None = None,
    ) -> Self:
        """Validate a pydantic model instance."""
        if isinstance(obj, Mapping) and "type" in obj:
            assert "type" in obj, f"Missing type field: {obj}"
            obj = dict(obj).copy()
            typ = obj.pop("type")
            sub = next((s for s in subclasses(cls) if s.type_ == typ), None)
            if sub is not None:
                return sub.model_validate(
                    obj,
                    strict=strict,
                    from_attributes=from_attributes,
                    context=context,
                )

        return super().model_validate(
            obj,
            strict=strict,
            from_attributes=from_attributes,
            context=context,
        )


def subclasses(cls: type[SubclassTyped]) -> list[type[SubclassTyped]]:
    """Determine subclasses of a class."""
    subs = []

    for sub in cls.__subclasses__():
        subs.append(sub)
        subs.extend(subclasses(sub))

    return subs


class Scalable(pydantic.BaseModel, abc.ABC):
    """Abstract base class for models with scalar multiplication.

    :param _scalers: A dictionary of functions defining how the function parameters
        change upon scalar multiplication, by attribute name
    """

    _scalers: ClassVar[Scalers | None] = None

    def __mul__(self, c: ArrayLike) -> Self:
        """Scalar multiplication.

        :param c: Scalar value to multiply by
        :return: Scaled object
        """
        data = self.model_dump()

        if self._scalers is not None:
            for key, scaler in self._scalers.items():
                data[key] = scaler(c, getattr(self, key))
        else:
            msg = "Scalar multiplication not implemented."
            raise NotImplementedError(msg)

        return self.model_validate(data)

    def __truediv__(self, c: ArrayLike) -> Self:
        """Scalar division.

        :param c: Scalar value to divide by
        :return: Scaled object
        """
        return self.__mul__(1 / c)

    __rmul__ = __mul__


# Annotated types for pydantic
Unit_ = Annotated[
    pydantic.SkipValidation[pint.Unit],
    pydantic.BeforeValidator(lambda x: pint.Unit(x)),
    # Use abbreviated unit names upon serialization
    pydantic.PlainSerializer(lambda x: format(x, "~")),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(str)),
    ),
]

Quantity_ = Annotated[
    pydantic.SkipValidation[pint.Quantity],
    pydantic.BeforeValidator(lambda x: pint.Quantity(x)),
    # Use abbreviated unit names upon serialization
    pydantic.PlainSerializer(lambda x: format(x, "~")),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(str)),
    ),
]

NDArray_ = Annotated[
    pydantic.SkipValidation[np.ndarray],
    pydantic.BeforeValidator(lambda x: np.array(x, dtype=np.float64)),
    pydantic.PlainSerializer(lambda x: np.array(x).tolist()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(list[float])),
    ),
]

Formula_ = Annotated[
    Formula,
    pydantic.BeforeValidator(form.normalize_input),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(Formula)),
    ),
]
