"""Type utilities."""

import abc
from collections.abc import Callable, Mapping
from typing import Annotated, Any, ClassVar, Self, TypeAlias, TypeVar

import numpy
import pint
import pydantic
from numpy.typing import ArrayLike, NDArray
from pydantic_core import core_schema

# Type variables
T = TypeVar("T")
S = TypeVar("S", bound="SubclassTyped")


# Type aliases
Number: TypeAlias = float | int
Scalers: TypeAlias = dict[str, Callable[[ArrayLike, ArrayLike], NDArray[numpy.float64]]]


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
    from typing import Final

    class A(SubclassTyped):
        x: int

    class B(A):
        _type: str = "b"

    class C(A):
        _type: str = "c"

    A.from_data({"x": 1, "type": "b"})

    # output: B(x=5)
    ```
    """

    _type: str

    @pydantic.computed_field
    def type(self) -> str:
        """Subclass type."""
        return self._type

    @classmethod
    def model_validate(
        cls,
        obj: Any,
        *,
        strict: bool | None = None,
        from_attributes: bool | None = None,
        context: Any | None = None,
    ) -> Self:
        """Validate a pydantic model instance."""
        if isinstance(obj, Mapping) and "type" in obj:
            assert "type" in obj, f"Missing type field: {obj}"
            obj = dict(obj).copy()
            _type = obj.pop("type")
            sub = next((s for s in cls.__subclasses__() if s._type == _type), None)
            if sub is not None:
                return sub.model_validate(
                    obj, strict=strict, from_attributes=from_attributes, context=context
                )

        return super().model_validate(
            obj, strict=strict, from_attributes=from_attributes, context=context
        )


class Scalable(pydantic.BaseModel, abc.ABC):
    """Abstract base class for models with scalar multiplication.

    :param _scalers: A dictionary of functions defining how the function parameters
        change upon scalar multiplication, by attribute name
    """

    _scalers: ClassVar[Scalers | None] = None

    def __mul__(self, c: ArrayLike):
        """Scalar multiplication.

        :param c: Scalar value to multiply by
        :return: Scaled function
        """
        data = self.model_dump()

        if self._scalers is not None:
            for key, scaler in self._scalers.items():
                data[key] = scaler(c, getattr(self, key))

        return self.model_validate(data)

    __rmul__ = __mul__


# Annotated types for pydantic
Unit_ = Annotated[
    pydantic.SkipValidation[pint.Unit],
    pydantic.BeforeValidator(lambda x: pint.Unit(x)),
    # Use abbreviated unit names upon serialization
    pydantic.PlainSerializer(lambda x: format(x, "~")),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(str))
    ),
]

Quantity_ = Annotated[
    pydantic.SkipValidation[pint.Quantity],
    pydantic.BeforeValidator(lambda x: pint.Quantity(x)),
    # Use abbreviated unit names upon serialization
    pydantic.PlainSerializer(lambda x: format(x, "~")),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(str))
    ),
]

NDArray_ = Annotated[
    pydantic.SkipValidation[numpy.ndarray],
    pydantic.BeforeValidator(lambda x: numpy.array(x, dtype=numpy.float64)),
    pydantic.PlainSerializer(lambda x: numpy.array(x).tolist()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(list[float]))
    ),
]

# # Annotated type for xarray
# DataArray_ = Annotated[
#     pydantic.SkipValidation[xarray.DataArray],
#     pydantic.BeforeValidator(lambda x: xarray.DataArray(x)),
#     pydantic.PlainSerializer(lambda x: xarray.DataArray(x).to_dict()),
#     pydantic.GetPydanticSchema(
#         lambda _, handler: core_schema.with_default_schema(handler(dict[str, object]))
#     ),
# ]
