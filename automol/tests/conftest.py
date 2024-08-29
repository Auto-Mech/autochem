"""Pytest fixtures."""

from pathlib import Path

import pytest


@pytest.fixture
def data_directory_path():
    return Path(__file__).parent / "data"
