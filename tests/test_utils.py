import subprocess

from fixtures import *
from utils import get_value_from_line, run_cmd


def test_get_value_from_line():
    line = "Test number: 100.0"
    value = get_value_from_line(line, "number")

    assert isinstance(value, float)
    assert value == pytest.approx(100)


def test_get_value_from_line_raises_value_error_if_string_not_in_line():
    line = "Test number: 100.0"

    with pytest.raises(ValueError):
        get_value_from_line(line, "bananas")


def test_run_cmd():
    run_cmd(["ls"])


def test_run_cmd_raises_exception():
    with pytest.raises(Exception):
        run_cmd(["thisIsNotAValidCommand"])
