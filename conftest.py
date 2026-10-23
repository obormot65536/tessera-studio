"""Repository-wide pytest hooks."""

import os

import pytest


def pytest_runtest_logreport(report: pytest.TestReport) -> None:
    """On GitHub Actions, turn each failure into an annotation on the check run.

    Annotations are readable through the checks API even where the job log is not.
    """
    if not os.environ.get("GITHUB_ACTIONS") or not report.failed:
        return
    path, line, _ = report.location
    text = str(report.longrepr)[-3000:]  # the tail holds the assertion and the error
    message = _escape(text)
    title = _escape(report.nodeid).replace(":", "%3A").replace(",", "%2C")
    print(f"\n::error file={path},line={(line or 0) + 1},title={title}::{message}")


def _escape(text: str) -> str:
    # ASCII only: the Windows runner's console encoding cannot print every character.
    text = text.encode("ascii", "backslashreplace").decode("ascii")
    return text.replace("%", "%25").replace("\r", "").replace("\n", "%0A")
