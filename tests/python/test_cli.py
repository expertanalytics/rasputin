"""The console script declared in pyproject must actually resolve and run."""

from typer.testing import CliRunner

from tin_engine.cli import app

runner = CliRunner()


def test_entry_point_is_importable() -> None:
    # Regression: [project.scripts] declared tin_engine.cli:app while the module
    # did not exist, so `pip install .` succeeded and `rasputin` then died with
    # ModuleNotFoundError.
    assert app is not None


def test_version_reports_the_installed_distribution() -> None:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip()


def test_bare_invocation_shows_help_rather_than_failing() -> None:
    result = runner.invoke(app, [])
    # Click's no_args_is_help convention: print usage, exit 2. The point of
    # this test is that a bare `rasputin` shows help instead of a traceback,
    # not that it succeeds -- no arguments is a usage error.
    assert result.exit_code == 2
    assert "Usage" in result.stdout
    assert result.exception is None or isinstance(result.exception, SystemExit)
