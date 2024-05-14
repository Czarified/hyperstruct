"""Command-line interface."""

from importlib.metadata import version

import click


__version__ = version("hyperstruct")


@click.command()
@click.version_option()
def main() -> None:
    """Hyperstruct."""


if __name__ == "__main__":
    main(prog_name="hyperstruct")  # pragma: no cover
