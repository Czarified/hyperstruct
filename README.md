# Hyperstruct

[![PyPI](https://img.shields.io/pypi/v/hyperstruct.svg)][pypi_]
[![Status](https://img.shields.io/pypi/status/hyperstruct.svg)][status]
[![Python Version](https://img.shields.io/pypi/pyversions/hyperstruct)][python version]
[![License](https://img.shields.io/pypi/l/hyperstruct)][license]

[![Read the documentation at https://hyperstruct.readthedocs.io/](https://img.shields.io/readthedocs/hyperstruct/latest.svg?label=Read%20the%20Docs)][read the docs]
[![Tests](https://github.com/czarified/hyperstruct/workflows/Tests/badge.svg)][tests]
[![Codecov](https://codecov.io/gh/czarified/hyperstruct/branch/main/graph/badge.svg)][codecov]

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)][pre-commit]
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)][black]

[pypi_]: https://pypi.org/project/hyperstruct/
[status]: https://pypi.org/project/hyperstruct/
[python version]: https://pypi.org/project/hyperstruct
[read the docs]: https://hyperstruct.readthedocs.io/
[tests]: https://github.com/czarified/hyperstruct/actions?workflow=Tests
[codecov]: https://app.codecov.io/gh/czarified/hyperstruct
[pre-commit]: https://github.com/pre-commit/pre-commit
[black]: https://github.com/psf/black

HyperStruct is a modern implementation of the SWEEP software originally
produced by Rockwell International in 1974. The Structural Weight Estimation
Program for Aircraft (SWEEP) was originally implemented in FORTRAN. The
source code for this program is no longer publicly available, but
the US Defence Technical Information Center has published all the original
documentation from Rockwell (which is extensive).

This project aims to meet the same objectives as SWEEP: devlop analytical
computer programs capable of predicting the weight of aircraft structural
components suitable for use during the conceptual phase of the aircraft
design cycle. It builds a very high level of fidelity in the weight
estimates, with flexiblity to a plethora of design variables. All of this
is predominantly derived from fundamental structural methods, or historically
predicated empirical formulas.

## Features

- TODO

## Requirements

- TODO

## Installation

You can install _Hyperstruct_ via [pip] from [PyPI]:

```console
$ pip install hyperstruct
```

## Usage

Please see the [Command-line Reference] for details.

## Contributing

Contributions are very welcome.
To learn more, see the [Contributor Guide].

## License

Distributed under the terms of the [MIT license][license],
_Hyperstruct_ is free and open source software.

## Issues

If you encounter any problems,
please [file an issue] along with a detailed description.

## Credits

This project was generated from [@cjolowicz]'s [Hypermodern Python Cookiecutter] template.

[@cjolowicz]: https://github.com/cjolowicz
[pypi]: https://pypi.org/
[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
[file an issue]: https://github.com/czarified/hyperstruct/issues
[pip]: https://pip.pypa.io/

<!-- github-only -->

[license]: https://github.com/czarified/hyperstruct/blob/main/LICENSE
[contributor guide]: https://github.com/czarified/hyperstruct/blob/main/CONTRIBUTING.md
[command-line reference]: https://hyperstruct.readthedocs.io/en/latest/usage.html
