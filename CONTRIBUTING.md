# Contributor Guide

Thank you for your interest in improving this project.
This project is open-source under the [MIT license] and
welcomes contributions in the form of bug reports, feature requests, and pull requests.

Here is a list of important resources for contributors:

- [Source Code]
- [Documentation]
- [Issue Tracker]
- [Code of Conduct]

[mit license]: https://opensource.org/licenses/MIT
[source code]: https://github.com/czarified/hyperstruct
[documentation]: https://hyperstruct.readthedocs.io/
[issue tracker]: https://github.com/czarified/hyperstruct/issues

## How to report a bug

Report bugs on the [Issue Tracker].

When filing an issue, make sure to answer these questions:

- Which operating system and Python version are you using?
- Which version of this project are you using?
- What did you do?
- What did you expect to see?
- What did you see instead?

The best way to get your bug fixed is to provide a test case,
and/or steps to reproduce the issue.

## How to request a feature

Request features on the [Issue Tracker].

## How to set up your development environment

You need Python 3.10+ and the following tools:

- [Poetry]
- [poetry-export-plugin]
- [Nox]
- [nox-poetry]

The general steps are:

1.  Make sure you have the above tools installed.
2.  Clone the repository locally.
3.  Navigate to the cloned directory in your terminal.
4.  Install the package with development requirements:

    ```console
    $ poetry install --with dev
    ```

5.  Verify you have the latest version

You can now run an interactive Python session,
or the command-line interface:

```console
$ poetry run python
$ poetry run hyperstruct
```

[poetry]: https://python-poetry.org/
[poetry-export-plugin]: https://github.com/python-poetry/poetry-plugin-export
[nox]: https://nox.thea.codes/
[nox-poetry]: https://nox-poetry.readthedocs.io/

## How to test the project

Run the partial test suite:

```console
$ poetry run pytest
```

Run the full test suite:

```console
$ poetry run nox
```

List the available Nox sessions:

```console
$ poetry run nox --list-sessions
```

You can also run a specific Nox session.
For example, invoke the unit test suite like this:

```console
$ poetry run nox --session=tests
```

Running the full nox suite by default will rebuild all virtual environments. This is pretty slow, and most of the time unnecessary.
So, if you want to reuse the previously-built environments, you can add the -r flag:

```console
$ poetry run nox -r
```

Unit tests are located in the _tests_ directory,
and are written using the [pytest] testing framework.

[pytest]: https://pytest.readthedocs.io/

## How to submit changes

Open a [pull request] to submit changes to this project.

Your pull request needs to meet the following guidelines for acceptance:

- The Nox test suite must pass without errors and warnings.
- Include unit tests. This project aims for 100% code coverage.
- If your changes add functionality, update the documentation accordingly.

Feel free to submit early, though we can always iterate on this.

To run linting and code formatting checks before committing your change, you can install pre-commit as a Git hook by running the following command:

```console
$ poetry run nox --session=pre-commit -- install
```

It is recommended to open an issue before starting work on anything.
This will allow a chance to talk it over with the owners and validate your approach.

[pull request]: https://github.com/czarified/hyperstruct/pulls

<!-- github-only -->

[code of conduct]: CODE_OF_CONDUCT.md
