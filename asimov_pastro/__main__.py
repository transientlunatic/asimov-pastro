"""
Entry point for asimov-pastro CLI.

This avoids circular imports by not importing the main package.
"""

if __name__ == '__main__':
    from asimov_pastro.cli import cli
    cli()
