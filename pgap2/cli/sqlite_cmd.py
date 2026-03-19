"""CLI definitions for the ``pgap2 sqlite`` subcommand.

Generates a self-contained SQLite database from PGAP2 output files.

Usage::

    pgap2 sqlite -i /path/to/output
    pgap2 sqlite -i /path/to/output -o /path/to/output/pgap2.db
    pgap2 sqlite -i /path/to/output --name "My Project" --species "E. coli"
"""

import os
import argparse
from argparse import ArgumentParser, _SubParsersAction
from pathlib import Path


def launch(args: argparse.Namespace):
    from pgap2.utils.sqlite_builder import build_sqlite, parse_readme

    indir = os.path.abspath(args.indir)

    # Resolve output path: default is <indir>/pgap2.db
    if args.output:
        output_path = os.path.abspath(args.output)
    else:
        output_path = os.path.join(indir, "pgap2.db")

    # Try to derive project name / species from README.md if not provided
    project_name = args.name
    species = args.species

    readme_path = Path(indir) / "README.md"
    if readme_path.exists():
        readme_info = parse_readme(readme_path)
        if not project_name:
            project_name = readme_info.get("project_name", "")
        if not species:
            species = readme_info.get("species", "")

    if not project_name:
        project_name = Path(indir).name
    if not species:
        species = "Unknown"

    info = build_sqlite(
        data_dir=indir,
        output_path=output_path,
        project_name=project_name,
        species=species,
        description=args.description or "",
    )

    return 0


def sqlite_portal(args):
    """Entry point called by argparse (dependency checks + launch)."""
    launch(args)


def sqlite_cmd(subparser: _SubParsersAction):
    """Register the ``sqlite`` subcommand."""
    sp: ArgumentParser = subparser.add_parser(
        "sqlite",
        help="Generate a SQLite database from PGAP2 output files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sp.add_argument(
        "--indir", "-i", required=True,
        help="PGAP2 output directory (contains .pav, .detail.tsv, etc.)")
    sp.add_argument(
        "--output", "-o", required=False, default=None,
        help="Output SQLite file path (default: <indir>/pgap2.db)")
    sp.add_argument(
        "--name", "-n", required=False, default=None,
        help="Project name (auto-detected from README.md if present)")
    sp.add_argument(
        "--species", "-s", required=False, default=None,
        help="Species name (auto-detected from README.md if present)")
    sp.add_argument(
        "--description", "-d", required=False, default="",
        help="Optional project description")
    sp.set_defaults(func=sqlite_portal)
