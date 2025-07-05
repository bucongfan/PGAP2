#!/usr/bin/env python
# coding:utf-8
from time import time
from loguru import logger
from argparse import ArgumentParser, _SubParsersAction, Namespace, ArgumentDefaultsHelpFormatter

from pgap2.utils.partition import partition_cmd
from pgap2.utils.preprocess import preprocess_cmd
from pgap2.utils.postprocess import postprocess_cmd

from pgap2 import __version__, __author__


def main():
    starttime = time()
    parser = ArgumentParser(description="Pan-Genome Analysis Pipeline",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"%(prog)s {__version__}  (author: {__author__})"
    )
    subparser: _SubParsersAction = parser.add_subparsers()

    preprocess_cmd(subparser)
    partition_cmd(subparser)
    postprocess_cmd(subparser)
    args: Namespace = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        exit(1)
    endtime = time()
    logger.success("Total time used: {:.2f}s".format(
        endtime - starttime))


if __name__ == "__main__":
    main()
