#!/usr/bin/env python
# coding:utf-8
from time import time
from loguru import logger
from argparse import ArgumentParser, _SubParsersAction, Namespace, ArgumentDefaultsHelpFormatter

from pgap2.utils.partition import partition_cmd
from pgap2.utils.preprocess import preprocess_cmd
from pgap2.utils.postprocess import postprocess_cmd


def main():
    starttime = time()
    parser = ArgumentParser(description="Pan-Genome Analysis Pipeline",
                            formatter_class=ArgumentDefaultsHelpFormatter)
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
