from pgap2.postprocess.stat import stat_cmd
from pgap2.postprocess.profile import profile_cmd
from pgap2.postprocess.baps import baps_cmd
from pgap2.postprocess.singletree import singletree_cmd
from pgap2.postprocess.tajimas_d import tajimas_d_cmd
from argparse import ArgumentParser, _SubParsersAction
import argparse


def postprocess_cmd(subparser: _SubParsersAction):
    main_subparser_postprocess: ArgumentParser = subparser.add_parser(
        'post', help='Postprocess the output files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub_subparsers = main_subparser_postprocess.add_subparsers()
    stat_cmd(sub_subparsers)
    profile_cmd(sub_subparsers)
    singletree_cmd(sub_subparsers)
    baps_cmd(sub_subparsers)
    tajimas_d_cmd(sub_subparsers)
