import os
import shutil
import subprocess

from tqdm import tqdm
from loguru import logger
from datetime import datetime


class _tqdm_prefix:
    def __init__(self):
        self.total_steps = None

    def set_total_step(self, total: int):
        self.total_steps = total

    def step(self, step):
        if step == -1:
            step = f'step- {self.total_steps}/{self.total_steps}'
        else:
            step = f'step- {step}/{self.total_steps}'
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_prefix = f"{current_time} | INFO | {step} "
        return log_prefix


tqdm_ = _tqdm_prefix()


class _Externals():
    def __init__(self) -> None:
        self.externals = {
            'cdhit': 'cd-hit',
            'mmseqs2': 'mmseqs',
            'muscle': 'muscle',
            'mafft': 'mafft',
            'tcoffee': 't_coffee',
            'raxml': 'raxml-ng',
            'iqtree': 'iqtree',
            'fasttree': 'fasttree',
            'miniprot': 'miniprot',
            'clipkit': 'clipkit',
            'fastbaps': 'run_fastbaps',
            'cfml': 'ClonalFrameML',
            'maskrc': 'maskrc-svg.py',
            'prodigal': 'prodigal',
            'diamond': 'diamond',
            'blastp': 'blastp',
            'makeblastdb': 'makeblastdb',
            'seqtk': 'seqtk',
            'mcxdeblast': 'mcxdeblast',
            'mcl': 'mcl',
            'draw_post_stat': 'draw_post_stat.r',
            'draw_post_profile': 'draw_post_profile.r',
            'draw_prep': 'draw_prep.r',
        }

    @staticmethod
    def _find_software(sfw_name):
        # 优先检查 ../dependencies 目录
        default_path = os.path.join(os.path.dirname(
            os.path.realpath(__file__)), '../dependencies', sfw_name)
        if os.path.exists(default_path) and os.access(default_path, os.X_OK):
            return default_path

        # 如果 ../dependencies 中没有，检查 PATH 中的可执行文件
        sfw_path = shutil.which(sfw_name)
        if sfw_path:
            if os.path.exists(sfw_path) and os.access(sfw_path, os.X_OK):
                return sfw_path
            else:
                logger.error(
                    f"{sfw_name} found at '{sfw_path}' but it is not executable.")
                return None

        # 如果都找不到，记录错误
        logger.error(
            f"{sfw_name} not found in '{default_path}' or in PATH, and neither is executable.")
        return None

    def check_dependency(self, dependency):
        dependency_path = self._find_software(self.externals[dependency])
        if not dependency_path:
            logger.error(f"Dependency [{dependency}] not managed.")
            raise ValueError(f'Dependency [{dependency}] not managed.')

        dependency_path = os.path.abspath(dependency_path)
        setattr(self, dependency, dependency_path)
        logger.info(f"Dependency {dependency} is checked: {dependency_path}")
        return dependency_path


sfw = _Externals()


def set_golbal(logger_):
    global logger
    logger = logger_


def set_verbosity_level(outdir, verbose=False, debug=False, tag='main'):
    log_format = "{time:YYYY-MM-DD HH:mm:ss} | {level} | {module} : {message}"

    logger.remove()
    logger.add(lambda msg: tqdm.write(msg, end=''),
               format=log_format, level="DEBUG" if verbose else "INFO")
    set_golbal(logger)
    logger.info("Set the output directory to {}".format(outdir))
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    logger.info(f"Set verbosity level: {'DEBUG' if verbose else 'INFO'}")
    if os.path.exists(f"{outdir}/{tag}.log"):
        os.remove(f"{outdir}/{tag}.log")

    if debug:
        logger.add(f"{outdir}/{tag}.log", level="TRACE")
    else:
        if verbose:
            logger.add(f"{outdir}/{tag}.log", level="DEBUG")
        else:
            logger.add(f"{outdir}/{tag}.log", level="INFO")

    logger.info(
        f"Please check the intact log in {tag}.log.")
    return outdir


def run_command(cmd):
    logger.info(f'Run: [{cmd}]')
    code, output = subprocess.getstatusoutput(cmd=cmd)
    if int(code) == 0:
        ...
    else:
        logger.error(f'Error in running command: {cmd}')
        logger.error(f'{output}')
        raise Exception(f'Error in running command: {cmd}')

    return output
