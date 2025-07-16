from loguru import logger
import numpy as np
import pickle

"""

Pickle check class for validating and loading pickled data.
This class is used to store and manage pickled data, ensuring that the data can be loaded correctly
and that the parameters match the expected values.
It provides methods for loading data, checking parameters, and dumping the data back to a pickle file.

"""


class PklCheck():
    def __init__(self, outdir: str, name: str) -> None:
        self.main_data = {}
        self.parameter = {}
        self.name = name
        self.outdir = outdir
        self.status = 0

    # def update_step(self,step:int):
    #     '''
    #     status:
    #         0 for not start,
    #         1 for preprocess finished load data and finished clust,
    #         2 for partition finished,
    #         3 for finished
    #     '''
    #     status_desc = {0: 'not start', 1: 'preprocess finished', 2: 'partition finished', 3: 'postprocess finished'}
    #     if step < self.status:
    #         logger.warning(f'If you want to re-start program, please delete the previous output files. Current step is {status_desc[self.status]}, however, you want to start from {status_desc[step]} that mignt cause some error')
    #     self.status = step

    def load(self, section: str, main_data: dict, parameter: dict = {}):
        self.main_data[section] = main_data
        if parameter:
            self.parameter.update(parameter)

    def pickle_(self):
        logger.info(f'Pickle the data to {self.outdir}/{self.name}.pkl')
        for k in self.main_data.keys():
            logger.info(f'Object: {k}')
        logger.info('-----------------')
        for k, v in self.parameter.items():
            logger.info(f'Key: {k} = {v}')
        with open(f'{self.outdir}/{self.name}.pkl', 'wb') as f:
            pickle.dump(self, f)

    def decode(self, **kwargs):
        all_keys = set(self.parameter.keys())
        for k, v in kwargs.items():
            if k in self.parameter:
                all_keys.remove(k)
                value1 = v
                value2 = self.parameter[k]
                if isinstance(value1, np.ndarray) and isinstance(value2, np.ndarray):
                    if not np.array_equal(value1, value2):
                        logger.warning(
                            f'The parameter {k} is not the same as previous data')
                        return False
                elif value1 != value2:
                    if k == 'retrieve':
                        if value1 is not True and value2 is True:  # if retrieve is false in partition step it is ok
                            continue
                    logger.warning(
                        f'The parameter {k} is not the same as previous data: {value1} vs {value2}')
                    return False
            else:
                logger.warning(
                    f'Cannot find the paramenter [{k}] in previous data')
                raise ValueError(
                    f'Cannot find the paramenter [{k}] in previous data')
        if all_keys:
            logger.warning(
                f'The following keys are not used in the decode process: {all_keys}')
            raise ValueError(
                f'The following keys are not used in the decode process: {all_keys}. Tell me in the github issue')
        return True

    def data_dump(self, name):
        if name in self.main_data:
            return self.main_data[name]
        else:
            logger.error(
                f'This pkl {self.name} is not you expected one {name}')
            raise ValueError(
                f'This pkl {self.name} is not you expected one {name}')
