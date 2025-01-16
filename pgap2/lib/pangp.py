import numpy as np
import pandas as pd
import random
import math
import numpy as np

from loguru import logger
from itertools import combinations


class PanGP():
    def __init__(self, pav: pd.DataFrame, S, R, K, method) -> None:
        self.pav = pav
        self.coef_k = K
        self.strain_num = pav.shape[1]
        self.sampling_num = S
        self.repeat = R
        self.known_distance = np.full((pav.shape[1], pav.shape[1]), np.nan)
        np.fill_diagonal(self.known_distance, 1)
        self.sampling_method = method
        self.infalted_sampling_num = self.sampling_num*self.coef_k
        self.deflated_sampling_num = 0
        pass

    @staticmethod
    def lazy_random_sample(n, k):
        seen = set()
        for _ in range(k):
            x = random.randint(0, n-1)
            while x in seen:
                x = random.randint(0, n-1)
            seen.add(x)
            yield x

    def _get_sample_combinations(self, k, sampling_num):
        expect_combination_num = math.comb(self.strain_num, k)
        self.deflated_sampling_num = expect_combination_num if expect_combination_num < self.sampling_num else self.sampling_num
        sampling_num = expect_combination_num if expect_combination_num < sampling_num else sampling_num
        if sampling_num == expect_combination_num:
            # return all combinations
            return list(combinations(range(self.strain_num), k))
        else:
            # selected_indices = random.sample(range(expect_combination_num), sampling_num)
            # selected_indices = self.lazy_random_sample(
            #     expect_combination_num, sampling_num)  # random sample and avoid memory error
            # result = []
            # for index in selected_indices:
            #     n = self.strain_num
            #     this_comb = []
            #     for i in range(k):
            #         for j in range(n):
            #             comb = math.comb(n - j - 1, k - i - 1)
            #             if comb > index:
            #                 this_comb.append(j)
            #                 n = n - j - 1
            #                 break
            #             index -= comb
            #     this_comb = [x + i for i, x in enumerate(this_comb)]
            #     result.append(this_comb)
            result = []
            for _ in range(sampling_num):
                # 每次从 range(N) 中抽取 K 个不重复的样本
                sample = np.array(random.sample(range(self.strain_num), k))
                result.append(sample)
            return result

    def sampling(self, strain_i: int):
        assert strain_i <= self.strain_num, logger.error(
            f"Strain index {strain_i} out of range. Please tell me the detail in github")
        if self.sampling_method == 'TR':

            selected_combinations = self._get_sample_combinations(
                strain_i, self.sampling_num)

        elif self.sampling_method == 'DG':
            total_combinations = self._get_sample_combinations(
                strain_i, self.infalted_sampling_num)
            diversity_scores = []
            for each_combination in total_combinations:
                dev_gene_cluster = 0
                if len(each_combination) == 1:
                    dev_gene_cluster = 1
                else:
                    dev_gene_cluster = self.calculate_dev_gene_cluster(
                        each_combination)
                diversity_scores.append((each_combination, dev_gene_cluster))
            sorted_combinations = sorted(
                diversity_scores, key=lambda x: x[1], reverse=True)
            selected_combinations = [sorted_combinations[i][0] for i in range(0, len(
                sorted_combinations), len(sorted_combinations) // self.deflated_sampling_num)]
        return selected_combinations

    def get_diff_clusters(self, strain_i: int, strain_j: int):
        diff_cluster_num = 0
        if not np.isnan(self.known_distance[strain_i][strain_j]):
            diff_cluster_num = self.known_distance[strain_i][strain_j]
        else:
            # strain_i_name = self.pav.columns[strain_i]
            strain_i_name = self.pav[:, strain_i]
            # strain_j_name = self.pav.columns[strain_j]
            strain_j_name = self.pav[:, strain_j]
            # set1 = set(self.pav.index[self.pav[strain_i_name] >= 1])
            set1 = set(np.where(self.pav[:, strain_i_name] >= 1)[0])
            # set2 = set(self.pav.index[self.pav[strain_j_name] >= 1])
            set2 = set(np.where(self.pav[:, strain_j_name] >= 1)[0])
            diff_cluster_num = len(set1.symmetric_difference(set2))
            self.known_distance[strain_i][strain_j] = diff_cluster_num
            self.known_distance[strain_j][strain_i] = diff_cluster_num
        return diff_cluster_num

    def calculate_dev_gene_cluster(self, combination):
        total_diff_clusters = 0
        pair_count = 0
        for i in range(len(combination)):
            for j in range(i + 1, len(combination)):
                strain_i = combination[i]
                strain_j = combination[j]
                diff_cluster_num = self.get_diff_clusters(strain_i, strain_j)
                total_diff_clusters += diff_cluster_num
                pair_count += 1
        return total_diff_clusters / pair_count
