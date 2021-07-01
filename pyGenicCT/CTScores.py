from miscSupports import load_yaml, directory_iterator, validate_path, flip_list, terminal_time, FileOut
from csvObject import CsvObject, write_csv
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from pathlib import Path
import pandas as pd
import numpy as np


class CTScores:
    def __init__(self, yaml_file):

        # Load yaml args and set logger
        self.args = load_yaml(yaml_file)
        self.logger = FileOut(self.args["write_path"], "logger", "log", True)

        # Set parse indexes
        self._chr_index = self.args["chromosome_index"]
        self._snp_index = self.args["snp_index"]
        self._coe_index = self.args["coefficient_index"]
        self._p_v_index = self.args["p_value_index"]

        # Validate genetic directories and base name of the genetic files
        custom_meta_path(validate_path(self.args["meta_path"]))
        self._root = validate_path(self.args["gen_path"])
        self._base_name = self.args["base_name"]

        # Load the valid snps from plink
        self._iid = Bgen(self.get_file_name(1)).iid

        # Set the headers based on the iid fid combination
        if len(self._iid[0]) == 2:
            self._headers = ["FID", "IID"]
        elif len(self._iid[1]) == 1:
            self._headers = ["IID"]
        else:
            raise IndexError("Headers expect IID / FID")

        # Setup a scores header as an attribute so we can call it rather than pass it back and forth via methods
        self._scores_holder = []

    def link_resources(self):
        """
        Link the valid and the values files via pandas
        """
        # Load the valid snps
        valid_snps = pd.read_csv(self.args["Valid"], sep=" ")

        # Load the snp data
        data = pd.read_csv(self.args["Values"], sep=" ")

        # Merge the two, drop duplicate, column name, then write to file
        df = valid_snps.merge(data, left_on=self.args["valid_snp_name"], right_on=self.args["values_snp_name"])
        df = df.drop(columns=self.args["values_snp_name"])
        df.to_csv(Path(self.args["write_directory"], "Snps.csv"), index=False)

    def create_score_levels(self):
        """
        Create score from the sum of the dosage effect * snp beta for each p level threshold specified by the user

        :return:
        """

        # set the p value thresholds scores for each threshold
        [self.create_score(threshold) for threshold in self.args["threshold"][::-1]]

        # Write the scores out
        self._scores_holder = [list(iid) + values for iid, values in zip(self._iid, flip_list(self._scores_holder))]
        write_csv(self.args["write_path"], self.args["write_name"], self._headers, self._scores_holder)
        self.logger.write(f"Finished at {terminal_time()}")
        self._scores_holder = []

    def create_score(self, threshold):
        """
        Create score from the sum of the dosage effect * snp beta

        :param threshold: Which p value threshold to extract snps on
        :type threshold: float

        :return: Nothing, append snps tp pit then stop out then stop
        """

        total_scores = np.zeros(len(self._iid))
        total = 0
        for i in range(1, 23):
            # Load the gen reference, extract snps that are within this chromosome and meet this threshold
            gen_file = Bgen(self.get_file_name(i))
            snp_effects = self.chromosome_thresholds(threshold, i)

            if len(snp_effects) > 0:
                # Index the gen file with the snp names
                indexes_snps = self.get_snp_indexes(snp_effects, gen_file)
                self.logger.write(f"Found {len(indexes_snps)} out of {len(snp_effects)} for chromosome {i} at "
                                  f"{terminal_time()}")

                # Extract the dosage and then multiple those values by the coefficients, add it to score
                gen_file = gen_file[:, [index for _, _, index in indexes_snps]]
                for (_, effect, _), snp_dosage in zip(indexes_snps, self._extract_dosage(gen_file)):
                    total_scores += effect * snp_dosage
                total += len(indexes_snps)

            else:
                self.logger.write(f"No snps for chromosome {i} valid at this threshold {threshold} at "
                                  f"{terminal_time()}")

        if sum(total_scores) != 0:
            self.logger.write(f"Finished threshold {threshold} with {total} snps at {terminal_time()}\n\n")
            self._headers.append(str(threshold))
            self._scores_holder.append(total_scores.tolist())
        else:
            self.logger.write(f"Found No valid snps for threshold {threshold} at {terminal_time()}\n\n")

    def get_file_name(self, target):
        """
        Get the bgen file for this chromosome

        :param target: Target chromosome
        :type target: int

        :return: The path to this file
        :rtype: Path
        """
        for file in directory_iterator(self._root):
            try:
                chromosome = int(file.split(self._base_name)[1].split(".")[0])
                if (chromosome == target) and (".bgen" in file) and (".bgi" not in file):
                    return Path(self._root, file)

            except (TypeError, ValueError):
                pass

    def chromosome_thresholds(self, threshold, i):
        """
        Extract the snps that are in this chromosome that also meets the p-value threshold

        :param threshold: P value threshold to validate snp against
        :type threshold: float

        :param i: Chromosome index
        :type i: index

        :return: List of snp_name, coefficient value
        :rtype: list[str, float]
        """
        valid_snps = CsvObject(Path(self.args["write_path"], "Snps.csv"))

        snp_effects = []
        for row in valid_snps.row_data:
            if int(float(row[self._chr_index])) == i and (0 < float(row[self._p_v_index]) <= threshold):
                snp_effects.append([row[self._snp_index], float(row[self._coe_index])])

        return snp_effects

    @staticmethod
    def get_snp_indexes(snp_effects, gen_file):
        """
        Get the index for each snp within the gen file

        :param snp_effects: The snp-effects list
        :type snp_effects: list[str, float]

        :param gen_file: The genetic file
        :type gen_file: Bgen

        :return: A list of Snp, effect, index
        :rtype: list[str, float, int]
        """

        # Get {variant id: rs_id} for each snp, so we can reconstruct 'variant_id,rs_id' for py-snp-tools
        v_dict = {snp[1]: snp[0] for snp in [snp.split(",") for snp in gen_file.sid]}

        # Reconstruct 'variant_id,rs_id' and if it exists add to snp_indexes as [snp, (beta) effect, snp_index]
        snp_indexes = []
        for snp, effect in snp_effects:
            try:
                snp_indexes.append([snp, effect, gen_file.sid_to_index([f"{v_dict[snp]},{snp}"]).tolist()[0]])
            except KeyError:
                pass

        return snp_indexes

    @staticmethod
    def _extract_dosage(gen_file):
        """Extract the dosage"""
        return sum(np.array([snp * i for i, snp in enumerate(gen_file.read(dtype=np.int8).val.T)], dtype=np.int8))
