#!/usr/bin/env python3

import argparse
import os
import logging
import sys


class ModuleDetectionArgs:

    def __init__(self, args: list = None):
        """
        Initiator serving as a main function
        :param args: Optional list of arguments; e.g. ["-i", "inputValue", ...]
        """
        # Collect version and set log
        self._module_path = os.sep.join(os.path.dirname(os.path.realpath(__file__)).split(os.sep)[:-2])
        self.__version__ = self._get_version()
        self._log = logging.getLogger('root')
        # Collect args
        self._args = vars(self._parse_arguments(args))
        # Set and validate args
        self._read = self._get_value("read")
        self._read_type = self._get_read_type(self._read)
        self._virus_genome = self._get_value("virusGenome")
        self._human_genome = self._get_value("humanGenome")
        self._out = self._get_value("out")
        self._virus_threshold = self._get_value("virusThreshold")
        self._virus_chr_ref = self._get_value("virusChrRef")
        self._thread = self._get_value("thread")
        self._read_files_command = self._get_value("readFilesCommand")
        self._single_cell_barcode = self._get_value("singleCellBarcode")
        self._seq_resolution = self._get_resolution(self._single_cell_barcode)
        self._single_unique_mol_ident = self._get_value("singleUniqueMolIdent")
        self._single_inclusion_list = self._get_value("singleInclusionList")
        self._sensitivity = self._get_value("sensitivity")
        self._human_star_parameters = self._get_value("humanSTARparameters")
        self._virus_star_parameters = self._get_value("virusSTARparameters")

    def _get_version(self) -> str:
        """
        Read version file in the main directory and return the value
        :return: Version number v#.#.#
        """
        with open(self._module_path + '/__version__', 'r') as ver:
            return ver.readline().strip()

    @staticmethod
    def _set_log_level(log_level: str) -> None:
        """
        Decide which level of log is displayed in the terminal
        :param log_level: Key of the the minimal level of log the user wish to receive
        :return: None
        """
        log = logging.getLogger('root')
        log_level = log_level.upper()
        # use for python 3.8
        if sys.version_info.major * 10 + sys.version_info.minor >= 38:
            logging.basicConfig(format='%(asctime)s %(module)s:%(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S',
                                level=logging.getLevelName(log_level),
                                force=True
                                )
        else:
            # use for python < 3.6
            log.setLevel(logging.getLevelName(log_level))
        log.info(f"Log level has been set to: {log_level}")

    def _get_value(self, key: str):
        """
        Extract values from self._args. If key doesn't exist (user did not provide value), return None
        :param key: Key name
        :return: Value of key or None
        """
        try:
            return self._args[key]
        except KeyError:
            return None

    @staticmethod
    def _get_read_type(read: list) -> str:
        """
        Based on the length of the read value return if it is single_end or paired_end reads
        :param read: List of reads
        :return: String with either 'single_end' or 'paired_end'
        """
        return "single_end" if len(read) == 1 else "paired_end"

    @staticmethod
    def _get_resolution(single_sample_barcode: list) -> str:
        """
        Returns the resolution of the sequence
        :param single_sample_barcode: list of single sample barcode
        :return: if single_sample_barcode is defined return 'single_cell' else 'bulk'
        """
        return "single_cell" if single_sample_barcode else "bulk"

    def _parse_arguments(self, args: list = None) -> argparse.Namespace:
        """
        Parses input values
        :param args: Optional list of arguments; e.g. ["-i", "inputValue", ...]
        :return: argparse.Namespace object
        """
        parser = argparse.ArgumentParser(description="VENUS, a subtractive analysis software: " +
                                                     "Virus dEtecting in humaN bUlk and Single cell rna sequencing")

        parser.add_argument("--read", type=str, required=True, nargs='+',
                            help="read of RNA-seq \n(single-cell) first read should be cDNA, second should be CB+UMI")

        parser.add_argument("--virusGenome", type=str, required=True,
                            help="directory path of virus genome index")

        parser.add_argument("--humanGenome", type=str, required=True,
                            help="directory path of human genome index")

        parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                            help="directory path of output dir")

        parser.add_argument("--virusThreshold", type=str, required=False, default=0,
                            help="viral load threshold to filter out negligible viruses")

        parser.add_argument("--virusChrRef", type=str, required=True,
                            help="tsv file to map NC_* id to virus species name, e.g. NC_001802.1 --> HIV-1")

        parser.add_argument("--thread", type=int, required=False, default=1,
                            help="number of parallel threads")

        parser.add_argument("--readFilesCommand", type=str, required=False,
                            help="uncompression command")

        parser.add_argument("--singleCellBarcode", type=str, required=False, nargs=2,
                            help="(single-cell) barcode specifications, e.g. '1 16' = start at 1, length is 16")

        parser.add_argument("--singleUniqueMolIdent", type=str, required=False, nargs=2,
                            help="(single-cell) umi specifications, e.g. '17 10' = start at 17, length is 10")

        parser.add_argument("--singleInclusionList", type=str, required=False,
                            help="(single-cell) barcode inclusion")

        parser.add_argument("--sensitivity", type=str, required=False, default="low",
                            help="sensitivity vs accuracy trade off, default option is low sensitivity, high accuracy")

        parser.add_argument("--humanSTARparameters", type=str, required=False,
                            help="a way to pass an argument to STAR in the human genome mapping stage " +
                                 "(e.g. --humanSTARparameters '--limitOutSJcollapsed 2000000')")

        parser.add_argument("--virusSTARparameters", type=str, required=False,
                            help="a way to pass an argument to STAR in the virus genome mapping stage " +
                                 "(e.g. --virusSTARparameters '--limitOutSJcollapsed 2000000')")

        parser.add_argument('-V', '-version', '--version', action='version', version=self.__version__)

        log_settings = parser.add_argument_group("Log settings")
        log_settings.add_argument('--logLevel', action='store', default="INFO",
                                  choices=["DEBUG", "WARNING", "INFO", "ERROR"],
                                  help='If set will determine the level of logging to be used. eg DEBUG, WARNING, INFO')

        return parser.parse_args(args=args)

    @property
    def read(self) -> list:
        return self._read

    @property
    def read_type(self) -> str:
        return self._read_type

    @property
    def seq_resolution(self) -> str:
        return self._seq_resolution

    @property
    def virus_genome(self) -> str:
        return self._virus_genome

    @property
    def human_genome(self) -> str:
        return self._human_genome

    @property
    def out(self) -> str:
        return self._out

    @property
    def virus_threshold(self) -> str:
        return self._virus_threshold

    @property
    def virus_chr_ref(self) -> str:
        return self._virus_chr_ref

    @property
    def thread(self) -> int:
        return self._thread

    @property
    def read_files_command(self) -> str:
        return self._read_files_command

    @property
    def single_cell_barcode(self) -> list:
        return self._single_cell_barcode

    @property
    def single_unique_mol_ident(self) -> list:
        return self._single_unique_mol_ident

    @property
    def single_inclusion_list(self) -> str:
        return self._single_inclusion_list

    @property
    def sensitivity(self) -> str:
        return self._sensitivity

    @property
    def human_star_parameters(self) -> str:
        return self._human_star_parameters

    @property
    def virus_star_parameters(self) -> str:
        return self._virus_star_parameters

    @property
    def log(self) -> logging.Logger:
        return self._log
