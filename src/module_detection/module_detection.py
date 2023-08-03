#!/usr/bin/env python3

# ----------------------------------------------------------------------------------------------------------------------
# Venus: Virus infection detection and integration site discovery method using single-cell RNA-seq
# Detection module
#
# (C) 2022 Che Yu Lee, Irvine, California
# Released under GNU Public License (GPL)
# email cheyul1@uci.edu
#
# ----------------------------------------------------------------------------------------------------------------------
#
# Modified by Alex Janse, Princess Máxima Center, Utrecht, the Netherlands
# email a.janse-3@prinsesmaximacentrum.nl
# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# Import Libraries
########################################################################################################################
import os
import pandas as pd
import pathlib
from src.module_detection.module_detection_args import ModuleDetectionArgs


class ModuleDetection:

    def __init__(self, args: list = None):
        self._args = ModuleDetectionArgs(args)
        self._log = self._args.log
        self.quality_control()
        self.map_human()
        self.map_virus()
        self.output_infection()

    def _run_command(self, cmd: str, log_tool: str) -> None:
        """
        Will run the command string in the background
        :param cmd: Command line to run
        :param log_tool: Tool to mention in the log
        :return: None
        """
        self._log.debug(f"The following command will be performed:\n{cmd}")
        self._log.info(f"Starting {log_tool}")
        # Run the command
        os.system(cmd)
        self._log.info(f"{log_tool} has finished")

    def quality_control(self) -> None:
        """Trims bad quality sequences"""
        output_dir = self._args.out + "/quality_control/"
        self._log.info("Preparing Quality Control")
        pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Creating command line to call trim_galore in the background
        cmd = "trim_galore " \
              "--trim-n " \
              "--quality 5 " \
              "--phred33 " \
              "--length 20 " \
              f"--output_dir {output_dir} " \
              "--gzip "

        cmd += "--cores "
        cmd += "7 " if self._args.thread >= 7 else f"{self._args.thread} "

        if self._args.read_type == "single_end":
            cmd += self._args.read[0] + " "
        elif self._args.read_type == "paired_end":
            cmd += f"--paired {self._args.read[0]} {self._args.read[1]} "

        self._run_command(cmd=cmd, log_tool="Quality Control")

    def map_human(self) -> None:
        """Maps to the human genome"""
        self._log.info("Preparing mapping to human reference")

        # Build command line to run in the background
        cmd = "STAR " \
              f"--runThreadN {self._args.thread} " \
              f"--outFileNamePrefix {self._args.out}/human/ " \
              f"--genomeDir {self._args.human_genome} " \
              "--outSAMtype None " \
              "--outReadsUnmapped Fastx "

        if self._args.read_files_command:
            cmd += f"--readFilesCommand {self._args.read_files_command} "

        if self._args.read_type == "single_end":
            cmd += f"--readFilesIn {self._args.out}/quality_control/*_trimmed.fq.gz "
        elif self._args.read_type == "paired_end":
            cmd += f"--readFilesIn {self._args.out}/quality_control/*_val_1.fq.gz " \
                   f"{self._args.out}/quality_control/*_val_2.fq.gz "

        if self._args.seq_resolution == "single_cell":
            cmd += "--soloType CB_UMI_Simple " \
                   f"--soloCBwhitelist {self._args.single_inclusion_list} " \
                   "--soloBarcodeReadLength 0 " \
                   f"--soloCBstart {self._args.single_cell_barcode[0]} " \
                   f"--soloCBlen {self._args.single_cell_barcode[1]} " \
                   f"--soloUMIstart {self._args.single_unique_mol_ident[0]} " \
                   f"--soloUMIlen {self._args.single_unique_mol_ident[1]} "

        if bool(self._args.human_star_parameters):
            cmd += self._args.human_star_parameters + " "

        self._run_command(cmd=cmd, log_tool="Human Mapping")

    def prep_human(self):
        """ Renames mate# files """
        if self._args.read_type == "single_end":
            os.rename(f"{self._args.out}/human/Unmapped.out.mate1", f"{self._args.out}/human/Unmapped.out.mate1.fastq")
        elif self._args.read_type == "paired_end":
            os.rename(f"{self._args.out}/human/Unmapped.out.mate1", f"{self._args.out}/human/Unmapped.out.mate1.fastq")
            os.rename(f"{self._args.out}/human/Unmapped.out.mate2", f"{self._args.out}/human/Unmapped.out.mate2.fastq")

    def map_virus(self) -> None:
        """Maps the leftover human reads to the virus genome"""
        self._log.info("Preparing mapping to virus references")
        self.prep_human()   # Make sure the files have the correct file name

        # Build the command
        cmd = "STAR " \
              f"--runThreadN {self._args.thread} " \
              f"--outFileNamePrefix {self._args.out}/virus/ " \
              f"--genomeDir {self._args.virus_genome} " \
              "--outFilterMultimapNmax 1 "

        if self._args.read_type == "single_end":
            cmd += f"--readFilesIn {self._args.out}/human/Unmapped.out.mate1.fastq "
        elif self._args.read_type == "paired_end":
            cmd += f"--readFilesIn {self._args.out}/human/Unmapped.out.mate1.fastq " \
                   f"{self._args.out}/human/Unmapped.out.mate2.fastq "

        if self._args.seq_resolution == "single_cell":
            cmd += "--soloType CB_samTagOut " \
                   f"--soloCBwhitelist {self._args.single_inclusion_list} " \
                   "--soloCBmatchWLtype 1MM " \
                   "--soloBarcodeReadLength 0 " \
                   f"--soloCBstart {self._args.single_cell_barcode[0]} " \
                   f"--soloCBlen {self._args.single_cell_barcode[1]} " \
                   f"--soloUMIstart {self._args.single_unique_mol_ident[0]} " \
                   f"--soloUMIlen {self._args.single_unique_mol_ident[1]} " \
                   "--outSAMtype BAM Unsorted " \
                   "--outSAMattributes NH HI nM AS CR UR "
        elif self._args.seq_resolution == "bulk":
            cmd += "--outSAMtype SAM "

        if self._args.sensitivity == "low":
            cmd += "--outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 "  # Default
        else:
            cmd += "--outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 "

        if bool(self._args.virus_star_parameters):
            cmd += self._args.virus_star_parameters + " "

        self._run_command(cmd=cmd, log_tool="Virus Mapping")

    def output_infection(self) -> pd.DataFrame:
        """Produces the detection output file for mega-virus mode"""
        self._log.info("Creating virus infection report")

        if self._args.seq_resolution == "single_cell":
            os.system(f"samtools view -h -o {self._args.out}/virus/Aligned.out.sam "
                      f"{self._args.out}/virus/Aligned.out.bam")

        sam_file = f"{self._args.out}/virus/Aligned.out.sam"
        with open(sam_file, "r") as file:
            self._log.debug(f"Processing {sam_file}")
            virus_species = []
            for line in file:
                if not line.startswith("@"):
                    virus_species.append(line.split("\t")[2])
            total = len(virus_species)
            virus_species = list(set(virus_species))

        self._log.debug(f"{total} hits and {len(virus_species)} unique hits")

        species_name = []
        species_count = []
        species_ratio = []
        species_barcodes = []
        reference = pd.read_csv(self._args.virus_chr_ref, sep='\t', names=["id", "name"])
        reference = reference.set_index("id")

        # # Testing
        # stopper = 0
        for species in virus_species:
            virus_count = 0
            virus_barcodes = []
            with open(f"{self._args.out}/virus/Aligned.out.sam", "r") as file:
                for line in file:
                    if (not line.startswith("@")) and species in line:
                        virus_count += 1
                        if self._args.seq_resolution == "single_cell":
                            virus_barcodes.append(line.split("\t")[-2].split(sep=":")[-1])
            if virus_count >= int(self._args.virus_threshold):
                species_count.append(virus_count)
                species_name.append(reference.loc[species, "name"])
                species_ratio.append(virus_count / total * 100.0)
                species_barcodes.append(virus_barcodes)

        output_dict = {"Name": species_name,
                       "Count": species_count,
                       "Percentage": species_ratio}

        if self._args.seq_resolution == "single_cell":
            output_dict.update({"Barcodes": species_barcodes})

        output = pd.DataFrame(data=output_dict)

        output.sort_values(by="Percentage", ascending=False, inplace=True)
        output_report = f"{self._args.out}/detection_output.tsv"
        self._log.info(f"Writing output to {output_report}")
        output.to_csv(path_or_buf=output_report, sep="\t", index=False)

        return output


def main():
    md = ModuleDetection()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
