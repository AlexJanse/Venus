# ----------------------------------------------------------------------------------------------------------------------
# Venus: Virus infection detection and integration site discovery method using single-cell RNA-seq
# Integration module
#
# (C) 2022 Che Yu Lee, Irvine, California
# Released under GNU Public License (GPL)
# email cheyul1@uci.edu
# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# Import Libraries
########################################################################################################################
import argparse
import os
import pysam


########################################################################################################################
# Main
########################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="VENUS, a subtractive analysis software: " + \
                                                 "Virus dEtecting in humaN bUlk and Single cell rna sequencing")

    parser.add_argument("--read", type=str, required=True,
                        help="read of RNA-seq \n(single-cell) first read should be cDNA, second should be CB+UMI")

    parser.add_argument("--virusGenome", type=str, required=True,
                        help="directory path of virus genome index")

    parser.add_argument("--hybridGenome", type=str, required=True,
                        help="directory path of hybrid genome index")

    parser.add_argument("--guideFASTA", type=str, required=False,
                        help="Fasta file for the guide sequence(s)")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory path of output dir")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    parser.add_argument("--readFilesCommand", type=str, required=False,
                        help="uncompression command")

    args = parser.parse_args()

    def quality_control():
        """Trims bad quality sequences"""
        cmd = "trim_galore " \
              + "--trim-n " \
              + "--quality 5 " \
              + "--phred33 " \
              + "--length 20 " \
              + "--output_dir " + args.out + "/quality_control/" + " " \
              + "--gzip " \
              + "--cores " + args.thread + " " \
              + args.read

        return cmd

    def map_virus():
        """Maps to the virus genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus/ " \
              + "--genomeDir " + args.virusGenome + " "

        if args.readFilesCommand is not None:
            cmd = cmd + "--readFilesCommand " + args.readFilesCommand + " "

        cmd = cmd + "--readFilesIn " + args.read + " "

        cmd = cmd \
              + "--outSAMmultNmax 1 " \
              + "--alignIntronMax 1 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 " \
                    + "--scoreDelBase -1 --scoreInsBase -1 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 "

        cmd = cmd \
              + "--outSAMtype BAM Unsorted "

        return cmd

    def map_hybrid():
        """Maps the viral mapped reads to the human genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/hybrid/ " \
              + "--genomeDir " + args.hybridGenome + " "

        cmd = cmd + "--readFilesIn " + args.read + " "

        cmd = cmd \
              + "--outSAMmultNmax 1 " \
              + "--alignIntronMin 20 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 " \
                    + "--scoreDelBase -1 --scoreInsBase -1 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 "

        cmd = cmd \
              + "--chimOutType WithinBAM SoftClip --chimSegmentMin 12 --chimJunctionOverhangMin 8 " \
              + "--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 " \
                    + "--chimNonchimScoreDropMin 10 " \
              + "--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 " \
                    + "--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 " \
              + "--alignInsertionFlush Right "

        cmd = cmd \
              + "--outSAMtype BAM Unsorted "

        return cmd

    def index_guide():
        """Generate the small genome index for guide sequence"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/guide/index/ " \
              + "--runMode genomeGenerate " \
              + "--genomeDir " + args.out + "/guide/index/ " \
              + "--genomeFastaFiles " + args.guideFASTA + " " \
              + "--genomeSAindexNbases 2 "

        return cmd

    def map_guide():
        """Classifies the integration sites based on guide sequences"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/guide/ " \
              + "--genomeDir " + args.out + "/guide/index/ " \
              + "--outFilterMultimapNmax 1 "

        cmd = cmd + "--readFilesIn " + args.read + " "

        cmd = cmd \
              + "--alignIntronMax 1 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNmin 25 --outFilterMatchNminOverLread 0 " \
                    + "--outFilterScoreMin 20 --outFilterScoreMinOverLread 0 "

        # cmd = cmd \
        #       + "--outSAMmultNmax 1 " \
        #       + "--alignIntronMax 1 --winBinNbits 7 " \
        #       + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 " \
        #             + "--scoreDelBase -1 --scoreInsBase -1 " \
        #       + "--outFilterMultimapNmax -1 -winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
        #             + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
        #       + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 "

        cmd = cmd \
              + "--outReadsUnmapped None "

        return cmd

    ##################################################################################################
    # Action Steps
    ##################################################################################################
    print(quality_control())
    print(map_virus())
    print(map_hybrid())
    print(index_guide())
    print(map_guide())
    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
