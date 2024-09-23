import glob
import os
import sys
import warnings
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from ecoli_methods import Methods


__author__ = 'castonguaya'
__version__ = '0.1'


#  Requirements
"""
conda create -n ecoli -y -c bioconda -c etetoolkit -c hcc -c conda-forge samtools=1.16 bcftools=1.16 git ete3=3.1.2 psutil=5.9.4 pandas=1.5.3 rebaler=0.2.0 perl-bioperl=1.7.8 perl=5.32.1
"""


class Bacon(object):
    def __init__(self, args):
        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.ref_genome = args.genome

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        # Check input file compatibility
        Methods.check_input(self.input)

        ############################################################

        # Step completion report files
        done_extract = self.done_extract + '/done_extract'
        done_result = self.done_extract + '/done_result'

        # Output folders to create
        extract_folder = self.output_folder + '/1_extract/'
        result_folder = self.output_folder + '/2_result/'

        # Create output folder
        Methods.make_folder(self.output_folder)

        print('\tAll good!')

        ##################
        #
        # 1- Position primer
        #
        ##################

        if not os.path.exists(done_extract):
            print('test1...')
            Methods.test1()
            Methods.flag_done(done_extract)
        else:
            print('Skipping extract. Already done.')

        # Update sample_dict after extracting
        self.position_genome = Methods.get_files(extract_folder)

        ##################
        #
        # 2- result position
        #
        ##################

        if not os.path.exists(done_result):
            print('test2...')
            Methods.test2()
            Methods.flag_done(done_result)
        else:
            print('Skipping result. Already done.')

        print('DONE!')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Extract, assemble and compare Nanopore reads matching a reference sequence.')
    parser.add_argument('-g', '--genome', metavar='/path/to/reference_organelle/genome.fasta',
                        required=True,
                        help='Reference genome for primer mapping. Mandatory.')
    parser.add_argument('-i', '--input', metavar='/path/to/input/folder/ or /path/to/my_fasta',
                        required=True,
                        help='Folder that contains the fasta files or individual fasta file. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Bacon(arguments)