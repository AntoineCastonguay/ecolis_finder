import os
from argparse import ArgumentParser
from ecoli_methode import Methods


__author__ = 'castonguaya'
__version__ = '0.1'


class Ecoli(object):
    def __init__(self, args):
        # I/O
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.ref_genome = args.genome
        self.run()

    def run(self):
        print('Checking a few things...')
        Methods.check_input(self.input)

        done_extract = self.output_folder + '/done_extract'
        done_result = self.output_folder + '/done_result'
        extract_folder = os.path.join(self.output_folder, '1_extract')
        result_folder = os.path.join(self.output_folder, '2_result')

        Methods.make_folder(self.output_folder)
        print('\tAll good!')

        if not os.path.exists(done_extract):
            print('test1...')
            Methods.test1()
            Methods.flag_done(done_extract)
        else:
            print('Skipping extract. Already done.')

        #self.position_genome = Methods.get_files(extract_folder)

        if not os.path.exists(done_result):
            print('test2...')
            Methods.test2()
            Methods.flag_done(done_result)
        else:
            print('Skipping result. Already done.')

        print('DONE!')

if __name__ == "__main__":
    parser = ArgumentParser(description='Extract, assemble and compare Nanopore reads matching a reference sequence.')
    parser.add_argument('-g', '--genome', required=True, help='Reference genome for primer mapping. Mandatory.')
    parser.add_argument('-i', '--input', required=True, help='Folder that contains the fasta files or individual fasta file. Mandatory.')
    parser.add_argument('-o', '--output', required=True, help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    
    arguments = parser.parse_args()
    Ecoli(arguments)