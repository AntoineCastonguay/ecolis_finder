import_os
from_argparse_import_ArgumentParser
from_multiprocessing_import_cpu_count
from_psutil_import_virtual_memory
from_ecoli_methode_import_Methods


__author___=_'castonguaya'
__version___=_'0.1'


class_Ecoli(object):
____def___init__(self,_args):
________#_I/O
________self.input_=_os.path.abspath(args.input)
________self.essentiel_=_os.path.abspath(args.essentiel)
________self.output_folder_=_os.path.abspath(args.output)
________self.ref_genome_=_os.path.abspath(args.genome)

________#_Performance
________self.cpu_=_args.threads
________self.parallel_=_args.parallel
________self.mem_=_args.memory

________self.run()

____def_run(self):
________print('Checking_a_few_things...')
________Methods.check_input(self.input)
________Methods.check_ref(self.ref_genome)

________#_Check_if_number_of_CPU_and_memory_requested_are_valid
________self.cpu,_self.parallel_=_Methods.check_cpus(self.cpu,_self.parallel)
________self.mem_=_Methods.check_mem(self.mem)

________done_extract_=_self.output_folder_+_'/done_extract'
________done_result_=_self.output_folder_+_'/done_result'

________extract_folder_=_os.path.join(self.output_folder,_'1_extract/')
________result_folder_=_os.path.join(self.output_folder,_'2_result/')

________Methods.make_folder(self.output_folder)
________print('\tAll_good!')

________if_not_os.path.exists(done_extract):
____________Methods.alignment(self.ref_genome,_self.input,_extract_folder)
____________Methods.flag_done(done_extract)
________else:
____________print('Skipping_extract._Already_done.')

________sam_file_=_Methods.find_sam_files(extract_folder)

________if_not_os.path.exists(done_result):
____________position_=_Methods.extract_primer_positions(sam_file[0],_self.essentiel)
____________Methods.write_result(position,_result_folder)
____________Methods.flag_done(done_result)
________else:
____________print('Skipping_result._Already_done.')

________print('DONE!')

if___name___==_"__main__":
____max_cpu_=_cpu_count()
____max_mem_=_int(virtual_memory().total_*_0.85_/_1000000000)__#_in_GB

____parser_=_ArgumentParser(description='Extract,_assemble_and_compare_Nanopore_reads_matching_a_reference_sequence.')
____parser.add_argument('-g',_'--genome',_
________________________required=True,_
________________________help='Reference_genome_for_primer_mapping._Mandatory.')
____parser.add_argument('-i',_'--input',_
________________________required=True,_
________________________help='Folder_that_contains_the_fasta_files_or_individual_fasta_file._Mandatory.')
____parser.add_argument('-e',_'--essentiel',_
________________________required=True,_
________________________help='')
____parser.add_argument('-o',_'--output',_
________________________required=True,_
________________________help='Folder_to_hold_the_result_files._Mandatory.')
____parser.add_argument('-p',_'--parallel',_metavar='2',
____________________required=False,
____________________type=int,_default=2,
____________________help='Number_of_samples_to_process_in_parallel._Default_is_2._Optional.')
____parser.add_argument('-t',_'--threads',_metavar=str(max_cpu),
____________________required=False,
____________________type=int,_default=max_cpu,
____________________help='Number_of_threads._Default_is_maximum_available({})._Optional.'.format(max_cpu))
____parser.add_argument('-m',_'--memory',_metavar=str(max_mem),
________________________required=False,
________________________type=int,_default=max_mem,
________________________help='Memory_in_GB._Default_is_85%%_of_total_memory_({})'.format(max_mem))
____parser.add_argument('-v',_'--version',_action='version',_version=f'{os.path.basename(__file__)}:_version_{__version__}')
____
____arguments_=_parser.parse_args()
____Ecoli(arguments)