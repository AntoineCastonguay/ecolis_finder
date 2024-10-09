import_subprocess
import_os
import_sys
from_concurrent_import_futures
import_pathlib
from_shutil_import_move
from_psutil_import_virtual_memory
from_multiprocessing_import_cpu_count
import_gzip
from_itertools_import_groupby
from_glob_import_glob
import_pysam
import_warnings


class_Methods(object):
____accepted_extensions_=_['.fq',_'.fq.gz',
___________________________'.fastq',_'.fastq.gz',
___________________________'.fasta',_'.fasta.gz',
___________________________'.fa',_'.fa.gz',
___________________________'.fna',_'.fna.gz']

____@staticmethod
____def_check_input(my_input):
________if_not_os.path.exists(my_input):
____________raise_Exception('Please_select_an_existing_file_or_folder_as_input.')

________#_Check_if_folder
________if_os.path.isdir(my_input):
____________file_list_=_os.listdir(my_input)__#_List_content_of_folder
________else:__#_elif_os.path.isfile(my_input):
____________file_list_=_[my_input]

________#_if_folder_is_not_empty_and_all_files_have_the_accepted_extensions
________if_not_all([f.endswith(tuple(Methods.accepted_extensions))_for_f_in_file_list]):
____________raise_Exception('Make_sure_files_in_input_folder_end_with_{}'.format(Methods.accepted_extensions))

____@staticmethod
____def_check_ref(ref):
________if_not_os.path.isfile(ref):
____________raise_Exception('The_reference_file_provided_does_not_exist')

________with_gzip.open(ref,_'rt')_if_ref.endswith('.gz')_else_open(ref,_'r')_as_f:
____________first_header_=_f.readline()
____________first_character_=_first_header[0]
____________if_first_character_!=_'>':
________________raise_Exception('The_reference_file_provided_does_not_appear_to_be_a_valid_fasta_file.')

____@staticmethod
____def_check_cpus(requested_cpu,_n_proc):
________total_cpu_=_cpu_count()

________if_1_>_requested_cpu_>_total_cpu:
____________requested_cpu_=_total_cpu
____________sys.stderr.write("Number_of_threads_was_set_to_{}".format(requested_cpu))
________if_1_>_n_proc_>_total_cpu:
____________n_proc_=_total_cpu
____________sys.stderr.write("Number_of_samples_to_parallel_process_was_set_to_{}".format(total_cpu))

________return_requested_cpu,_n_proc

____@staticmethod
____def_check_mem(requested_mem):
________max_mem_=_int(virtual_memory().total_*_0.85_/_1000000000)__#_in_GB
________if_requested_mem:
____________if_requested_mem_>_max_mem:
________________requested_mem_=_max_mem
________________sys.stderr.write("Requested_memory_was_set_higher_than_available_system_memory_({})".format(max_mem))
________________sys.stderr.write("Memory_was_set_to_{}".format(requested_mem))
________else:
____________requested_mem_=_max_mem

________return_requested_mem

____@staticmethod
____def_make_folder(folder):
________#_Will_create_parent_directories_if_don't_exist_and_will_not_return_error_if_already_exists
________pathlib.Path(folder).mkdir(parents=True,_exist_ok=True)
____
____@staticmethod
____def_find_sam_files(folder):
________sam_files_=_[]
________for_filename_in_os.listdir(folder):
____________if_filename.endswith('.sam'):
________________sam_files.append(os.path.join(folder,_filename))
________return_sam_files

____@staticmethod
____def_list_to_file(my_list,_output_file):
________with_open(output_file,_'wt')_as_f:
____________for_l_in_my_list:
________________f.write('{}\n'.format(l))

____@staticmethod
____def_list_files_in_folder(folder,_extension):
________return_glob(folder_+_'/*'_+_extension)

____@staticmethod
____def_flag_done(flag_file):
________with_open(flag_file,_'w')_as_f:
____________pass
________
____@staticmethod
____def_alignment(genome,_primer,_output):
________print('Alignment_processing...')

________file_=_Methods.list_files_in_folder(primer,_'fa')
________mon_dict_=_{}

________for_f_in_file:
____________name_file_=_os.path.basename(f)
____________name_=_os.path.splitext(name_file)[0]
____________number_=_name[-1]
____________base_=_name[:-1]
____________if_base_not_in_mon_dict:
________________mon_dict[base]_=_{}
____________mon_dict[base][number]_=_f

________#_Crée_le_dossier_de_sortie_si_nécessaire
________os.makedirs(output,_exist_ok=True)

________#_Alignment_BWA
________BWA_index_cmd_=_['bwa',_'index',_genome]
________subprocess.run(BWA_index_cmd,_stdout=subprocess.DEVNULL,_stderr=subprocess.STDOUT)

________for_key,_sub_dict_in_mon_dict.items():
____________BWA_cmd_=_['bwa',_'mem',_genome,_sub_dict['1'],_sub_dict['2']]
____________with_open(f'{output}/BWA_output.sam',_'w')_as_outfile,_open(os.devnull,_'w')_as_errfile:
________________subprocess.run(BWA_cmd,_stdout=outfile,_stderr=errfile)

____@staticmethod
____def_extract_primer_positions(sam_file,_txt_file):
________print('Extrat_position_primer...')
________primer_positions_=_{}
________with_open(txt_file,_'r')_as_fichier:
____________liste_gene_essentiel_=_[ligne.strip()_for_ligne_in_fichier]

________with_open(sam_file,_'r')_as_file:
____________for_line_in_file:
________________#_Ignorer_les_lignes_d'en-tête
________________if_line.startswith('@'):
____________________continue

________________#_Séparer_la_ligne_en_colonnes
________________columns_=_line.strip().split('\t')
________________read_id_=_columns[0]__#_ID_de_la_lecture
________________flag_=_int(columns[1])__#_Flag_de_la_lecture
________________position_=_int(columns[3])__#_Position_d'alignement
________________qualite_=_columns[5]__#qualité_alignment
________________postion_mate_=_int(columns[7])
________________length_=_int(columns[8])

________________if_read_id_in_liste_gene_essentiel:
____________________essentiel_=_True
________________else:
____________________essentiel_=_False

________________list_var_=_[position,qualite,postion_mate,length,essentiel]

________________#_Vérifier_si_la_lecture_est_alignée_(flag_!=_4)
________________#_Ajouter_l'ID_de_lecture_et_la_position_au_dictionnaire
________________if_read_id_not_in_primer_positions:
____________________primer_positions[read_id]_=_{}
________________primer_positions[read_id][flag]_=_list_var

________return_primer_positions
____
____@staticmethod
____def_write_result(data,output):
________print('Creation_of_result_file...')

________Methods.make_folder(output)
________with_open(f'{output}/output.txt',_'w')_as_f:
____________for_read_id,_sub_dict_in_data.items():
________________for_flag,_list_var_in_sub_dict.items():
____________________if_flag_==_99_or_flag_==_147:
________________________f.write(f"{read_id}\t{flag}\t{list_var[0]}\t{list_var[2]}\t{list_var[3]}\t{list_var[1]}\t{list_var[4]}\n")
____________________else:
________________________f.write(f"{read_id}\t{flag}\t{list_var[0]}\t{list_var[2]}\t{list_var[3]}\t{list_var[1]}\t{list_var[4]}\n")