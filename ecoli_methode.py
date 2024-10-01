import subprocess
import os
import sys
from concurrent import futures
import pathlib
from shutil import move
from psutil import virtual_memory
from multiprocessing import cpu_count
import gzip
from itertools import groupby
from glob import glob
import pysam
import warnings


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_input(my_input):
        if not os.path.exists(my_input):
            raise Exception('Please select an existing file or folder as input.')

        # Check if folder
        if os.path.isdir(my_input):
            file_list = os.listdir(my_input)  # List content of folder
        else:  # elif os.path.isfile(my_input):
            file_list = [my_input]

        # if folder is not empty and all files have the accepted extensions
        if not all([f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]):
            raise Exception('Make sure files in input folder end with {}'.format(Methods.accepted_extensions))

    @staticmethod
    def check_ref(ref):
        if not os.path.isfile(ref):
            raise Exception('The reference file provided does not exist')

        with gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref, 'r') as f:
            first_header = f.readline()
            first_character = first_header[0]
            if first_character != '>':
                raise Exception('The reference file provided does not appear to be a valid fasta file.')

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = total_cpu
            sys.stderr.write("Number of samples to parallel process was set to {}".format(total_cpu))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
    
    @staticmethod
    def find_sam_files(folder):
        sam_files = []
        for filename in os.listdir(folder):
            if filename.endswith('.sam'):
                sam_files.append(os.path.join(folder, filename))
        return sam_files

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def list_files_in_folder(folder, extension):
        return glob(folder + '/*' + extension)

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass
        
    @staticmethod
    def minimap2(genome, primer, output):
        print('minimap2 processing...')

        # Crée le dossier de sortie si nécessaire
        os.makedirs(output, exist_ok=True)

        # Prépare la commande sans redirection
        #minimap_cmd = ['minimap2', '-ax', 'map-ont', genome, primer]

        #BWA_index_cmd = ['bwa', 'index', genome]
        #subprocess.run(BWA_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        #BWA_cmd = ['bwa', 'mem', genome, primer]
        #with open(f'{output}/output.sam', 'w') as outfile, open(os.devnull, 'w') as errfile:
        #    subprocess.run(BWA_cmd, stdout=outfile, stderr=errfile)

        # Commande pour indexer le génome
        index_cmd = ['bowtie2-build', genome, 'genome_index']
        subprocess.run(index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Commande pour aligner les reads au format FASTA
        bowtie_cmd = ['bowtie2', '-x', 'genome_index', '-f', primer, '-k', '1', '--no-mixed', '--no-discordant']
        with open(f'{output}/output.sam', 'w') as outfile, open(os.devnull, 'w') as errfile:
            subprocess.run(bowtie_cmd, stdout=outfile, stderr=errfile)

        # Redirection de la sortie standard vers un fichier et la sortie d'erreur vers /dev/null
        #with open(f'{output}/output.sam', 'w') as outfile, open(os.devnull, 'w') as errfile:
        #    subprocess.run(minimap_cmd, stdout=outfile, stderr=errfile)

    @staticmethod
    def extract_primer_positions(sam_file):
        print('Extrat position primer...')
        primer_positions = {}

        with open(sam_file, 'r') as file:
            for line in file:
                # Ignorer les lignes d'en-tête
                if line.startswith('@'):
                    continue

                # Séparer la ligne en colonnes
                columns = line.strip().split('\t')
                read_id = columns[0]  # ID de la lecture
                flag = int(columns[1])  # Flag de la lecture
                reference_name = columns[2]  # Nom de la séquence de référence
                position = columns[3]  # Position d'alignement

                # Vérifier si la lecture est alignée (flag != 4)
                if flag != 4:
                    # Ajouter l'ID de lecture et la position au dictionnaire
                    primer_positions[read_id] = int(position)
                else:
                    primer_positions[read_id] = '*'

        return primer_positions
    
    @staticmethod
    def write_result(data,output):
        print('Creation of result file...')

        Methods.make_folder(output)

        grouped_data = {}
        for key, value in data.items():
            primer = key[:-2]  # Supprime '_l' ou '_r' pour obtenir le primer
            if primer not in grouped_data:
                grouped_data[primer] = []
            grouped_data[primer].append(value)
        with open(f'{output}/output.txt', 'w') as f:
            for primer, positions in grouped_data.items():
                f.write(f"{primer}\t{'\t'.join(map(str, positions))}\n")