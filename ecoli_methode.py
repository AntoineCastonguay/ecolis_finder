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
from ete3 import Tree, TreeStyle
from ete3.parser.newick import NewickError
from glob import glob
import pysam
import warnings
import pandas as pd
import io


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_input(my_input):
        error_message_list = ['Please provide a folder as input.',
                              'The input folder provided does not exist.',
                              'Make sure files in input folder end with {}'.format(Methods.accepted_extensions)]

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
            first_character = split(first_header)[0]
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
    def check_version(log_file):
        # Not being used right now because versions are captured in the requirements.txt file
        with open(log_file, 'w') as f:
            # Python
            p = subprocess.Popen(['python', '--version'])
            stderr, stdout = p.communicate()

            # Porechop
            p = subprocess.Popen(['porechop', '--version'])
            stderr, stdout = p.communicate()

            # bbmap suite
            p = subprocess.Popen(['bbduk.sh', '--version'])
            stderr, stdout = p.communicate()

            # Filtlong
            p = subprocess.Popen(['filtlong', '--version'])
            stderr, stdout = p.communicate()

            # SAMtools
            p = subprocess.Popen(['samtools', '--version'])
            stderr, stdout = p.communicate()

            # Minimap2
            p = subprocess.Popen(['minimap2', '--version'])
            stderr, stdout = p.communicate()

    @staticmethod
    def test1():
        print('test1 win!')

    @staticmethod
    def test2():
        print('test1 win!')

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder):
        sample_dict = dict()

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].replace('_pass', '').replace('_filtered', '')
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]
                    sample_dict[sample] = file_path
        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

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
    def get_fastq_from_bam(sample, bam_file, fastq_file, output_folder):
        read_dict = dict()
        with pysam.AlignmentFile(bam_file, 'rb') as f:
            for read in f.fetch():
                read_name = read.qname
                read_dict[read_name] = ''  # push reads into dict to avoid duplicates
                # if read_name not in read_dict:
                #     read_dict[read_name] = ''  # push reads into dict to avoid duplicates

        extracted_fastq = output_folder + sample + '.fastq.gz'

        # Parse fastq file to dictionary
        line_counter = 0
        with gzip.open(extracted_fastq, 'wb') as out_fh:
            with gzip.open(fastq_file, 'rt') if fastq_file.endswith('.gz') else open(fastq_file, 'r') as f:
                fastq_entry = list()
                for line in f:
                    line_counter += 1
                    fastq_entry.append(line)
                    if line_counter == 4:  # last line of fastq entry

                        # Ditch the leading "@" and everything after 1st space
                        seq_id = fastq_entry[0].split()[0][1:]

                        # Write to
                        if seq_id in read_dict:
                            out_fh.write(''.join(fastq_entry).encode('ascii'))

                        # Prepare for new fastq entry
                        line_counter = 0
                        fastq_entry = list()

    @staticmethod
    def gzipped_file_size(gzipped_file):
        with gzip.open(gzipped_file, 'rb') as f:
            return f.seek(0, whence=2)

    @staticmethod
    def run_minimap2(sample, ref, fastq_file, cpu, output_folder, keep_bam):
        print('\t{}'.format(sample))

        output_bam = output_folder + sample + '.bam'

        minimap2_cmd = ['minimap2',
                        '-a',
                        '-x', 'map-ont',
                        '-t', str(cpu),
                        '--MD',
                        '--secondary=no',
                        ref,
                        fastq_file]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4', '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '--reference', ref,
                             '-']
        samtools_markdup_cmd = ['samtools', 'markdup',
                                '-r',
                                '-@', str(cpu),
                                '-',
                                output_bam]
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index',
                              output_bam]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_markdup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p4.communicate()

        # Index bam file
        if os.path.exists(output_bam):
            if os.stat(output_bam).st_size != 0:  # bam file exists and not empty
                subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

                # Convert bam to fastq
                Methods.get_fastq_from_bam(sample, output_bam, fastq_file, output_folder)

                # Remove bam
                if not keep_bam:
                    bam_list = glob(output_folder + '*.bam*')
                    for bam in bam_list:
                        os.remove(bam)
        else:
            warnings.warn('No reads were extracted for {}!'.format(sample))

    @staticmethod
    def run_minimap2_parallel(output_folder, ref, sample_dict, cpu, parallel, keep_bam):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref, path, int(cpu / parallel), output_folder, keep_bam)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_minimap2(*x), args):
                pass