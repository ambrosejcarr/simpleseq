__author__ = 'ambrose'

import os
import ftplib
import gzip
import bz2
from subprocess import Popen, PIPE, call, check_output, CalledProcessError
from shutil import rmtree, copyfileobj
from collections import defaultdict


# download links for supported genomes on GEO
_download_links = dict(hg38={
    'genome': ('ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/'
               'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-80/gtf/homo_sapiens/'
            'Homo_sapiens.GRCh38.80.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/'
              'ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz'),
             ('ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/cdna/'
              'Homo_sapiens.GRCh38.cdna.all.fa.gz')]
}, mm38={
    'genome': ('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/'
               'Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/'
            'Mus_musculus.GRCm38.76.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/ncrna/'
              'Mus_musculus.GRCm38.ncrna.fa.gz'),
             ('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/cdna/'
              'Mus_musculus.GRCm38.cdna.all.fa.gz')]
}, ci2={
    'genome': ('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_intestinalis/dna/'
               'Ciona_intestinalis.KH.dna.toplevel.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-81/gtf/ciona_intestinalis/'
            'Ciona_intestinalis.KH.81.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_intestinalis/cdna/'
              'Ciona_intestinalis.KH.cdna.all.fa.gz')]
}, cs2={
    'genome': ('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_savignyi/dna/'
               'Ciona_savignyi.CSAV2.0.dna.toplevel.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-81/gtf/ciona_savignyi/'
            'Ciona_savignyi.CSAV2.0.81.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_savignyi/cdna/'
              'Ciona_savignyi.CSAV2.0.cdna.all.fa.gz')]
})


# define expected file names: efn[organism][filetype] = [names]  (list)
def _define_file_names(download_links):
    file_names = defaultdict(dict)
    for organism, file_types in download_links.items():
        for file_type, link in file_types.items():
            if isinstance(link, list):  # deal with lists of cDNA files
                file_name = [l.split('/')[-1] for l in link]
            else:
                file_name = link.split('/')[-1]
            file_names[organism][file_type] = file_name
    return dict(file_names)

_file_names = _define_file_names(_download_links)


def _check_type(obj, type_):
    if not isinstance(obj, type_):
        raise TypeError('%s must be of type %s' % (obj, type_))


class STAR:

    @staticmethod
    def verify_organism(organism):
        valid_organisms = ['mm38', 'hg38', 'mmhg38', 'cs2', 'ci2']
        if organism:
            if not all(org in valid_organisms for org in organism):
                raise ValueError('Invalid organism value. Supported organisms: %r' %
                                 valid_organisms)

    @staticmethod
    def verify_temp_dir(directory):
        if not os.path.isdir(directory):
            os.makedirs(directory)
        if not directory.endswith('/'):
            directory += '/'
        return directory

    @staticmethod
    def default_alignment_args(fastq_records, n_threads, index, temp_dir):
        default_align_args = {
            '--runMode': 'alignReads',
            '--runThreadN': str(n_threads),
            '--genomeDir': index,
            '--outFilterType': 'BySJout',
            '--outFilterMultimapNmax': '10',  # require <= 10 matches
            '--limitOutSJcollapsed': '2000000',  # deal with many splice variants
            '--alignSJDBoverhangMin': '8',
            '--outFilterMismatchNoverLmax': '0.04',
            '--alignIntronMin': '20',
            '--alignIntronMax': '1000000',
            '--readFilesIn': fastq_records,
            '--outSAMunmapped': 'Within',
            '--outSAMprimaryFlag': 'AllBestScore',  # all equal-scoring reads are primary
            '--outFileNamePrefix': temp_dir,
        }
        return default_align_args

    @classmethod
    def _append_phiX_to_fasta(cls, fasta, cdna=False):

        # download phiX genome
        genome_link = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
                       'Enterobacteria_phage_phiX174_sensu_lato_uid14015/NC_001422.fna')
        phix_genome = cls._download_ftp_file(genome_link, './')

        # add it to the fasta file
        with open(phix_genome, 'rb') as fin:
            data = fin.read()

            if cdna:  # change phix header to reflect our fake transcript id!
                data = data.decode()
                data = data.split('\n')
                data[0] = ('>ENSPXT00000000001 chromosome:NC_001422.1:1:5386:1 '
                           'gene:ENSPXG00000000001 gene_biotype:genome '
                           'transcript_biotype:genome')
                data = ('\n'.join(data)).encode()

            if fasta.endswith('.gz'):
                with gzip.open(fasta, 'ab') as fout:
                    fout.write(data)
            elif fasta.endswith('.bz2'):
                raise ValueError('appending to bz2 files is not currently supported')
            else:
                with open(fasta, 'ab') as fout:
                    fout.write(data)

        os.remove(phix_genome)

    @staticmethod
    def _append_phiX_to_gtf(gtf):
        # fake up some gtf records for a gene and transcript
        gtf_data = '\n'.join([
            '\t'.join([
                'NC_001422.1', 'RefSeq', 'gene', '1', '5386', '.', '+', '.',
                'ID "id0"; Dbxref "taxon:374840"; Is_circular "true"; gbkey "Src"; '
                'genome "genomic"; mol_type "genomic DNA"; nat-host "Escherichia coli"; '
                'gene_id "ENSPXG00000000001";'
            ]), '\t'.join([
                'NC_001422.1', 'RefSeq', 'transcript', '1', '5386', '.', '+', '.',
                'ID "id0"; Dbxref "taxon:374840"; Is_circular "true"; gbkey "Src"; '
                'genome "genomic"; mol_type "genomic DNA"; nat-host "Escherichia coli"; '
                'transcript_id "ENSPXT00000000001"; gene_id "ENSPXG00000000001";'
            ])
        ])

        if gtf.endswith('.gz'):
            with gzip.open(gtf, 'ab') as fout:
                fout.write(gtf_data.encode())
        elif gtf.endswith('.bz2'):
            raise ValueError('appending to bz2 files is not currently supported')
        else:
            with open(gtf, 'ab') as fout:
                fout.write(gtf_data.encode())

    @staticmethod
    def _download_ftp_file(link, prefix, clobber=False):
        """downloads ftp_file available at 'link' into the 'prefix' directory"""

        # check link validity
        if not link.startswith('ftp://'):
            raise ValueError(
                'link must start with "ftp://". Provided link is not valid: %s' % link)

        # create prefix directory if it does not exist
        if not os.path.isdir(prefix):
            os.makedirs(prefix)

        # make sure prefix has a trailing '/':
        if not prefix.endswith('/'):
            prefix += '/'

        ip, *path, file_name = link.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
        path = '/'.join(path)

        # check if file already exists
        if os.path.isfile(prefix + file_name):
            if not clobber:
                return  # file is already present, nothing to do.

        ftp = ftplib.FTP(ip)
        try:
            ftp.login()
            ftp.cwd(path)
            with open(prefix + file_name, 'wb') as fout:
                ftp.retrbinary('RETR %s' % file_name, fout.write)
        finally:
            ftp.close()

        return prefix + file_name

    @staticmethod
    def _gunzip_file(gzipped_file):
        # todo want to replace with gzip module, but it is too slow!
        # with gzip.open(gzipped_file, 'rb') as f:
        #     fout = gzipped_file.replace('.gz', '')
        #     data = f.read()
        #     with open(fout, 'wb') as fo:
        #         fo.write(data)
        unzipped = gzipped_file.replace('.gz', '')
        call(['gunzip -c %s > %s' % (gzipped_file, unzipped)], shell=True)

    @staticmethod
    def _gzip_files(*file_names):
        # todo want to replace with gzip module; but it is too slow!
        for f in file_names:
            call(['gzip', f])

    @staticmethod
    def _merge_files(fout, *file_names):
        with open(fout, 'wb') as wfd:
            for f in file_names:
                if f.endswith('.gz'):
                    fd = gzip.open(f, 'rb')
                elif f.endswith('.bz2'):
                    fd = bz2.open(f, 'rb')
                else:
                    fd = open(f, 'rb')
                copyfileobj(fd, wfd, 1024**2 * 128)

    @classmethod
    def align(cls, fastq_file, index, n_threads, temp_dir, reverse_fastq_file=None,
              **kwargs):

        # make the directory if it doesn't exist
        if not os.path.isdir(temp_dir):
            os.makedirs(temp_dir, exist_ok=True)

        # check if file exists; if it does, return the filename
        if os.path.isfile(temp_dir + 'Aligned.out.sam'):
            if os.path.getsize(temp_dir + 'Aligned.out.sam') > 0:
                return temp_dir + 'Aligned.out.sam'

        runtime_args = cls.default_alignment_args(
            fastq_file, n_threads, index, temp_dir)

        for k, v in kwargs.items():  # overwrite or add any arguments passed from cmdline
            if not isinstance(k, str):
                try:
                    k = str(k)
                except ValueError:
                    raise ValueError('arguments passed to STAR must be strings')
            if not isinstance(v, str):
                try:
                    v = str(v)
                except ValueError:
                    raise ValueError('arguments passed to STAR must be strings')
            runtime_args[k] = v

        # construct command line arguments for STAR
        cmd = ['STAR']
        if reverse_fastq_file:
            for key, value in runtime_args.items():
                if key == '--readFilesIn':
                    cmd.extend((key, value))
                    cmd.append(reverse_fastq_file)
                else:
                    cmd.extend((key, value))
        else:
            for pair in runtime_args.items():
                cmd.extend(pair)

        aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
        out, err = aln.communicate()
        if err:
            raise ChildProcessError(err)

        return temp_dir + 'Aligned.out.sam'

    @classmethod
    def align_multiple_files(cls, fastq_files, index, n_threads, working_dir,
                             reverse_fastq_files=None, **kwargs):

        # use shared memory to map each of the individual cells
        kwargs['--genomeLoad'] = 'LoadAndKeep'

        # make sure forward and reverse fastq file lists match in length
        if reverse_fastq_files:
            if not len(fastq_files) == len(reverse_fastq_files):
                raise ValueError('unequal number of forward and reverse fastq files '
                                 'provided.')

        # make temporary directories for each file
        runs = []
        samfiles = []
        for fastq_file in fastq_files:
            alignment_dir = (working_dir +
                             fastq_file.split('/')[-1].replace('.fastq', '/'))

            # only add the file to runs if the samfile does not already exist
            if os.path.isfile(alignment_dir + 'Aligned.out.sam'):
                if os.path.getsize(alignment_dir + 'Aligned.out.sam') > 0:
                    samfiles.append(alignment_dir + 'Aligned.out.sam')
                    continue

            # file doesn't exist, add it to runs
            try:
                os.mkdir(alignment_dir)
            except FileExistsError:
                pass
            runs.append(cls.default_alignment_args(
                fastq_file, n_threads, index, alignment_dir)
            )
            samfiles.append(alignment_dir + 'Aligned.out.sam')

        # get output filenames

        for runtime_args in runs:
            for k, v in kwargs.items():  # overwrite or add any arguments from cmdline
                if not isinstance(k, str):
                    try:
                        k = str(k)
                    except ValueError:
                        raise ValueError('arguments passed to STAR must be strings')
                if not isinstance(v, str):
                    try:
                        v = str(v)
                    except ValueError:
                        raise ValueError('arguments passed to STAR must be strings')
                runtime_args[k] = v

        # construct command line arguments for STAR runs
        cmds = []
        for i, runtime_args in enumerate(runs):
            cmd = ['STAR']
            if reverse_fastq_files:
                for key, value in runtime_args.items():
                    if key == '--readFilesIn':
                        cmd.extend((key, value))
                        cmd.append(reverse_fastq_files[i])
                    else:
                        cmd.extend((key, value))
            else:
                for pair in runtime_args.items():
                    cmd.extend(pair)
            cmds.append(cmd)

        # load shared memory
        cls.load_index(index)
        try:
            for cmd in cmds:
                aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
                out, err = aln.communicate()
                if err:
                    raise ChildProcessError(err.decode())
        finally:
            cls.remove_index(index)

        return samfiles

    @staticmethod
    def load_index(index):
        # set shared memory; this may bug out with larger indices; test!
        try:
            _ = check_output(['sysctl', '-w', 'kernel.shmmax=36301783210'])
            _ = check_output(['sysctl', '-w', 'kernel.shmall=36301783210'])
        except CalledProcessError:
            pass  # not available on OS X
        star_args = ['STAR', '--genomeDir', index, '--genomeLoad',
                     'LoadAndExit']
        star = Popen(star_args, stderr=PIPE, stdout=PIPE)
        return star.communicate()

    @staticmethod
    def remove_index(index):
        # remove index
        star_args = ['STAR', '--genomeDir', index, '--genomeLoad',
                     'Remove']
        star = Popen(star_args, stderr=PIPE, stdout=PIPE)
        out, err = star.communicate()
        if err:  # don't remove temp files, they have logging info
            return err

        # remove temporary files
        call(['rm', '-r', '_STARtmp/'])
        call(['rm', './Aligned.out.sam', './Log.out', './Log.progress.out'])
        return

    @staticmethod
    def clean_up(directory):
        rmtree(directory)


def alignment_metadata(log_final_out, meta):
    """store a summary of the alignment run written in log_final_out"""

    with open(log_final_out, 'r') as f:
        lines = f.readlines()

        # alignment rates
        meta['mmap_rate'] = float(lines[24].strip().split('\t')[-1][:-1])
        meta['uniq_rate'] = float(lines[9].strip().split('\t')[-1][:-1])
        unmapped_rate = (float(lines[28].strip().split('\t')[-1][:-1]) +
                         float(lines[29].strip().split('\t')[-1][:-1]) +
                         float(lines[30].strip().split('\t')[-1][:-1]))
        meta['unmapped_rate'] = unmapped_rate

        # alignment numbers
        total_reads = int(lines[5].strip().split('\t')[-1])
        meta['total_reads'] = total_reads
        meta['unique_reads'] = int(lines[8].strip().split('\t')[-1])
        meta['mmapped_reads'] = int(lines[23].strip().split('\t')[-1])
        meta['unmapped_reads'] = round(unmapped_rate * total_reads / 100)

        # error rates:
        meta['mismatch_rate'] = float(lines[17].strip().split('\t')[-1][:-1])
        meta['deletion_rate'] = float(lines[18].strip().split('\t')[-1][:-1])
        meta['insertion_rate'] = float(lines[20].strip().split('\t')[-1][:-1])

        # error magnitudes:
        meta['deletion_size'] = float(lines[19].strip().split('\t')[-1])
        meta['insertion_size'] = float(lines[21].strip().split('\t')[-1])

    return meta
