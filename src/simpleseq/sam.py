__author__ = 'ambrose'

import os
import ftplib
import gzip
import bz2
from subprocess import Popen, PIPE, call, check_output, CalledProcessError
from shutil import rmtree, copyfileobj
import numpy as np
import simpleseq
import fileinput
from itertools import islice
from collections import Counter, defaultdict
from scipy.sparse import coo_matrix


class SamAnnotation:
    """Container for annotations added to fastq records by simpleseq.fastq.merge()"""

    __slots__ = ['_fields']

    def __init__(self, name_field):
        fields, *name = name_field.split(b';')
        self._fields = fields.split(b':') + [b';'.join(name)]

    def __repr__(self) -> str:
        return '<SamAnnotation: {0}>'.format(str(self))

    def __str__(self) -> str:
        return bytes(self).decode()

    def __bytes__(self) -> bytes:
        return b':'.join(self._fields)[:-1] + b';' + self._fields[-1]

    @property
    def pool(self) -> bytes:
        return self._fields[0]

    @property
    def cell(self) -> bytes:
        return self._fields[1]

    @property
    def rmt(self) -> bytes:
        return self._fields[2]

    @property
    def n_poly_t(self) -> int:
        try:
            return int(self._fields[3])
        except ValueError:  # empty field denotes zero
            return 0

    # @property
    # def barcode_quality(self) -> int:
    #     return int(self._fields[4])

    @property
    def qname(self) -> bytes:
        return self._fields[4]

    @property
    def encoded_pool(self) -> int:
        return simpleseq.encodings.DNA3Bit.encode(self._fields[0])

    @property
    def encoded_cell(self) -> int:
        return simpleseq.encodings.DNA3Bit.encode(self._fields[1])

    @property
    def encoded_rmt(self) -> int:
        return simpleseq.encodings.DNA3Bit.encode(self._fields[2])

    def nearest_TTS(self, annotation):
        raise NotImplementedError  # todo implement


class SamRecord:
    """Record object for iterating over samfiles"""

    __slots__ = ['_data']

    def __init__(self, record: bytes):
        self._data = record.strip(b'\n').split(b'\t')

    def __repr__(self) -> str:
        return '<SamRecord: {0}>'.format(str(self))

    def __str__(self) -> str:
        return b'\t'.join(self._data).decode()

    def __bytes__(self) -> bytes:
        return b'\t'.join(self._data)

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def qname(self) -> bytes:
        return self._data[0]

    @property
    def flag(self) -> int:
        return int(self._data[1])

    @property
    def rname(self) -> bytes:
        return self._data[2]

    @property
    def pos(self) -> int:
        return int(self._data[3])

    @property
    def mapq(self) -> int:
        return int(self._data[4])

    @property
    def cigar(self) -> bytes:
        return self._data[5]

    @property
    def rnext(self) -> bytes:
        return self._data[6]

    @property
    def pnext(self) -> int:
        return int(self._data[7])

    @property
    def tlen(self) -> int:
        return int(self._data[8])

    @property
    def seq(self) -> bytes:
        return self._data[9]

    @property
    def qual(self) -> bytes:
        return self._data[10]

    @property
    def optional_fields(self) -> dict:
        flags_ = {}
        for f in self._data[11:]:  # are all optional fields
            k, type_, v = f.split(b':')
            if type_ == b'i':
                flags_[k] = int(v)
            else:
                flags_[k] = v
        return flags_

    @property
    def strand(self) -> bytes:
        return b'-' if self.flag & 16 else b'+'

    @property
    def annotations(self) -> SamAnnotation:
        return SamAnnotation(self.qname)

    @property
    def is_mapped(self) -> bool:
        return False if self.flag & 4 else True

    @property
    def is_unmapped(self) -> bool:
        return True if self.flag & 4 else False

    @property
    def is_multimapped(self) -> bool:
        return True if self.optional_fields[b'NH'] > 1 and self.is_mapped else False

    @property
    def is_uniquely_mapped(self) -> bool:
        return True if self.optional_fields[b'NH'] == 1 and self.is_mapped else False

    @property
    def average_quality(self) -> int:
        """"""
        return np.mean(np.frombuffer(self.qual, dtype=np.int8, count=len(self)))\
            .astype(int) - 33

    @property
    def data(self) -> [bytes]:
        return self._data

    @property
    def dust_score(self) -> int:
        """"""  # todo
        counts = {}
        for i in range(len(self.seq) - 2):
            kmer = self.seq[i:i + 3]
            counts[kmer] = counts.get(kmer, 0) + 1

        # Calculate dust score
        score = np.sum([i * (i - 1) / 2 for i in counts.values()]) / (len(self.seq) - 3)

        # Scale score (Max score possible is no. of 3mers/2)
        score = np.int8(score / ((len(self.seq) - 2) / 2) * 100)

        return score


class MultiAlignment:

    __slots__ = ['_data']

    def __init__(self, records):
        """initialized from a list of bytes records"""
        self._data = [record.data for record in records]

    @classmethod
    def from_bytes(cls, records):
        """initialized from a list of SamRecords objects"""
        return cls(SamRecord(r) for r in records)

    def __repr__(self) -> str:
        return '<MultiAlignment: {0}>'.format(str(self))

    def __str__(self) -> str:
        return b'\n'.join(b'\t'.join(r) for r in self._data).decode()

    def __bytes__(self) -> bytes:
        return b'\n'.join(b'\t'.join(r) for r in self._data)

    def __len__(self) -> int:
        return len(self._data)

    @property
    def qname(self) -> bytes:
        return self._data[0][0]

    @property
    def flags(self) -> [int]:
        return [int(r[1]) for r in self._data]

    @property
    def rnames(self) -> [bytes]:
        return [r[2] for r in self._data]

    @property
    def positions(self) -> [int]:
        return [int(r[3]) for r in self._data]

    @property
    def mapqs(self) -> [int]:
        return [int(r[4]) for r in self._data]

    @property
    def cigars(self) -> [bytes]:
        return [r[5] for r in self._data]

    @property
    def rnext(self) -> [bytes]:
        return [r[6] for r in self._data]

    @property
    def pnext(self) -> [int]:
        return [int(r[7]) for r in self._data]

    @property
    def tlen(self) -> [int]:
        return [int(r[8]) for r in self._data]

    @property
    def seq(self) -> bytes:
        return self._data[0][9]

    @property
    def qual(self) -> bytes:
        return self._data[0][10]

    @property
    def optional_fields(self) -> [dict]:
        annotations = []
        for record in self._data:
            flags_ = {}
            for f in record[11:]:  # are all optional fields
                k, type_, v = f.split(b':')
                if type_ == b'i':
                    flags_[k] = int(v)
                else:
                    flags_[k] = v
            annotations.append(flags_)
        return annotations

    @property
    def strands(self) -> [bytes]:
        return [b'-' if flag & 16 else b'+' for flag in self.flags]

    @property
    def annotations(self) -> SamAnnotation:
        return SamAnnotation(self.qname)

    @property
    def is_mapped(self) -> bool:
        return False if self.flags[0] & 4 else True

    @property
    def is_unmapped(self) -> bool:
        return True if self.flags[0] & 4 else False

    @property
    def is_multimapped(self) -> bool:
        return True if len(self._data) > 1 else False

    @property
    def is_uniquely_mapped(self) -> bool:
        return True if len(self._data) == 1 and self.is_mapped else False

    def best_alignment(self):
        # todo implement; should filter based on AlignScore & CIGAR & TTS distance
        raise NotImplementedError

    @property
    def average_quality(self) -> int:
        """"""
        return np.mean(np.frombuffer(self.qual, dtype=np.int8, count=len(self)))\
            .astype(int) - 33

    @property
    def dust_score(self) -> int:
        """"""
        counts = {}
        for i in range(len(self.seq) - 2):
            kmer = self.seq[i:i + 3]
            counts[kmer] = counts.get(kmer, 0) + 1

        # Calculate dust score
        score = np.sum([i * (i - 1) / 2 for i in counts.values()]) / (len(self.seq) - 3)

        # Scale score (Max score possible is no. of 3mers/2)
        score = np.int8(score / ((len(self.seq) - 2) / 2) * 100)

        return score

    def to_alignment_row(
            self, genome_annotation: simpleseq.gtf.Annotation) -> (tuple, [int], [int]):

        anno = self.annotations

        cell = anno.encoded_cell
        rmt = anno.encoded_rmt
        poly_t = anno.n_poly_t
        dust_score = self.dust_score
        barcode_quality = anno.barcode_quality
        genomic_quality = self.average_quality
        max_alignment_score = max(f[b'AS'] for f in self.optional_fields)

        if not self.is_mapped or self.is_multimapped:
            features, positions = [], []
        else:
            features = genome_annotation.translate(self.strands[0], self.rnames[0],
                                                   self.positions[0])
            if features:  # todo should be unique
                positions = self.positions[0]
            else:
                positions = []

        row = (cell, rmt, poly_t, dust_score, barcode_quality, genomic_quality,
               max_alignment_score)

        return row, features, positions


class SamReader(simpleseq.reader.Reader):
    """"""

    def __iter__(self):
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:

            # get rid of header lines
            file_iterator = iter(f)
            first_record = next(file_iterator)
            while first_record.startswith(b'@'):
                first_record = next(file_iterator)
                continue
            yield SamRecord(first_record)

            for record in file_iterator:  # now, run to exhaustion
                yield SamRecord(record)

    def iter_multialignments(self):
        """yields lists of all alignments for each fastq record"""
        sam_iter = iter(self)
        multialignment = [next(sam_iter)]
        for record in sam_iter:
            if record.qname == multialignment[0].qname:  # is this expensive?
                multialignment.append(record)
            else:
                yield MultiAlignment(multialignment)
                multialignment = [record]
        yield MultiAlignment(multialignment)

    def iter_chunk(self, size: int=int(1e9)):
        """
        Yield chunks of approximately size, guaranteeing that chunks do not disrupt
        multi-alignments; default size=1GB

        args:
        -----
        size: approximate number of bytes per chunk
        """

        def find_multialignment_boundary(chunk):
            """find the position that the final complete multialignment terminates at"""

            def final_record_start_index(bytestring):
                return bytestring.rindex(b'\n') + 1

            incomplete_record_start = final_record_start_index(chunk)
            last_complete = final_record_start_index(chunk[:incomplete_record_start])
            record = SamRecord(chunk[last_complete:incomplete_record_start])
            while True:

                prev_start = final_record_start_index(chunk[:last_complete])
                compare = SamRecord(chunk[prev_start:last_complete])

                # keep backing up by one record until the name no longer matches.
                if record.qname != compare.qname:
                    return last_complete
                else:
                    record = compare
                    last_complete = prev_start

        # function main
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:
            partial_chunk = b''
            while True:
                record_chunk = f.read(size)
                if not record_chunk:
                    if partial_chunk:
                        yield partial_chunk
                    break
                index = find_multialignment_boundary(record_chunk)
                yield partial_chunk + record_chunk[:index]
                partial_chunk = record_chunk[index:]

    @property
    def empty(self) -> bool:
        try:
            next(iter(self))
            return False
        except StopIteration:
            return True

    def estimate_record_number(self, n_records: int=10000) -> int:
        """
        estimate the number of records in the files passed to self

        args:
        -----
        n_records: the number of records from which to calculate the mean record size

        returns:
        --------
        expected number of records
        """
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:
            data = list(islice(f, n_records))
        mean_record_size = np.mean([len(r) for r in data])
        return int(self.size / mean_record_size)

    def molecule_counts_subset(self, samfile, gtf):

        # read the samfile, get a list of all the cells
        cells = set()
        with open(samfile, 'rb') as f:
            fiter = iter(f)
            line = next(fiter)
            while line.startswith(b'@'):
                line = next(fiter)

            # get cells
            for line in fiter:
                cell = b''.join(line.split(b':')[:2])
                cells.add(cell)

        # assign each cell an id
        cells = dict(zip(cells, range(len(cells))))

        # read the gtf to get all the genes
        genes = set()
        rd = simpleseq.gtf.GTFReader(gtf)
        for gene in rd.iter_genes():
            genes.add(gene.gene_name)

        # create an annotation to read all the genes
        gtf_anno = simpleseq.gtf.Annotation(gtf)
        gtf_anno.create_interval_tree()

        # assign each gene an id
        genes = dict(zip(genes, range(len(genes))))

        n = int(1e8 + 1)

        # construct the arrays to hold the sparse matrix data; assumes < 1 mol / read
        row = np.zeros(n, dtype=np.uint32)
        col = np.zeros(n, dtype=np.uint16)
        molecule_data = np.ones(n, dtype=np.uint16)
        i = 0

        # track all observed molecules
        observed_molecules = set()

        riter = iter(self)
        j = 0
        while i < 1e8:
            record = next(riter)
            if record.is_uniquely_mapped:
                sam_anno = record.annotations
                gene = gtf_anno.translate(record.strand, record.rname, record.pos)
                if not gene:
                    continue
                molecule = (sam_anno.cell, sam_anno.rmt, gene)
                if molecule in observed_molecules:
                    continue
                else:
                    row[i] = genes[gene]
                    col[i] = cells[sam_anno.pool + sam_anno.cell]
                    observed_molecules.add(molecule)
                    i += 1
            j += 1

        # construct molecules
        mrow = row[:i]
        mcol = col[:i]
        mdata = molecule_data[:i]
        molecules = coo_matrix((mdata, (mrow, mcol)), shape=(len(genes), len(cells)))

        return molecules, genes, cells

    def molecule_counts(self, samfile, gtf, alignment_summary):

        cells = set()
        # read the samfile, get a list of all the cells
        with open(samfile, 'rb') as f:
            fiter = iter(f)
            line = next(fiter)
            while line.startswith(b'@'):
                line = next(fiter)

            # get cells
            for line in fiter:
                cell = b''.join(line.split(b':')[:2])
                cells.add(cell)

        # assign each cell an id
        cells = dict(zip(cells, range(len(cells))))

        # read the gtf to get all the genes
        genes = set()
        rd = simpleseq.gtf.GTFReader(gtf)
        for gene in rd.iter_genes():
            genes.add(gene.gene_name)

        # create an annotation to read all the genes
        gtf_anno = simpleseq.gtf.Annotation(gtf)
        gtf_anno.create_interval_tree()

        # assign each gene an id
        genes = dict(zip(genes, range(len(genes))))

        # get the maximum number of reads/molecules
        metadata = get_alignment_metadata(alignment_summary)
        n = metadata['unique_reads']

        # construct the arrays to hold the sparse matrix data; assumes < 1 mol / read
        row = np.zeros(n, dtype=np.uint32)
        col = np.zeros(n, dtype=np.uint16)
        molecule_data = np.ones(n, dtype=np.uint8)
        i = 0

        # track all observed molecules
        observed_molecules = set()

        for record in self:
            if record.is_uniquely_mapped:
                sam_anno = record.annotations
                gene = gtf_anno.translate(record.strand, record.rname, record.pos)
                if not gene:
                    continue
                molecule = (sam_anno.cell, sam_anno.rmt, gene)
                if molecule in observed_molecules:
                    continue
                else:
                    row[i] = genes[gene]
                    col[i] = cells[sam_anno.pool + sam_anno.cell]
                    observed_molecules.add(molecule)
                    i += 1

        # construct molecules
        mrow = row[:i]
        mcol = col[:i]
        mdata = molecule_data[:i]
        molecules = coo_matrix((mdata, (mrow, mcol)), shape=(len(genes), len(cells)))

        return molecules, genes, cells



    def pileup(self, gtf, alignment_summary=None):
        """aggregate data from samfile"""
        unmapped_data = {
            'cell': Counter(),
            'rmt': Counter(),
            'average_quality': Counter(),
            'pool': Counter(),
        }
        mapped_data = {
            'cell': Counter(),
            'rmt': Counter(),
            'average_quality': Counter(),
            'pool': Counter(),
            'unique_position': defaultdict(Counter),  # by chromosome
            'multi_positions': defaultdict(Counter),
            'TTS_dist': Counter()
        }

        for ma in self.iter_multialignments():
            anno = ma.annotations
            unmapped_data['cell'][anno.cell] += 1
            unmapped_data['rmt'][anno.rmt] += 1
            unmapped_data['pool'][anno.pool] += 1
            aq = ma.average_quality
            unmapped_data['average_quality'][aq] += 1
            mapped_data['cell'][anno.cell] += 1
            mapped_data['rmt'][anno.rmt] += 1
            mapped_data['pool'][anno.pool] += 1
            mapped_data['average_quality'][aq] += 1

            # todo implement TTS finding + translation into gene space
            raise NotImplementedError


            # separate unique and multi-positions
            if ma.is_uniquely_mapped:
                mapped_data['unique_position'][ma.rnames[0]][ma.positions[0]] += 1
            else:
                for i in range(len(ma)):
                    mapped_data['multi_positions'][ma.rnames[i]][ma.positions[i]] += 1

        if alignment_summary:
            mapped_data['alignment_summary'] = get_alignment_metadata(alignment_summary)

        return mapped_data, unmapped_data


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


def get_alignment_metadata(log_final_out, meta=None):
    """store a summary of the alignment run written in log_final_out"""

    if meta is None:
        meta = {}

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

