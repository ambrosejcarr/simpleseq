import fileinput
from itertools import islice
from collections.abc import Iterable, Mapping
import seqc_clean
import numpy as np


class FastqRecord:

    __slots__ = ['_data']

    def __init__(self, record: [bytes, bytes, bytes, bytes]):
        self._data = record

    @property
    def name(self) -> bytes:
        return self._data[0]

    @name.setter
    def name(self, value: bytes):
        self._data[0] = value

    @property
    def sequence(self) -> bytes:
        return self._data[1]  # todo this was the error .. ?

    @sequence.setter
    def sequence(self, value: bytes):
        self._data[1] = value

    @property
    def name2(self) -> bytes:
        return self._data[2]

    @name2.setter
    def name2(self, value: bytes):
        self._data[2] = value

    @property
    def quality(self) -> bytes:
        return self._data[3]

    @quality.setter
    def quality(self, value: bytes):
        self._data[3] = value

    def __bytes__(self) -> bytes:
        return b''.join(self._data)

    def __str__(self) -> str:
        return bytes(self).decode()

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def annotations(self) -> list:
        """
        returns:
        --------
        list of annotations present in the fastq header
        """
        try:
            end = self.name.index(b';')
            return self.name[:end].split(b':')
        except ValueError:
            return []

    @property
    def metadata(self) -> dict:
        """
        returns:
        --------
        dictionary of annotations and fields, if any are present"""
        try:
            start = self.name.rindex(b'|')
        except ValueError:
            return {}
        fields = {}
        for field in self.name[start + 1:].split(b':'):
            k, v = field.split(b'=')
            fields[k] = v
        return fields

    def add_annotation(self, values: Iterable) -> None:
        """prepends a list of annotations to the name field of self.name"""
        self._data[0] = self.name[:-1] + b';' + b':'.join(values) + b'\n'


    def add_metadata(self, values: Mapping) -> None:
        """appends a list of metadata fields to the name field of self.name"""
        self.name += b'|' + b':'.join(k + '=' + v for k, v in values.items())

    def average_quality(self) -> int:
        """"""
        return np.mean(np.frombuffer(self.quality, dtype=np.int8, count=len(self)))\
            .astype(int) - 33


class FastqReader(seqc_clean.file_io.Reader):

    def __iter__(self):
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:
            while True:
                yield FastqRecord(list(islice(f, 4)))


def merge_fastq(merge_function, fout, genomic, barcode=None):
    """
    annotate genomic fastq with barcode information from reverse read

    args:
    -----
    merge_function: function from merge_functions.py
    fout: merged output file name
    genomic: fastq containing genomic data
    barcode: fastq containing barcode data
    """
    genomic = FastqReader(genomic)
    barcode = FastqReader(barcode)
    with open(fout, 'wb') as f:
        for g, b in zip(genomic, barcode):
            r = merge_function(g, b)
            f.write(bytes(r))

    return fout
