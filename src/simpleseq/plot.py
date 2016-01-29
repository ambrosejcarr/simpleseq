import pickle
import scipy.sparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns

# set a bunch of defaults
sns.set_style('ticks')

# set default fonts
matplotlib.rc('font', **{'family': 'serif',
                         'serif': ['Computer Modern Roman'],
                         'monospace': ['Computer Modern Typewriter']
                         })

matplotlib.rc('figure', **{'figsize': (4, 4),
                           'dpi': 150})

matplotlib.rc('patch', **{'facecolor': 'royalblue',
                          'edgecolor': 'none'})

matplotlib.rc('lines', **{'color': 'royalblue'})

cmap = matplotlib.cm.viridis


def qualitative_colors(n):
    return sns.color_palette('husl', n)

# plotting defaults
defargs = {'alpha': 0.8, 's': 50}


def get_fig(fig=None, ax=None):
    """fills in any missing axis or figure with the currently active one"""
    if not fig:
        fig = plt.gcf()
    if not ax:
        ax = plt.gca()
    return fig, ax


class SparseCounts():

    def __init__(self, coo, genes, cells):
        """

        args:
        -----
        genes: dictionary mapping integer column ids to read names
        cells: dictionary mapping integer row cell ids to cell barcode values
        coo: Coordinate sparse matrix
        """
        csr = coo.tocsr()
        csr[csr < 0] = 256
        coo = csr.tocoo()

        self.coo = coo
        self.genes = genes
        self.cells = cells


    @classmethod
    def from_pickle(cls, filename):
        with open(filename, 'rb') as f:
            return cls(**pickle.load(f))

    def to_pickle(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(vars(self), f)

    def sample_cellsums(self, n):
        colsums = np.ravel(self.coo.tocsr().sum(axis=1))
        total = np.sum(colsums)
        if n > total:
            print('warning: {!s} exceeds requested sample number: {!s}. Upsampling.'
                  ''.format(n, total))
        return np.random.multinomial(n)

    def cell_counts_histogram(self, fig=None, ax=None, title='', log=True, **kwargs):
        fig, ax = get_fig(fig, ax)
        if log:
            ax.set_xlabel('log10(molecules)')
        else:
            ax.set_xlabel('molecules')
        ax.set_ylabel('cells')
        ax.set_title(title)

        cellsums = np.ravel(self.coo.tocsr().sum(axis=1))
        plt.hist(cellsums, **kwargs)
