
import pickle
import scipy.sparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
import simpleseq

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
        self.coo = coo
        self.coo.data = self.coo.data.astype(np.uint32)
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
        return np.random.multinomial(n, colsums / np.sum(colsums))

    def cell_counts_histogram(self, fig=None, ax=None, title='', log=True, **kwargs):
        fig, ax = get_fig(fig, ax)
        if log:
            ax.set_xlabel('log10(molecules)')
        else:
            ax.set_xlabel('molecules')
        ax.set_ylabel('cells')
        ax.set_title(title)
        cellsums = np.ravel(self.coo.tocsr().sum(axis=1))
        plt.hist(cellsums, log=log, **kwargs)
        labels = ax.get_xticklabels()
        plt.setp(labels, rotation=90)
        sns.despine(ax=ax)

    def yield_curve(self, threshold, count_vector, fig=None, ax=None, **kwargs):
        """plot yield (cells above threshold) for various numbers of total counts."""
        fig, ax = get_fig(fig, ax)
        passing = [np.sum(self.sample_cellsums(c) > threshold) for c in count_vector]
        ax.plot(count_vector, passing, **kwargs)
        return fig, ax


def compare_yield(count_matrices):
    f = plt.figure(figsize=(6.5, 5))
    ax = f.add_axes([0.1, 0.1, 0.55, 0.55])
    colors = qualitative_colors(len(count_matrices))
    for c, color in zip(count_matrices, colors):
        sc = SparseCounts.from_pickle(c)
        sc.yield_curve([1000], [100000, 500000, 1000000, 2000000, 5000000], c=color,
                       fig=f, ax=ax, label=c.replace('.p', ''))
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.xlabel('total mols.')
    plt.ylabel('cells w. > 1000 mols.')
    plt.title('comparative yield')
    sns.despine()
    plt.tight_layout()
    return f, ax


def compare_yield_sam(samfiles):
    """yield should properly be compared by subsampling the sam files;

    never see more than about 5M molecules; can straight up count molecules?"""
    raise NotImplementedError


def compare_rmt_distributions(rmt_frequencies):
    raise NotImplementedError


class compare_cell_distributions():

    def __init__(self, count_matrices):

        self.count_matrices = count_matrices

        with open(count_matrices[0], 'rb') as f:
            data = pickle.load(f)
            columns = [count_matrices.replace('.p', '')]
            seed = pd.DataFrame(np.ravel(data['coo'].tocsr().sum(axis=1)),
                                     columns=columns)

        for c in count_matrices[1:]:
            with open(c, 'rb') as f:
                data = pickle.load(f)
                columns = [count_matrices.replace('.p', '')]
                right = pd.DataFrame(np.ravel(data['coo'].tocsr().sum(axis=1)),
                                         columns=columns)
                seed = pd.merge(seed, right, left_index=True, right_index=True,
                                how='outer')

        self.cellsums = seed

    def boxplot(self):
        raise NotImplementedError

    def to_pickle(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(vars(self), f)

    @classmethod
    def from_pickle(cls, filename):
        with open(filename, 'rb') as f:
            cls(**pickle.load(f))


def rmt_histogram(rmt_counts, fig=None, ax=None, bins=15, log=True, title='RMT Histogram',
                  **kwargs):
    fig, ax = get_fig(fig, ax)
    ax.set_xlabel('number sequences')
    if log:
        ax.set_ylabel('log(RMTs)')
    else:
        ax.set_ylabel('RMTs')

    ax.set_title(title)

    # delete outliers // failed keys
    keylen = len(next(iter(rmt_counts.keys())))
    del rmt_counts[b'T' * keylen]
    del rmt_counts[b'A' * keylen]

    counts = list(rmt_counts.values())
    bin_counts, bin_edges = np.histogram(counts, bins=bins)
    left = np.arange(len(bin_counts))
    if log:
        bin_counts = np.log(bin_counts)

    ax.bar(left=left, height=bin_counts, width=1)
    ax.set_xticks(np.arange(len(bin_counts) + 1))
    ax.set_xticklabels(bin_edges)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=90)
    sns.despine(ax=ax)

    return fig, ax


def comparative_alignment(alignment_summaries):
    """plot various comparisons of alignment results"""
    plt.figure(figsize=(8, 5.5))
    alndata = []
    for summary in alignment_summaries:
        alndata.append(simpleseq.sam.get_alignment_metadata(summary))

    # below is specific for nested plots
    names = [a.split('/')[0] for a in alignment_summaries]

    # uniq_rate mmap_rate unmapped_rate
    unique = np.array([s['uniq_rate'] for s in alndata])
    multim = np.array([s['mmap_rate'] for s in alndata])

    i = np.argsort(unique)
    i = i[::-1]  # reverse the sort
    left = np.arange(len(unique))

    plt.barh(left, unique[i], 1, color='seagreen', label='unique')
    plt.barh(left, multim[i], 1, unique[i], color='royalblue', label='multi')
    plt.barh(left, np.array([100] * len(left)) - (unique[i] + multim[i]),
             1, unique[i] + multim[i],
             color='indianred', label='unmapped')
    plt.yticks(left + 0.5, np.array(names)[i], fontsize=10)
    plt.xlabel('percentage')
    plt.title('comparative alignment summary')
    plt.ylim((0, len(left)))
    plt.legend(frameon=True)
    sns.despine()
    plt.tight_layout()


def plot_unique_by_date(alignment_summaries, metadata):
    plt.figure(figsize=(8, 5.5))
    df_meta = pd.DataFrame.from_csv(metadata)
    df_meta['Date Produced'] = pd.to_datetime(df_meta['Date Produced'])

    alndata = []
    for summary in alignment_summaries:
        alndata.append(simpleseq.sam.get_alignment_metadata(summary))

    unique = pd.Series(np.array([s['uniq_rate'] for s in alndata]),
                       index=alignment_summaries)

    # plot unique alignments
    index = df_meta.index.intersection(unique.index)
    order = df_meta.loc[index].sort(columns='Date Produced', ascending=False).index
    left = np.arange(len(index))
    height = unique.ix[order]
    width = 0.9
    plt.barh(left, height, width)
    plt.yticks(left + 0.5, order, fontsize=10)
    ymin, ymax = 0, len(left)
    plt.ylim((ymin, ymax))
    plt.xlabel('percentage')
    plt.title('comparative alignment summary')
    plt.ylabel('time (descending)')

    # plot klein in-drop line
    plt.vlines(unique['Klein_in_drop'], ymin, ymax, color='indianred', linestyles='--')

    sns.despine()
    plt.tight_layout()


def plot_unique_by_meta_col(alignment_summaries, metadata, color_col):
    plt.figure(figsize=(8, 5.5))
    df_meta = pd.DataFrame.from_csv(metadata)
    df_meta['Date Produced'] = pd.to_datetime(df_meta['Date Produced'])

    alndata = []
    for summary in alignment_summaries:
        alndata.append(simpleseq.sam.get_alignment_metadata(summary))

    unique = pd.Series(np.array([s['uniq_rate'] for s in alndata]),
                       index=alignment_summaries)

    # order by date produced, get only values that have metadata
    index = df_meta.index.intersection(unique.index)
    order = df_meta.loc[index].sort(columns='Date Produced', ascending=False).index

    # separate by column
    categories = set(df_meta[color_col])
    colors = ['royalblue', 'indianred', 'seagreen', 'orange']
    range_start = 0
    width = 0.9
    ylabels = []
    for c, cat in zip(colors, categories):
        in_cat = df_meta.ix[order].index[df_meta.ix[order][color_col] == cat]
        height = unique.ix[in_cat]
        left = np.arange(range_start, range_start + len(height))
        plt.barh(left, height, width, facecolor=c, label='{}={}'.format(color_col, cat))
        range_start += len(left)
        ylabels.extend(in_cat)
    ymin, ymax = 0, range_start
    plt.title('Unique alignments vs. {}'.format(color_col))
    plt.yticks(np.arange(len(ylabels)) + 0.5, ylabels, fontsize=10)
    plt.ylim((ymin, ymax))
    plt.xlabel('percentage')
    plt.ylabel('time (descending)')
    try:
        plt.legend(frameon=True)
    except:
        pass  # not sure what's causing the problem here.
    # plot klein in-drop line
    # plt.vlines(unique['Klein_in_drop'], ymin, ymax, color='indianred', linestyles='--')

    sns.despine()
    plt.tight_layout()


def unique_alignment_by_metadata(alignment_summaries, metadata, column, color_col):
    """scatterplot of unique alignment rate vs. date published, colored by type"""
    df_meta = pd.DataFrame.from_csv(metadata)
    df_meta['Date Produced'] = pd.to_datetime(df_meta['Date Produced'])


    alndata = []
    for summary in alignment_summaries:
        alndata.append(simpleseq.sam.get_alignment_metadata(summary))

    unique = pd.Series(np.array([s['uniq_rate'] for s in alndata]),
                       index=alignment_summaries)

    # get intersection between unique names and dataframe

    index = df_meta.index.intersection(unique.index)
    plt.scatter(df_meta.ix[index, column],
                unique.ix[index])


