import matplotlib
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
scatter_default_args = {'alpha': 0.8, 's': 50}
