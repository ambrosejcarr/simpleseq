from setuptools import setup

__author__ = 'Ambrose J. Carr'
__version__ = '0.0.1'

setup(name='simpleseq',
      version=__version__,
      description='simple single-cell sequencing tools',
      author='Ambrose J. Carr',
      author_email='mail@ambrosejcarr.com',
      package_dir={'': 'src'},
      packages=['simpleseq'],
      install_requires=[
          'numpy>=1.10.0',
          'pandas>=0.16.0',
          'matplotlib>=1.5.0',
          'seaborn',
          'regex',
          'boto3',
          'requests'],
      scripts=['src/scripts/merge_annotations',
               'src/scripts/extract_test_fastq',
               'src/scripts/align',
               'src/scripts/align_basespace',
               'src/scripts/align_s3',
               'src/scripts/count_molecules',
               'src/scripts/post_process_mars_seq',
               'src/scripts/rmt_distribution',
               'src/scripts/rmt_distribution_mars']
      )
