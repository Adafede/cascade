# Install from https://github.com/sorenwacker/ms-peakonly

import argparse
import pandas
import os

from ms_peakonly import PeakOnly
from glob import glob


# Get arguments 
parser = argparse.ArgumentParser(
	description='Finds peaks in a directory of mzML(s) file(s)'
	)
parser.add_argument(
	'--input',
	'-i',
	type=str,
	nargs=1,
	help='the directory containing the mzML file(s)'
	)
parser.add_argument(
	'--output',
	'-o',
	type=str,
	nargs=1,
	help='the path to write the output'
	)
args = parser.parse_args()

# Get a list of file names to process (I believe onl mzML files are supported)
fns = glob(args.input[0]+'/*.mzML')

# Instantiate the engine if the neural network weights
# are not already downloaded this will
# also download the models. 
po = PeakOnly(
	model_dir='https://github.com/sorenwacker/ms-peakonly/tree/main/ms_peakonly/models'
	)

# Simply pass the list of filenames to the `process` method.
table = po.process(fns)

# Export
os.makedirs(os.path.dirname(args.output[0]), exist_ok=True)
table.to_csv(args.output[0])