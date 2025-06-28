# Title: 02.csv2loom.py
# Author: ydgenomics
# Date: 20250525
# Image: GRN-SCENIC-database--01 /opt/conda/bin/python
# Description: This script is the preparation of runing pySCENIC.
# Pre-requiry:
# params: rds_path

import os,sys
import loompy as lp
import numpy as np
import scanpy as sc
import argparse

# Using argparse to handle command line arguments
parser = argparse.ArgumentParser(description="csv to loom")
parser.add_argument("-i", "--input_csv", type=str, default="scenic.data.csv", help="Input CSV file (default: scenic.data.csv)")
parser.add_argument("-o", "--output_loom", type=str, default="scenic.loom", help="Output loom file (default: scenic.loom)")
args = parser.parse_args()

x=sc.read_csv(args.input_csv);
row_attrs={"Gene":np.array(x.var_names),};
col_attrs={"CellID":np.array(x.obs_names)};
lp.create(args.output_loom,x.X.transpose(),row_attrs,col_attrs);