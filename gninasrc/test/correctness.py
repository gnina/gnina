#!/usr/bin/python

import argparse
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from plumbum import local
from plumbum.cmd import sdsorter, awk

plt.style.use('seaborn-white')

parser = argparse.ArgumentParser(description='Checks correctness and makes \
performance plots for gnina')
parser.add_argument('-i', '--input', required=False, nargs='+', help='Input \
binaries to test')
parser.add_argument('-c', '--cpu_off', dest='cpu_on', action='store_false',
        help='Include CPU in performance plots')
parser.add_argument('-t', '--threshold', type=float, default=.8,
        help='Correlation threshold before we start to get worried.')
parser.add_argument('-r', '--receptor', default='4pps.pdb', help=
'Receptor to use for testing')
parser.add_argument('-l', '--ligand', default='4pps_ligands.sdf', help=
'Ligands to use for testing')
parser.add_argument('-s', '--structures', required=False, nargs='+',
help='Provide structure files produced as output from running the program \
instead of binaries. If these are provided, the program will not be run; \
first file will be used as independent variable.')
parser.add_argument('--old', dest='old_format', action='store_true', help='Pre-CNN scoring format for \
extraction from sdsorter output')
parser.set_defaults(cpu_on=True, old_format=False)
args = parser.parse_args()

if args.old_format:
    column = 3
else: 
    column = 4
# Check correctness; generate plot if under the threshold
# TODO: print more specific info about how many ligands differed and the
# expected value of the deviation
fig,ax = plt.subplots()
if not args.structures:
    cpu_cmd = local[os.getcwd() + '/%s' % args.input[0]]
    cpu = cpu_cmd["-r", args.receptor, "-l", args.ligand, "-o",
            args.input[0].split(".")[0] + "_cpu_out.sdf", "--minimize"]
    cpu()
    get_cpu_out = sdsorter['-print','-omit-header',args.input[0].split(".")[0] +
    "_cpu_out.sdf"] | awk["{print $%s}" % column]
    cpu_out = np.array([float(x) for x in get_cpu_out().split()])
    
    for file in args.input:
        cmd = local[os.getcwd() + '/%s' % file]
        gpu = cmd["-r", args.receptor, "-l", args.ligand, "-o",
                file.split(".")[0] + "_gpu_out.sdf", "--minimize", "--gpu"]
        gpu()
        get_gpu_out = sdsorter['-print','-omit-header',file.split(".")[0] +
        "_gpu_out.sdf"] | awk["{print $%s}" % column]
        gpu_out = np.array([float(x) for x in get_gpu_out().split()])
        ax.plot(cpu_out, gpu_out, marker='o', alpha=0.5, ls='None', label='%s GPU, \
        R=%.2f' % (file, pearsonr(cpu_out, gpu_out)[0]))
    ax.set_xlabel('CPU-calculated energy')
    ax.set_ylabel('GPU-calculated energy')
else:
    outputs = {}
    get_cpu_out = sdsorter['-print','-omit-header',args.structures[0]] |\
    awk["{print $%s}" % column]
    cpu_out = np.array([float(x) for x in get_cpu_out().split()])
    outputs[args.structures[0].split('.')[0]] = cpu_out
    for file in args.structures[1:]:
        get_gpu_out = sdsorter['-print','-omit-header',file] |\
        awk["{print $%s}" % column]
        gpu_out = np.array([float(x) for x in get_gpu_out().split()])
        outputs[file.split('.')[0]] = gpu_out
    outname = ''
    for item in outputs.keys():
        if item != args.structures[0].split('.')[0]:
            outname += '%s, ' % item
        best = -1.0
        best_id = ''
        worst = 1.0
        worst_id = ''
        for other_item in outputs.keys():
            if item != other_item:
                correlation = pearsonr(outputs[other_item], outputs[item])[0]
                if correlation > best:
                    best = correlation
                    best_id = other_item
                if correlation < worst:
                    worst = correlation
                    worst_id = other_item
        ax.plot(cpu_out, gpu_out, marker='o', alpha=0.5, ls='None', label='%s, \
        best R=%.2f with %s, worst R=%.2f with %s' % (item, best, best_id,
            worst, worst_id))
    ax.set_xlabel(args.structures[0].split('.')[0])
    ax.set_ylabel('Additional smina runs')
    # ax.set_ylabel(outname)

ax.legend(loc='best', bbox_to_anchor=(1,0.5))
fig.savefig('correctness.pdf', bbox_inches='tight')
