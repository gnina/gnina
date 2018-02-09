#!/usr/bin/python

import argparse
import numpy as np
import itertools
import matplotlib.pyplot as plt
from plumbum import local
import os, gzip, sys, math
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD

plt.style.use('seaborn-white')

def get_time(out):
    '''
    Returns loop time for gnina run; does this by parsing the stdout to get the
    value computed by the internal loop timer. This might be a crappy way of
    doing things (e.g. what if we decide to disable timer printing?) so
    probably this should be changed to use /usr/bin/time or a python timing
    utility
    '''
    tstring = ''
    
    for char in out:
        if char == "\n":
            if tstring.startswith('Loop time'):
                return float(tstring.split()[-1])
                break
            else:
                tstring = ''
        else:
            tstring = tstring + char

def get_torsions(ligand):
    '''
    Returns a list of the number of torsions for all the ligands in an input
    file. File must be in sdf or gzipped sdf format for now. 
    '''
    torsions = []
    ext = os.path.splitext(ligand)
    if ext[-1] == '.sdf':
        f = open(ligand)
    elif ext[-1] == '.sdf.gz':
        f = gzip.open(ligand)
    else:
        print 'Ligands not provided in sdf format\n'
        sys.exit()
    suppl = Chem.ForwardSDMolSupplier(f)
    for mol in suppl:
        if mol is None: continue
        torsions.append(rdMD.CalcNumRotatableBonds(mol))
    f.close()
    return torsions

def run_command(cmd, *args):
    '''
    Runs a command and returns stdout. Includes somewhat convoluted existence
    checking that appears to be necessary (at least for plumbum), unfortunately. 
    '''
    try:
        local.which(cmd)
    except:
        try:
            local.path(cmd).exists()
        except:
            '%s not found\n' %cmd
            sys.exit()
        else:
            if not os.path.dirname(cmd):
                cmd = './' + cmd
    return local[cmd][[arg for arg in args]]()

parser = argparse.ArgumentParser(description='Generates gnina speedup figure')

parser.add_argument('-i', '--input', required=True, nargs='+', help='Input \
binaries to test')

parser.add_argument('-c', '--cpu', required=False, help='CPU binary to test \
against; if not provided, the CPU functionality of the first \
input binary will be used, unless paired_test is passed')

parser.add_argument('-p', '--paired_test', default='False',
        action='store_true', help='Calculate speedup for each binary compared \
with its own CPU implementation; less efficient but ensures that \
changes in CPU code across different versions do not affect speedup \
figure')

parser.add_argument('-r', '--receptor', required=False, nargs='+',
        help='Optionally provide one or more receptors for testing performance')

parser.add_argument('-l', '--ligands', required=False, nargs='+',
        help='Optionally provide one or more ligands for testing performance')

parser.add_argument('-a', '--args', required=False, nargs='+', help='Optionally \
        pass additional arguments to gnina')

parser.add_argument('-o', '--outprefix', default='', help='Optionally provide \
prefix for output plots; default is binary names connected with underscores')

args = parser.parse_args()

if args.receptor and not args.ligands:
    print 'When providing a receptor, must also provide at least one ligand\n'
    sys.exit()
if args.ligands and not args.receptor:
    print 'When providing ligands, must also provide at least one receptor\n'
    sys.exit()

if not args.ligands:
    args.ligands = []
    #TODO: this should be 16, but the new GPU version has a problem that I
    #assume is related to the buffer being too small
    for i in range(16):
        args.ligands.append('performance/a17/test/a17_%dtorsions.sdf' %i)
    assert not args.receptor
    args.receptor = 'performance/a17/a17_rec.pdb'

out_filebase = {}
#generate basenames for plot labels
for bin in args.input:
    if args.args:
        args_ext = '_' + '_'.join([arg.lstrip('-') for arg in args.args]) 
    else:
        args_ext = ''
    out_filebase[bin] = os.path.basename(bin).replace('.', '_') + args_ext

#set plot outprefix if not user-defined
if not args.outprefix:
    out_plotname = '_'.join(out_filebase.values())

#args.cpu takes precedence over args.paired_test
if args.cpu and args.paired_test:
    args.paired_test = False
    print 'CPU argument passed, which takes precedence over paired testing\n'

#iterate over input binaries to generate results.
#there is a single plot showing speedup as a function of the mean number of torsions
#for all the binaries; the default files have the same number of torsions
#for all ligands within a given file
torsions_fig,torsions_ax = plt.subplots()
#there is a grid of plots showing speedup as a function of torsion distribution
#for different workloads; speedup is used to determine the color of the density
#plot
total_plots = len(args.input)
grid_width = int(math.ceil(math.sqrt(total_plots)))
grid_length = int(math.ceil(float(total_plots)/grid_width))
density_fig, density_ax = plt.subplots()
#subplots for the grid plot, used to set axis labels and colors
subplot = {}
#and then compute mean/stderr of speedup
torsions_speedup = {}
torsions_err = {}
density_workloads = {}

for idx,bin in enumerate(args.input):
    #generate plot for speedup as function of ntorsions; for each input file
    #plot vs. mean number of torsions for ligands in file. for the default
    #torsions dataset, all the ligands in the file actually have the same 
    #number of torsions
    ntorsions = []
    density_workloads[idx] = []
    nruns = 3
    #run nruns times to get mean and stderr for speedup, make nruns cmdline arg?
    for run in range(nruns):
        for lignum,ligand in enumerate(args.ligands):
            torsions = get_torsions(ligand)
            density_workloads[idx].append(torsions)
            if idx==0 and run==0:
                ntorsions.append(sum(torsions)/float(len(torsions)))
            if run==0 and lignum==0: cpu_time = [[0] * 3 for
                    i in range(len(args.ligands))]
            if run==0 and lignum==0: gpu_time = [[0] * 3 for
                    i in range(len(args.ligands))]
            if args.cpu and idx == 0:
                out = run_command(args.cpu, '-r', args.receptor, '-l', ligand, '--minimize')
                cpu_time[lignum][run] = get_time(out)
            elif args.paired_test or not args.paired_test and idx == 0:
                out = run_command(bin, '-r', args.receptor, '-l', ligand,
                        '--minimize')
                cpu_time[lignum][run] = get_time(out)
            out = run_command(bin, '-r', args.receptor, '-l', ligand, '--minimize',
            '--gpu', '--cpu', '3')
            gpu_time[lignum][run] = get_time(out)

    torsions_speedup[bin] = []
    torsions_err[bin] = []
    for lignum in range(len(args.ligands)):
        pairs = itertools.product(cpu_time[lignum], gpu_time[lignum])
        vals = [float(i)/j for i,j in pairs]
        torsions_speedup[bin].append(sum(vals)/float(len(vals)))
        torsions_err[bin].append(np.std(vals, ddof=1)/math.sqrt(len(vals)))

    torsions_ax.errorbar(ntorsions, torsions_speedup[bin],
            yerr=torsions_err[bin], marker='o', alpha=0.5, ls='None',
            label=out_filebase[bin])

torsions_ax.set_xlabel('Number of torsions')
torsions_ax.set_ylabel('Speedup over CPU')
torsions_ax.legend(loc='best', bbox_to_anchor=(1,0.5))
torsions_fig.savefig('%s_torsions.pdf' %out_plotname, bbox_inches='tight')

#set up palettes for density plots TODO: make distance in color space
#proportional to relative speedup
# palette = sns.color_palette("RdBu", n_colors=10)
# for plot_num,workloads in density_workloads.items():
    # for trial,workload in enumerate(workloads):
        # subplot[idx] = sns.kdeplot(workload, ax=density_ax[plot_num/
            # grid_width, plot_num % grid_width],
            # color=round(density_speedups[plot_num][trial]),
            # label=out_filebase[plot_num])
    # if int(idx) / grid_width == grid_length-1:
        # subplot[idx].set_xlabel('Torsions Distribution for Workload')
    # if idx % grid_width == 0:
        # subplot[idx].set_ylabel('Speedup Over CPU')
# 
# density_ax.legend(loc='best', bbox_to_anchor=(1,0.5))
# density_fig.savefig('%s_densities.pdf' %out_plotname,bbox_inches='tight')
