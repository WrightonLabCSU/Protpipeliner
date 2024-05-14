#!/usr/bin/python

# Python script
# Written by: Richard Wolfe
# Original version by Brian C. Thomas
# Modified by Reed Woyda on May 10 2022

# NOTE - This version limits model selection to those compatible with pplacer.

import sys
import argparse
import os

# Create an argument parser object
parser = argparse.ArgumentParser(description='A script to run protpipeliner')

# Add the available arguments
parser.add_argument('-i', '--input_file', type=argparse.FileType('r'), help='fasta file', required=True)
parser.add_argument('-t', '--threads', type=int, help='number of threads', required=True)
parser.add_argument('-b', '--bootstraps', type=int, help='number of bootstraps for RAxML', required=True)
parser.add_argument('-m', '--mode', help='either high (n) or med (h) or low (a) or none (NO_GBLOCKS)', required=True)
parser.add_argument('-a', '--aligned_file', help='Is input file aligned fasta file either T or F', required=True)
parser.add_argument('--osc', help='Running on OSC', action="store_true")
parser.add_argument('--stop_after_prottest', help='Stop script after running prottest', action="store_true")
parser.add_argument('--skip_prottest', help='skip ProtTest and the model to use', default="USE_PROTEST")

# Get the args
args = parser.parse_args()

# Additional argument tests
if args.threads <= 0:
    print("Error: argument --threads <= 0")
    sys.exit(0)

if args.bootstraps <= 0:
    print("Error: argument --bootstraps <= 0")
    sys.exit(0)

mode = ""
if args.mode == "high":
    mode = "n"
elif args.mode == "med":
    mode = "h"
elif args.mode == "low":
    mode = "a"
elif args.mode == "none":
    mode = "NO_GBLOCKS"
else:
    print("Error -m needs to be high or med or low or none")
    sys.exit(1)

if args.aligned_file != "T":
    if args.aligned_file != "F":
        print("Error -a needs to be T or F")
        sys.exit(1)

# DELETE_FILES = "TRUE"
DELETE_FILES = "FALSE"

##############################################
# Set paths to programs on OSC -- use string that runs command on the system
muscle_path = "muscle"
fasta2phy_path = "/opt/scripts/bin/fasta_fastq/fasta2phy"
phyToFasta_path = "/opt/scripts/bin/Phylogeny_Protpipe/phyToFasta.pl"
gblocks_path = "Gblocks"
prottest_path = "/opt/prottest-3.4-20140123/prottest-3.4.jar"
raxml_path = "raxmlHPC-PTHREADS"

if args.osc:  # If running on OSC need to change paths
    muscle_path = "/users/PAS1018/osu9681/bin/protpipeliner/muscle"
    fasta2phy_path = "/users/PAS1018/osu9681/bin/protpipeliner/fasta2phy"
    phyToFasta_path = "/users/PAS1018/osu9681/bin/protpipeliner/phyToFasta.pl"
    gblocks_path = "/users/PAS1018/osu9681/bin/protpipeliner/Gblocks"
    prottest_path = "/users/PAS1018/osu9681/bin/protpipeliner/prottest-3.4-20140123/prottest-3.4.jar"
    raxml_path = "/users/PAS1018/osu9681/bin/protpipeliner/raxmlHPC-PTHREADS"

print("Script started ...")

# Get the filename from the -i path
file_name = os.path.basename(args.input_file.name)
print(file_name)

# Close the files so there are no problems
args.input_file.close()

# Check if the RAxML files with this name already exist
# If they do then RAxML will be in error
# If they exist then we will exit and let the user delete the files
if os.path.isfile("RAxML_bestTree." + file_name):
    print("ERROR...RAxML_bestTree." + file_name + " file already exists")
    sys.exit(1)
if os.path.isfile("RAxML_bipartitionsBranchLabels." + file_name):
    print("ERROR...RAxML_bipartitionsBranchLabels." + file_name + " file already exists")
    sys.exit(1)
if os.path.isfile("RAxML_bipartitions." + file_name):
    print("ERROR...RAxML_bipartitions." + file_name + " file already exists")
    sys.exit(1)
if os.path.isfile("RAxML_bootstrap." + file_name):
    print("ERROR...RAxML_bootstrap." + file_name + " file already exists")
    sys.exit(1)
if os.path.isfile("RAxML_info." + file_name):
    print("ERROR...RAxML_info." + file_name + " file already exists")
    sys.exit(1)

# Rename the fasta sequences
print("\n\n-- Renaming input fasta sequences")

with open(args.input_file.name, 'r') as infile:
    outfilename = file_name + '.rename'
    with open(outfilename, 'w') as outfile:
        scaffold = 0
        for line in infile:
            if line.startswith('>'):
                line = ">g_" + str(scaffold) + "\n"
                scaffold += 1
            outfile.write(line)

# Running muscle which aligns the sequences
if args.aligned_file == 'F':
    print("-- Starting muscle MSA on " + file_name + ".rename")

    cmd = muscle_path + ' -in ' + file_name + '.rename' + ' -out ' + file_name + '.al'
    retvalue = os.system(cmd)  # Returns 0 if no error in command

    if retvalue != 0:  # If command failed
        print("ERROR command did not return 0 and failed")
        sys.exit(1)  # Return 0 if success

    print("-- done with muscle")

if args.aligned_file == 'T':
    print("Input file is aligned fasta so copy input file to input file .al")
    cmd = 'cp ' + file_name + '.rename ' + file_name + '.al'
    retvalue = os.system(cmd)  # Returns 0 if no error in command

    if retvalue != 0:  # If command failed
        print("ERROR command did not return 0 and failed")
        sys.exit(1)  # Return 0 if success

# Convert to phylip format
print("-- Converting to phylip format and cleaning up")

cmd = "cp " + file_name + ".al " + file_name + ".tmp"
retvalue = os.system(cmd)  # Returns 0 if no error in command
if retvalue != 0:  # If command failed
    print("ERROR command did not return 0 and failed")
    sys.exit(1)  # Return 0 if success

cmd = 'grep -c ">" ' + file_name + '.al'
f = os.popen(cmd)
seqs = f.read().strip()  # Remove whitespace and endlines
print("The number of sequences in alignment file = ", seqs)

if mode != "NO_GBLOCKS":
    # Running Gblocks
    print("-- Running Gblocks to detect conserved regions in the MSA")

    cmd = gblocks_path + ' ' + file_name + '.tmp' + ' -t=p -p=n -e=.fst -b1=$(((' + str(seqs) + '/2)+1)) -b2=$(((' + str(seqs) + '/2)+1)) -b3=$(((' + str(seqs) + '/2))) -b4=2 -b5=' + mode

    print("cmd = ")
    print(cmd)

    retvalue = os.system(cmd)  # Returns 0 if no error in command
    print("-- done with Gblocks")

else:
    # Not running Gblocks so rename the alignment file
    cmd = 'mv ' + file_name + '.tmp ' + file_name + '.tmp.fst'
    retvalue = os.system(cmd)  # Returns 0 if no error in command

    if retvalue != 0:  # If command failed
        print("ERROR command did not return 0 and failed")
        sys.exit(1)  # Return 0 if success

# Convert to phy 
cmd = 'cat ' + file_name + '.tmp.fst | ' + fasta2phy_path + ' > ' + file_name + '.phy'
retvalue = os.system(cmd)  # Returns 0 if no error in command

if retvalue != 0:  # If command failed
    print("ERROR command did not return 0 and failed")
    sys.exit(1)  # Return 0 if success

# Run ProtTest
if args.skip_prottest == "USE_PROTEST":
    print("\n\n-- Starting ProtTest")

    cmd = 'java -jar ' + prottest_path + ' -i ' + file_name + '.phy -o ' + file_name + '.model -all-matrices -all-distributions -log disabled -threads ' + str(args.threads)
    f = os.popen(cmd)
    output = f.read()

    # Extract best model
    cmd = 'grep "Best model according to" ' + file_name + """.model | awk '{print$6}' | awk -F+ '{print $1'} | tr "[:lower:]" "[:upper:]" """
    f = os.popen(cmd)
    model = f.read().strip()  # Remove whitespace both ends

    # I added this because I got 2 best results when I tried a sample
    best = model.split()  # Will split on whitespace including endl
    model = best[0]  # Take first element

    # Ensure the model is supported by pplacer
    pplacer_supported_models = ["LG", "WAG", "JTT", "DAYHOFF", "BLOSUM62", "MTREV"]
    if model not in pplacer_supported_models:
        print(f"Model {model} is not supported by pplacer. Selecting a fallback model LG.")
        model = "LG"

    print("\n\nBest model from ProtTest was: " + model)

    if args.stop_after_prottest:
        sys.exit(0)

    # Run RAxML
    print("\n\n--Starting raxmlHPC")
    cmd = raxml_path + ' -f a -m PROTCAT' + model + ' -n ' + file_name + ' -N ' + str(args.bootstraps) + ' -p 1234 -s ' + file_name + '.phy -x 1234 -T ' + str(args.threads)
    f = os.popen(cmd)

    for line in f:
        print(line, end='')
    f.close()

else:
    if args.stop_after_prottest:
        sys.exit(0)

    # Run RAxML
    print("\n\n--Starting raxmlHPC")
    cmd = raxml_path + ' -f a -m ' + args.skip_prottest + ' -n ' + file_name + ' -N ' + str(args.bootstraps) + ' -p 1234 -s ' + file_name + '.phy -x 1234 -T ' + str(args.threads)
    f = os.popen(cmd)

    for line in f:
        print(line, end='')
    f.close()

# Rename the tree with original names
print("\n\nRenaming output tree to use original names")

fn = 'RAxML_bestTree.' + file_name
try:
    with open(fn, 'r') as r:
        raxml_file = r.read()

    with open(args.input_file.name, 'r') as f:
        index = 0
        for line in f:
            if line.startswith('>'):
                gene_label = 'g_' + str(index) + ':'
                line_list = line.split()
                orig_label = line_list[0].strip().replace(">", "").replace(":", "_").replace(")", "_").replace("(", "_").replace(",", "_").replace(";", "_")
                orig_label += ':'
                raxml_file = raxml_file.replace(gene_label, orig_label)
                index += 1

    out_filename = "bestTree." + file_name + "_mode_" + args.mode + ".renamed"
    with open(out_filename, 'w') as out_file:
        out_file.write(raxml_file)

except IOError:
    print("Error.. RAxML did not make file " + fn)
    sys.exit(1)

# Rename the Branch labels tree with original names
print("\n\nRenaming Branch labels output tree to use original names")

fn = 'RAxML_bipartitionsBranchLabels.' + file_name
try:
    with open(fn, 'r') as r:
        raxml_file = r.read()

    with open(args.input_file.name, 'r') as f:
        index = 0
        for line in f:
            if line.startswith('>'):
                gene_label = 'g_' + str(index) + ':'
                line_list = line.split()
                orig_label = line_list[0].strip().replace(">", "").replace(":", "_").replace(")", "_").replace("(", "_").replace(",", "_").replace(";", "_")
                orig_label += ':'
                raxml_file = raxml_file.replace(gene_label, orig_label)
                index += 1

    out_filename = "bipartitionsBranchLabels." + file_name + "_mode_" + args.mode + ".renamed"
    with open(out_filename, 'w') as out_file:
        out_file.write(raxml_file)

except IOError:
    print("Error.. RAxML did not make file " + fn)

# Rename the RAxML_bipartitions tree with original names
print("\n\nRenaming Branch labels RAxML_bipartitions output tree to use original names")

fn = 'RAxML_bipartitions.' + file_name
try:
    with open(fn, 'r') as r:
        raxml_file = r.read()

    with open(args.input_file.name, 'r') as f:
        index = 0
        for line in f:
            if line.startswith('>'):
                gene_label = 'g_' + str(index) + ':'
                line_list = line.split()
                orig_label = line_list[0].strip().replace(">", "").replace(":", "_").replace(")", "_").replace("(", "_").replace(",", "_").replace(";", "_")
                orig_label += ':'
                raxml_file = raxml_file.replace(gene_label, orig_label)
                index += 1

    out_filename = "bipartitions." + file_name + "_mode_" + args.mode + ".renamed"
    with open(out_filename, 'w') as out_file:
        out_file.write(raxml_file)

except IOError:
    print("Error.. RAxML did not make file " + fn)

# Remove intermediate files no longer needed
print("--removing unneeded files")

if DELETE_FILES == "TRUE":
    os.system('rm -f ' + file_name + '.rename')
    os.system('rm -f ' + file_name + '.al')
    os.system('rm ' + file_name + '.tmp')
    os.system('rm ' + file_name + '.tmp.fst')
    os.system('rm ' + file_name + '.phy')
    os.system('rm ' + file_name + '.phy.reduced')
    os.system('rm ' + file_name + '.model')
    os.system('rm ' + 'RAxML_bestTree.' + file_name)
    os.system('rm ' + 'RAxML_bipartitionsBranchLabels.' + file_name)
    os.system('rm ' + 'RAxML_bipartitions.' + file_name)
    os.system('rm ' + 'RAxML_bootstrap.' + file_name)
    os.system('rm ' + 'RAxML_info.' + file_name)
    os.system('rm -rf snapshot')

print("\n\nScript finished...")

