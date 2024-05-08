# Protpipeliner
This script combines MUSCLE, GBLOCKS, Prottest and RaxML for an efficient and streamlined method of building trees.

Example usage:
`protpipeliner.py -i dmso_fortree_aligned.fasta -t 10 -b 100 -m low -a T`

Example further stepd for creating a DRAM2 `.refpkg` for Trees:

1) Determine the number of positions:
```
seqkit stats dmso_refs.fasta.al
#Output:
#file                format  type     num_seqs  sum_len  min_len  avg_len  max_len
#dmso_refs.fasta.al  FASTA   Protein        86  160,046    1,861    1,861    1,861
```

2) Create the `.refpkg`
```
taxit create -P dmso_package -l 1861 -p dmso_refs.fasta.model -f RAxML_info.dmso_refs.fasta --stats-type RAxML
```

The main protpipeliner.py script does not rename the alignment file. Use `rename-alignment.py` to rename this. This is needed for DRAM2 trees `.refpkg`:

Example usage:

`python rename-alignment.py dmso_refs.fasta dmso_refs.fasta.al renamed_dmso_alignment.al`
