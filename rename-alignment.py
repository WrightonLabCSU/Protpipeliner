import sys

def load_original_names(input_fasta):
    """Load original names from the input FASTA file and create a mapping to their sequences."""
    original_names = {}
    current_name = ""
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_name = line.strip()[1:]  # remove the '>'
                original_names[current_name] = ""
            else:
                original_names[current_name] += line.strip()
    return original_names

def create_name_mapping(original_names):
    """Create a mapping from renamed sequences (g_0, g_1, ...) to original names."""
    name_mapping = {}
    for index, original_name in enumerate(original_names.keys()):
        renamed_name = f"g_{index}"
        name_mapping[renamed_name] = original_name
    return name_mapping

def replace_names(alignment_file, name_mapping, output_file):
    """Replace renamed sequences in the alignment file with original names."""
    with open(alignment_file, 'r') as f:
        lines = f.readlines()

    with open(output_file, 'w') as f:
        for line in lines:
            if line.startswith('>'):
                renamed_name = line.strip()[1:]  # remove the '>'
                original_name = name_mapping.get(renamed_name, renamed_name)
                f.write(f'>{original_name}\n')
            else:
                f.write(line)

def main(input_fasta, alignment_file, output_file):
    original_names = load_original_names(input_fasta)
    name_mapping = create_name_mapping(original_names)
    replace_names(alignment_file, name_mapping, output_file)
    
    print(f"Renamed sequences in {alignment_file} and saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_sequences.py <input_fasta> <alignment_file> <output_file>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    alignment_file = sys.argv[2]
    output_file = sys.argv[3]

    main(input_fasta, alignment_file, output_file)

