#!/usr/bin/env python 
def get_char_at_position(phylip_file, species_name, position):
    with open(phylip_file, 'r') as file:
        # Read the header (first line) to skip it
        file.readline()

        for line in file:
            # Each line in PHYLIP contains a species name followed by the sequence
            if line.startswith(species_name):
                # Remove the species name from the line
                sequence = line[len(species_name):].strip()

                # Check if the sequence is split across lines
                full_sequence = [sequence]
                for seq_line in file:
                    if not seq_line.strip():
                        break
                    full_sequence.append(seq_line.strip())

                # Join all sequence lines into one single sequence
                full_sequence = ''.join(full_sequence)

                # Return the character at the given position
                return full_sequence[position - 100:position+100]  # -1 because positions are 0-based in Python

        return None  # Return None if species not found

# Example usage:
phylip_file = '/media/lepidodactylus/2aa24196-95e9-4ebf-8899-7161cb272356/home/leptodactylus/genomes_2/busco_run_05_results/updated_concatenated_names.phylip'
species_name = 'Podarcis_muralis'
position = 2004674

character = get_char_at_position(phylip_file, species_name, position)
if character:
    print(f"Character at position {position} of {species_name} is: {character}")
else:
    print(f"Species {species_name} not found in the PHYLIP file.")
