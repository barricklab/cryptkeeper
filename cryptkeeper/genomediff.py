"""
This file contains code used to generate and present differences between a reference genome and a mutated genome.
"""
from collections import namedtuple
import pandas as pd
import re
from rich import print as rprint
from Bio import SeqIO
from .cryptkeeper import cryptkeeper
from copy import deepcopy
from .crystalization import build_unified_feature_list

mutation_data = namedtuple('mutation_data', 'type start end details')

def adjust_for_mutations(mutations, results):
    """Receives a list of mutations and a list of results and
    adjusts the results to account for INDELs that shift the start 
    and end positions of a given feature prediction"""

    for mutation in mutations:
        if mutation.type in ["INS", "DEL"]:
            pass

    pass

def prepare_genomes(reference, mutations):
    """Receives a reference genome and a list of mutations and
    returns a genome with the mutations applied. Also subset the original genome."""

    # Find the flanking mutations and subset the reference genome with 1000bp flanking space
    reference = reference.seq
    flanking_mutations = [None, None]
    for mutation in mutations.iterrows():
        mutation_data = {'population': mutation[1]['population'], 'type': mutation[1]['type'], 'start': mutation[1]['start'], 'end': mutation[1]['end'], 'details': mutation[1]['details']}
        if not flanking_mutations[0] and not flanking_mutations[1]:
            flanking_mutations[0] = mutation_data['start']
            flanking_mutations[1] = mutation_data['end']
            continue
        if mutation_data['start'] < flanking_mutations[0]:
            flanking_mutations[0] = mutation_data['start']
        if mutation_data['end'] > flanking_mutations[1]:
            flanking_mutations[1] = mutation_data['end']

    # Expand the window by 1000bp on each side or to the end, whichever is smaller
    if flanking_mutations[0] - 1000 < 0:
        flanking_mutations[0] = 0
    else:
        flanking_mutations[0] -= 1000

    if flanking_mutations[1] + 1000 > len(reference):
        flanking_mutations[1] = len(reference)
    else:
        flanking_mutations[1] += 1000

    # Make this 0-indexed exclusive
    flanking_mutations[0] -= 1

    # Subset the reference genome
    reference = reference[flanking_mutations[0]:flanking_mutations[1]]
    mutated_genome = deepcopy(reference)

    # Sort the mutations from 3' to 5'
    mutations.sort_values(by=['end'], inplace=True, ascending=False)

    oritinal_length = len(reference.seq)

    for mutation in mutations.iterrows():
        mutation_data = {'population': mutation[1]['population'], 'type': mutation[1]['type'], 'start': mutation[1]['start'], 'end': mutation[1]['end'], 'details': mutation[1]['details']}

        # Adjust the start and end positions of mutations to reflect the subset
        offset = flanking_mutations[0]
        mutation_data['start'] -= (offset+1) # subtract +1 for 0 indexing
        mutation_data['end'] -= (offset)     # no +1 for exclusive

        # Apply the mutation

        if mutation_data['type'] == "INS":
            # Assert that details[0] is the same as the mutated region
            original_extend = len(mutation_data['details'][0])-1
            assert mutation_data['details'][0] == mutated_genome[mutation_data['start']-original_extend:mutation_data['end']], f"Region {mutated_genome[mutation_data['start']:mutation_data['end']]} for INS does not match expected initial {mutation_data['details'][0]}"
            mutated_genome = mutated_genome[:mutation_data['start']-original_extend] + mutation_data['details'][1] + mutated_genome[mutation_data['start']+1:]
        elif mutation_data['type'] == "SNP":
            assert mutation_data['details'][0] == mutated_genome[mutation_data['start']], f"Region {mutated_genome[mutation_data['start']]} for SNP does not match expected initial {mutation_data['details']}"
            mutated_genome = mutated_genome[:mutation_data['start']] + mutation_data['details'][1] + mutated_genome[mutation_data['start']+1:]
        elif mutation_data['type'] == "DEL":
            # We can't check with an assertion because that info is not provided by the mutation file
            mutated_genome = mutated_genome[:mutation_data['start']] + mutated_genome[mutation_data['end']:]
        
    return reference, mutated_genome


def load_mutation_file(mutation_file: str) -> pd.DataFrame:
    """Receives a mutation file and returns a list of mutations"""
    raw_data = pd.read_csv(mutation_file, sep=",", header=0, index_col=0)
    processed_data = pd.DataFrame(columns=['population', 'type', 'start', 'end', 'details'])
    for index, row in raw_data.iterrows():
        population = f"{row['treatment']}_{row['population']}_{row['clone']}"
        mutationtype = row['type']
        start = row['start_position']
        end = row['end_position']
        html_mutation = row['html_mutation']
        details = []
        if mutationtype == "INS":
            # Example: (T)<sub>7&rarr;8</sub>
            details.extend(re.search('\((.*)\)<sub>', html_mutation).group(1))
            details.extend(re.search('<sub>(.*)&rarr;', html_mutation).group(1))
            details.extend(re.search('&rarr;(.*)</sub>', html_mutation).group(1))
            details = [details[0]*int(details[1]), details[0]*int(details[2])]
        elif mutationtype == "DEL":
            # Example: &Delta;2&nbsp;bp
            details.extend(re.search('&Delta;([0-9]*)&nbsp;bp', html_mutation).group(1))
        elif mutationtype == "SNP":
            # Example: C&rarr;T
            details.append(re.search('(.*)&rarr;', html_mutation).group(1))
            details.append(re.search('&rarr;(.*)', html_mutation).group(1))
        elif mutationtype == "CON":
            raise NotImplementedError("CON mutations are not yet supported")
        elif mutationtype == "INV":
            raise NotImplementedError("INV mutations are not yet supported")
        elif mutationtype == "SUB":
            raise NotImplementedError("SUB mutations are not yet supported")
        elif mutationtype == "AMP":
            raise NotImplementedError("AMP mutations are not yet supported")
        elif mutationtype == "MOB":  # Unlikely this will ever be supported
            raise NotImplementedError("MOB mutations are not yet supported")
        else:
            raise ValueError(f"Mutation type {mutationtype} not recognized")
        processed_data.loc[index] = [population, mutationtype, int(start), int(end), details]
    return processed_data
    

if __name__ == "__main__":
    mutationdata = "/home/croots/Documents/Code/cryptkeeper/examples/ltee/LTEE-Ecoli-data 2022-08-19 14 23 14.csv"
    mutations = load_mutation_file(mutationdata)
    rprint(mutations)
    reference_genome = "/home/croots/Documents/Code/cryptkeeper/examples/ltee/REL606.gbk"
    reference_genome = SeqIO.read(open(reference_genome, "r"), "genbank")
    reference_genome, mutated_genome = prepare_genomes(reference_genome, mutations)
    reference_results = cryptkeeper(input_file, output=None, circular=False, plot_only=False, name=None, rbs_score_cutoff=2.0, web=False)
    mutant_results = cryptkeeper(input_file, output=None, circular=False, plot_only=False, name=None, rbs_score_cutoff=2.0, web=False)




