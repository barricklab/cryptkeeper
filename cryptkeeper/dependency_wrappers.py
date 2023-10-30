import numpy as np
import os
from Bio import SeqIO
from collections import namedtuple
import tempfile
import subprocess
from logging import debug, info, warning, error, critical
from operator import itemgetter, attrgetter
from promoter_calculator import promoter_calculator


def ostir(inseq, threads=1):
    #start_codons = ['ATG', 'GTG']
    i=0
    sequence_length = 0
    forward_seq = ''
    reverse_seq = ''

    for this_seq in SeqIO.parse(inseq, "fasta"):
        i += 1
        assert i == 1, "Only one sequence allowed"
        sequence_length = len(this_seq)
        forward_seq = str(this_seq.upper().seq)
        reverse_seq = str(this_seq.reverse_complement().upper().seq)

    #Run OSTIR Calculator twice on entire sequences. Once for each strand.
    from ostir.ostir import run_ostir

    findings_forward = run_ostir(forward_seq, threads=threads, verbosity=0)
    findings_reverse = run_ostir(reverse_seq, threads=threads, verbosity=0)

    for finding in findings_forward:
        finding['strand'] = "+"
    for finding in findings_reverse:
        finding['strand'] = '-'

    result = findings_forward + findings_reverse


    clean_findings = []    # Remove duplicates
    [clean_findings.append(finding) for finding in result if finding not in clean_findings]
    processed_data = []
    rbs_hit = namedtuple('rbs_hit', 'position start_codon strand score score2')
    for finding in clean_findings:
        position = finding['start_position']
        pos_1 = int(position) + 1
        strand = finding['strand']
        position = pos_1 if not strand == "-" else sequence_length - pos_1 + 1
        start_codon = forward_seq[pos_1-2:pos_1+1] if not strand == "-" else reverse_seq[pos_1-2:pos_1+1]
        strand = finding['strand']
        score = finding['expression']
        score2 = finding['dG_total']
        if score > 0:
            processed_data.append(rbs_hit(position, start_codon, strand, float(score), score2))

    return processed_data

def transterm(infile, circular_length):
    with tempfile.TemporaryDirectory() as temp_dir:
        i=0
        for this_seq in SeqIO.parse(infile, "fasta"):
            i += 1
            assert i == 1, "Only one sequence per file is supported"
            main_seq = this_seq.upper()

        #unique step, create dummy cords files so that entire sequence is downstream of genes on both strands
        with open(temp_dir + '.dummy.coords','w') as dummy_coords_file:
            dummy_coords_file.write('gene1 1 2 ' + main_seq.id + '\n')
            dummy_coords_file.write('gene2 ' + str(len(main_seq)) + ' ' + str(len(main_seq)-1) + ' ' + main_seq.id + '\n')

        transterm_expdat_path = ''
        transterm_path = subprocess.check_output(['which', 'transterm']).strip()
        if "TRANSTERM_EXPDAT_PATH" in os.environ:
            transterm_expdat_path = os.getenv("TRANSTERM_EXPDAT_PATH")
        elif transterm_path:
            # Find the expterm.dat file
            transterm_path = subprocess.check_output(['which', 'transterm']).strip()
            transterm_expdat_path = os.path.normpath(transterm_path).decode('utf-8')
            transterm_expdat_path = transterm_expdat_path.split(os.sep)
            transterm_expdat_path = os.sep.join(transterm_expdat_path[:-2]) + os.sep +  "data" + os.sep + 'expterm.dat'

        if not transterm_expdat_path:
            warning("Could not detect transterm installation. Please verify that transterm is installed via conda or that the environment variable TRANSTERM_EXPDAT_PATH is set correctly.")
            return 2

        if isinstance(transterm_expdat_path, bytes):
            transterm_expdat_path = transterm_expdat_path.decode('utf-8')

        #run once for predictions on both strands
        #print('transterm --min-conf=70  -p ' + transterm_expdat_path + ' ' + options.i + options.o + '.dummy.coords > ' + options.o + '.predictions.txt')
        subprocess.call('transterm --min-conf=70  -p ' + transterm_expdat_path + ' ' + infile + ' ' + temp_dir + '.dummy.coords > ' + temp_dir + '.predictions.txt', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)



        #returns a list of dictionaries for the rows
    
        # Parse output and create one summary file
        final_list = _read_transterm_output(temp_dir +'.predictions.txt')

        # Sort list by coordinate
        final_list = sorted(final_list, key=itemgetter('start'))
        transtermhits = namedtuple("transterm_result", "start end strand conf hairpin_score tail_score seq_upstream seq_hairpin_open seq_hairpin_loop seq_hairpin_close seq_tail")
        predictions = []
        for entry in final_list:
            entry['start'] -= circular_length
            entry['end'] -= circular_length
            predictions.append(transtermhits(**entry))

    return predictions

def _read_transterm_output(input_file_name):
    entry_list = []

    last_line_was_prediction = False
    new_entry = {}

    lines = []
    with open(input_file_name) as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()

        # skip all but lines of interest
        if line == "":
            continue
        split_line = line.split()

        # Pick up the sequence line
        if last_line_was_prediction:
            new_entry["seq_upstream"] = split_line[0]
            new_entry["seq_hairpin_open"] = split_line[1]
            new_entry["seq_hairpin_loop"] = split_line[2]
            new_entry["seq_hairpin_close"] = split_line[3]
            new_entry["seq_tail"] = split_line[4]

            entry_list.append(new_entry)

            new_entry = {}
            last_line_was_prediction = False
            continue

        if len(split_line) < 10:
            continue
        if split_line[0] != 'TERM':
            continue

        # New terminator
        last_line_was_prediction = True
        new_entry = {}

        new_entry["strand"] = split_line[5]  # either '+' or '-'
        if new_entry["strand"] == '+':
            new_entry["start"] = int(split_line[2])
            new_entry["end"] = int(split_line[4])
        else:
            new_entry["start"] = int(split_line[4])
            new_entry["end"] = int(split_line[2])

        new_entry["conf"] = int(split_line[7])
        new_entry["hairpin_score"] = float(split_line[8])
        new_entry["tail_score"] = float(split_line[9])

    return entry_list

def promocalc(inseq, circular_length, threads=1):

    """ Runs the promoter-calculator on the input, returns results and writes them to the outfile """
    # Extract the sequence from the input

    # Run the promoter-calculator
    # Promoter calculator wants the input to be a string
    inseq = str(inseq.seq)
    results = promoter_calculator(inseq, threads=threads, verbosity=0)

    # Write the results to the outfile
    outdata = namedtuple('promoter_calculator_result', 'seq score strand TSSpos box35pos box35seq box10pos box10seq')
    final_list = []
    for result in results:
        # Convert the results to a namedtuple
        result = outdata(result.promoter_sequence,
                         result.Tx_rate, 
                         result.strand,
                         result.TSS-circular_length,
                         result.hex35_position,
                         result.hex35,
                         result.hex10_position,
                         result.hex10)
        final_list.append(result)
    
    # Sort the list by coordinate
    final_list = sorted(final_list, key=attrgetter('TSSpos'))

    # Send the results back to cryptkeeper
    return final_list
