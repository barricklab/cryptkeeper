"""
Simple wrapper for the promoter-calculator that takes the output, makes a file, and also formats it
in a way that cryptkeeper likes
"""
import csv
from collections import namedtuple
from operator import attrgetter
from promoter_calculator import promoter_calculator
from .helpers import persistant_cache

#@persistant_cache(filepath='tts_predict_cache.db')
def call_predict(inseq, quiet):
    """ Calls the promoter-calculator on the input, returns results """
    # Extract the sequence from the input
    inseq = str(inseq)
    results = promoter_calculator(inseq, quiet)
    return results

def tts_predict(inseq, outfile):
    """ Runs the promoter-calculator on the input, returns results and writes them to the outfile """
    # Extract the sequence from the input

    # Run the promoter-calculator
    # Promoter calculator wants the input to be a string
    inseq = str(inseq.seq)
    results = call_predict(inseq, quiet=True)
    print(len(results))

    # Write the results to the outfile
    outdata = namedtuple('tss_prediction', 'seq score strand TSSpos box35pos box35seq box10pos box10seq')
    final_list = []
    for result in results:
        # Convert the results to a namedtuple
        result = outdata(result.promoter_sequence,
                         result.Tx_rate, 
                         result.strand,
                         result.TSS,
                         result.hex35_position,
                         result.hex35,
                         result.hex10_position,
                         result.hex10)
        final_list.append(result)
    
    # Sort the list by coordinate
    final_list = sorted(final_list, key=attrgetter('TSSpos'))
    print(len(final_list))


    with open(outfile,'w') as final_predictions_file:
      writer = csv.DictWriter(
          final_predictions_file,
          fieldnames = ["promoter", "score", "strand", "TSSpos", "box35pos", "box35seq", "box10pos", "box10seq",]
        )
      writer.writeheader()
      writer.writerows([{"promoter": data.seq,
                        "score": data.score,
                        "strand": data.strand,
                        "TSSpos": data.TSSpos,
                        "box35pos": data.box35pos,
                        "box35seq": data.box35seq,
                        "box10pos": data.box10pos,
                        "box10seq": data.box10seq} for data in final_list])

    # Send the results back to cryptkeeper
    return final_list