#!/usr/bin/env python3
"""
Convert p5 barcode to p3 barcode
"""
import argparse
import json
import os
import logging
import gzip
import collections

import pyfastx
import utils
import sys

# stdout
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO, stream=sys.stdout)
logger = logging.getLogger(__name__)



def get_protocol_dict(assets_dir):
    """
    Return:
    protocol_dict. Key: protocol name, value: values from protocol.json

    >>> protocol_dict = get_protocol_dict("./assets/")
    >>> protocol_dict["AccuraSCOPE-V1"]["pattern_dict_p5"]
    {'C': [slice(0, 6, None)], 'U': [slice(6, 11, None)]}
    """
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        for x in ["bc_p3", "bc_p5"]:
            cur[x] = os.path.join(whitelist_dir, protocol, cur[x])
        cur["pattern_dict_p3"] = utils.parse_pattern(cur["pattern_p3"])
        cur["pattern_dict_p5"] = utils.parse_pattern(cur["pattern_p5"])
    return protocol_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--fq1', required=True)
    parser.add_argument('--fq2', required=True)
    parser.add_argument('--assets_dir', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--threshold', default=1e-3, type=float)
    args = parser.parse_args()

    protocol_dict = get_protocol_dict(args.assets_dir)
    protocol = protocol_dict[args.protocol]
    bc_p3 = utils.read_one_col(protocol["bc_p3"])
    bc_p5 = utils.read_one_col(protocol["bc_p5"])
    bc_p3_mismatch_dict = utils.get_mismatch_dict(bc_p3)
    bc_p5_mismatch_dict = utils.get_mismatch_dict(bc_p5)
    # p3 p5 bc correspondance
    bc_p5_p3_dict = {p5:p3 for p5,p3 in zip(bc_p5, bc_p3)}

    pattern_dict_p3 = protocol["pattern_dict_p3"]
    pattern_dict_p5 = protocol["pattern_dict_p5"]
    outdict = {}
    well_read_count = collections.defaultdict(lambda: collections.defaultdict(int))
    total_reads = 0

    fq1_list = args.fq1.split(',')
    fq2_list = args.fq2.split(',')
    """
    write seq2 to `{sample}_{bc}.fq.gz`. reverse complement seq2 if it is p5
    """
    for fq1,fq2 in zip(fq1_list, fq2_list):
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)
        for (name1, seq1, qual1), (name2,seq2,qual2) in zip(fq1, fq2):
            total_reads += 1
            prime = 'invalid'
            bc_p3 = utils.get_seq_str(seq1, pattern_dict_p3["C"])
            if bc_p3 in bc_p3_mismatch_dict:
                prime = 'p3'
                bc = bc_p3_mismatch_dict[bc_p3]
                umi = utils.get_seq_str(seq1, pattern_dict_p3["U"])
            else:
                bc_p5 = utils.get_seq_str(seq1, pattern_dict_p5["C"])
                if bc_p5 in bc_p5_mismatch_dict:
                    prime = 'p5'
                    bc = bc_p5_p3_dict[bc_p5_mismatch_dict[bc_p5]]
                    umi = utils.get_seq_str(seq1, pattern_dict_p5["U"]) 

            if prime in ['p3','p5']:
                well_read_count[bc][prime] += 1
                read_name = ":".join([prime, str(total_reads), umi])
                if bc not in outdict:
                    outdict[bc] = utils.openfile(f'{args.sample}_{bc}.fq.gz', 'wt')
                if prime == 'p5':
                    seq2 = utils.reverse_complement(seq2)
                    qual2 = qual2[::-1]
                outdict[bc].write(utils.str_fq(read_name, seq2, qual2))
    
    """
    move all bc fastq files to `pass` folder it has reads > total_reads * threshold
    move others to `below`
    """
    if not os.path.exists('pass'):
        os.mkdir('pass')
        os.mkdir('below')
    for bc in outdict:
        fn = f'{args.sample}_{bc}.fq.gz'
        if (well_read_count[bc]['p3'] + well_read_count[bc]['p5']) >= total_reads * args.threshold:
            os.rename(fn, f'pass/{fn}')
        else:
            os.rename(fn, f'below/{fn}')
    # dump well_read_count json
    with open(f'{args.sample}_well_read_count.json', 'w') as f:
        json.dump(well_read_count, f, indent=4)


if __name__ == "__main__":
    main()








