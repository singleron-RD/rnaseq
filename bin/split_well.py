#!/usr/bin/env python3
"""
Convert p5 barcode to p3 barcode
"""

import argparse
import json
import os
import logging
import collections

import pysam
import utils
import sys
import pandas as pd

# stdout
logging.basicConfig(
    format="%(asctime)s - %(message)s", level=logging.INFO, stream=sys.stdout
)
logger = logging.getLogger(__name__)


def get_protocol_dict(assets_dir):
    """
    Return:
    protocol_dict. Key: protocol name, value: values from protocol.json

    >>> protocol_dict = get_protocol_dict("./assets/")
    >>> protocol_dict["AccuraSCOPE-V1"]["pattern_dict_p5"]
    {'C': [slice(0, 6, None)]}
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


def get_well_barcode(
    bc_file: str,
) -> dict[int, str]:
    barcodes = utils.one_col_to_list(bc_file)
    return {i: x for i, x in enumerate(barcodes, start=1)}


def get_barcode_sample(
    bc_file: str,
    well_sample_file: str,
) -> dict[str, str]:
    well_barcode = get_well_barcode(bc_file)
    well_sample = utils.two_col_to_dict(well_sample_file)
    barcode_sample = {}
    for well, barcode in well_barcode.items():
        if well in well_sample:
            barcode_sample[barcode] = well_sample[well]
        else:
            barcode_sample[barcode] = f"noise_well_{well}"
    return barcode_sample


def get_barcode_in_sample(bc_file, well_sample_file):
    well_barcode = get_well_barcode(bc_file)
    well_sample = utils.two_col_to_dict(well_sample_file)
    barcode_in_sample = []
    for well, barcode in well_barcode.items():
        if well in well_sample:
            barcode_in_sample.append(barcode)
    return barcode_in_sample


def get_mismatch_dict(bc_file, well_sample_file):
    """
    Correct the barcodes to sample barcodes rather than noise barcodes as much as possible.
    """
    barcode_in_sample = get_barcode_in_sample(bc_file, well_sample_file)
    all_bc = utils.one_col_to_list(bc_file)
    mismatch_dict = utils.get_mismatch_dict(all_bc)
    sample_mismatch_dict = utils.get_mismatch_dict(barcode_in_sample)
    mismatch_dict.update(sample_mismatch_dict)
    return mismatch_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument(
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
        required=True,
    )
    args = parser.parse_args()

    protocol_dict = get_protocol_dict(args.assets_dir)
    protocol = protocol_dict[args.protocol]
    # mismatch_dict
    p3_mismatch_dict = get_mismatch_dict(protocol["bc_p3"], args.well_sample)
    p5_mismatch_dict = get_mismatch_dict(protocol["bc_p5"], args.well_sample)
    p3_linker = "ATACGCGGA"
    p3_linker_mismatch_dict = utils.get_mismatch_dict([p3_linker])
    # barcode sample
    p3_barcode_sample = get_barcode_sample(protocol["bc_p3"], args.well_sample)
    p5_barcode_sample = get_barcode_sample(protocol["bc_p5"], args.well_sample)
    # pattern
    p3_pattern_dict = protocol["pattern_dict_p3"]
    p5_pattern_dict = protocol["pattern_dict_p5"]
    # metrics and output
    if not os.path.exists("signal"):
        os.mkdir("signal")
        os.mkdir("noise")
    outdict = {}
    sample_read_count = collections.defaultdict(lambda: collections.defaultdict(int))
    total_reads = 0
    p3_reads = 0
    p5_reads = 0
    signal_reads = 0

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    for fq1, fq2 in zip(fq1_list, fq2_list):
        fq1 = pysam.FastxFile(fq1)
        fq2 = pysam.FastxFile(fq2)
        for entry1, entry2 in zip(fq1, fq2):
            seq1 = entry1.sequence
            seq2, qual2 = entry2.sequence, entry2.quality
            total_reads += 1
            prime = "invalid"
            p3_linker = utils.get_seq_str(seq1, p3_pattern_dict["L"])
            if p3_linker in p3_linker_mismatch_dict:
                bc_p3 = utils.get_seq_str(seq1, p3_pattern_dict["C"])
                if bc_p3 in p3_mismatch_dict:
                    prime = "p3"
                    bc = p3_mismatch_dict[bc_p3]
                    umi = utils.get_seq_str(seq1, p3_pattern_dict["U"])
                else:
                    continue
            else:
                bc_p5 = utils.get_seq_str(seq1, p5_pattern_dict["C"])
                if bc_p5 in p5_mismatch_dict:
                    prime = "p5"
                    bc = p5_mismatch_dict[bc_p5]
                    umi = "-"
                else:
                    continue

            if prime == "p3":
                sample = p3_barcode_sample[bc]
                p3_reads += 1
            else:
                sample = p5_barcode_sample[bc]
                p5_reads += 1
            sample_read_count[sample][prime] += 1
            read_name = ":".join([prime, str(total_reads), umi])
            folder = "noise" if sample.startswith("noise") else "signal"
            if folder == "signal":
                signal_reads += 1
            if sample not in outdict:
                outdict[sample] = utils.openfile(
                    f"{folder}/{args.sample}_{sample}.fq.gz", "wt", compresslevel=1
                )
            outdict[sample].write(utils.str_fq(read_name, seq2, qual2))

    df = pd.DataFrame.from_dict(sample_read_count, orient="index")
    df = df.reindex(columns=["p3", "p5"]).fillna(0).astype(int)
    df.sort_values(by=["p3", "p5"], ascending=False, inplace=True)
    df.to_csv(f"{args.sample}_read_count.tsv", index_label="sample", sep="\t")

    p3_percent = round(p3_reads / total_reads * 100, 2)
    p5_percent = round(p5_reads / total_reads * 100, 2)
    signal_percent = round(signal_reads / total_reads * 100, 2)
    metrics_dict = {
        "total_reads": total_reads,
        "p3_reads": f"{p3_reads} ({p3_percent}%)",
        "p5_reads": f"{p5_reads} ({p5_percent}%)",
        "signal_reads": f"{signal_reads} ({signal_percent}%)",
    }
    with open(f"{args.sample}_metrics.txt", "wt") as f:
        for key, value in metrics_dict.items():
            f.write(f"{key}:{value}\n")


if __name__ == "__main__":
    main()
