#!/usr/bin/env python3
import os
import argparse
import logging
import logging.config
import multiprocessing

from Bio import SeqIO

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence

ALLOWED_DETECTORS = {
    "HHrepID", 
    "T-REKS",
    "TRUST", 
    "XSTREAM",
    "TRED",
    "TRF"
}

def detect_trs(record, output_dir, seq_type, detectors):
    print("Started work on sequence {}".format(record.id))

    seq = sequence.Sequence(seq=str(record.seq), name=record.id)
    denovo_list = seq.detect(
        denovo=True, 
        sequence_type=seq_type,
        detection = {"detectors": detectors}     
    )
    denovo_list.repeats = sorted(denovo_list.repeats, key=lambda x: x.begin)
    seq.set_repeatlist(denovo_list, tag="denovo")

    print("Detected {} repeat(s) in sequence {} Writing results to {}".format(len(seq.get_repeatlist("denovo").repeats), seq.name, output_dir))

    output_file_name = os.path.join(output_dir, record.id + ".pickle")
    seq.get_repeatlist(tag="denovo").write(output_format="pickle", file=output_file_name)

def check_output_dir(output_dir):
    """
    Checks the output directory for files generated in previous runs, these can be skipped later by detect_trs()
    Checking is done quite naively, only looking for files ending in '.pickle' (so no support for .pcl, .pkl ...)
    Parameters:
    output_dir (str):   Directory to check for output from previous runs
    Returns:
    finished_sequences (set):   Set of genomic regions that can be skipped by detect_trs()
    """
    finished_sequences = {i.replace(".pickle", "") for i in os.listdir(output_dir) if i.endswith(".pickle")}
    return finished_sequences

def load_sequences(fasta_file, output_dir):
    finished_sequences = check_output_dir(output_dir)
    sequences = [seq for seq in SeqIO.parse(fasta_file, "fasta") if not seq.id in finished_sequences]
    return sequences

def get_n_cpus(cla_n_cpus, n_seqs):
    n_cpus = cla_n_cpus
    if n_cpus == -1:
        n_cpus = min(n_seqs, os.cpu_count())
    elif n_cpus >= 1:
        n_cpus = min(n_seqs, os.cpu_count(), n_cpus)
    else:
        raise ValueError("--processes must be positive or -1")    
    return n_cpus

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta_file", "-f", type=str, required=True, help="Path to input fasta file"
    )
    parser.add_argument(
        "--sequence_type", "-s", type=str, required=True, help="Type of sequence that is being analyzed (currently supports 'AA' and 'DNA')"
    )    
    parser.add_argument(
        "--output_dir", "-o", type=str, required=True, help="Path to directory where output TRs will be deposited"
    )
    parser.add_argument(
        "--processes", "-p", type=int, default=1, help="Number of repeat detection processes to launch in parallel. If -1, all available CPUs will be used. default: 1."
    )
    parser.add_argument(
        "--detectors", "-d", nargs="+", required=False, help="Whitespace separated list of repeat detectors to use"
    )

    return parser.parse_args()


def main():
    args = parse_cla()
    if args.detectors:
        if not all([i in ALLOWED_DETECTORS for i in args.detectors]):
            raise ValueError(f"Unrecognized detector specified, allowed detectors: \n\t{ALLOWED_DETECTORS}")

    sequences = load_sequences(args.fasta_file, args.output_dir)
    n_cpus = get_n_cpus(args.processes, len(sequences))

    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    print("""
        -------------------------------------------------------
        Launching {} TRAL job(s) with the following parameters:
        Sequence type: {}
        Tandem repeat detectors: 
            {}
        -------------------------------------------------------
        Parameters can be set in ~/.tral/config.ini or
        at the command line
        
        Commencing run
        -------------------------------------------------------
        """.format(n_cpus, args.sequence_type, args.detectors))     

    with multiprocessing.Pool(n_cpus) as pool:
        tasks = [(seq, args.output_dir, args.sequence_type, args.detectors) for seq in sequences]
        pool.starmap(detect_trs, tasks)


if __name__ == "__main__":
    main()
