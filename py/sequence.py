#!/usr/bin/env python3

from glob import glob
import os
import sys
import uuid
import random
import argparse
import functools
from multiprocessing import Pool
import gzip

from tqdm import tqdm
import tksm_badread


def set_tksm_models_dicts(env_var="TKSM_MODELS"):
    var = os.getenv(env_var)
    if var is None:
        return
    for model_dir in reversed(var.split(":")):
        print(f"Loading models from {model_dir}")
        for tool, kind, dictionary in [
            ("badread", "error", tksm_badread.ERROR_MODEL_PY.error_model_names),
            ("badread", "qscore", tksm_badread.QSCOREMODEL_PY.qscore_model_names),
            ("badread", "tail", tksm_badread.TAIL_NOISE_MODEL_PY.tail_model_names),
        ]:
            for path in glob(f"{model_dir}/{tool}/*.{kind}.gz"):
                name = os.path.basename(path)
                name = name[: -len(f".{kind}.gz")]
                dictionary[name] = path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sequencer module of tksm. Generates FASTQ/A file from MDF (molecule description file) and reference FASTA files. "
        + "Note: tksm looks for sequencing model names in colon sperated paths in $TKSM_MODELS environment variable. "
        + f"Current $TKSM_MODELS value: {os.getenv('TKSM_MODELS')}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="MDF file.")
    parser.add_argument(
        "-r",
        "--references",
        type=str,
        nargs="+",
        required=False,
        help="Reference FASTA files, space separated.",
    )
    parser.add_argument(
        "-o",
        "--badread",
        type=str,
        help="Badread reads output file. Default: no output.",
    )
    parser.add_argument(
        "--perfect", type=str, help="Perfect reads output file. Default: no output."
    )
    parser.add_argument(
        "--skip-qual-compute",
        default=False,
        action="store_true",
        help="Use all 'K' for q-scores. Default: compute quality score using Edlib alignment to input sequences.",
    )
    parser.add_argument(
        "-O",
        "--output-format",
        type=str,
        default=None,
        choices=["fastq", "fasta"],
        help="Output format. Default: infer from output file extension and fallback to FASTA. Accepts .gz files too.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads.",
    )
    parser.add_argument(
        "--badread-identity",
        type=str,
        default="84.0,99.0,5.5",
        help="Sequencing identity distribution (mean, max and stdev)",
    )
    parser.add_argument(
        "--badread-error-model",
        type=str,
        default=(
            "nanopore2020"
            if "nanopore2020" in tksm_badread.ERROR_MODEL_PY.error_model_names
            else "random"
        ),
        help="Badread tail model name or file path. "
        + f"Available model names: [{', '.join(tksm_badread.ERROR_MODEL_PY.error_model_names)}]",
    )
    parser.add_argument(
        "--badread-qscore-model",
        type=str,
        default=(
            "nanopore2020"
            if "nanopore2020" in tksm_badread.QSCOREMODEL_PY.qscore_model_names
            else "random"
        ),
        help="Badread qscore model name or file path. "
        + f"Available model names: [{', '.join(tksm_badread.QSCOREMODEL_PY.qscore_model_names)}]",
    )
    parser.add_argument(
        "--badread-tail-model",
        type=str,
        default="no_noise",
        help="Badread tail model name or file path. "
        + f"Available model names: [{', '.join(tksm_badread.TAIL_NOISE_MODEL_PY.tail_model_names)}]",
    )

    class ListPrinter(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            txt = "\n".join([getattr(k, "dest") for k in parser._actions])
            print(txt)
            parser.exit()

    parser.add_argument("--list", nargs=0, action=ListPrinter)

    args = parser.parse_args()

    if args.list:

        options = parser._actions
        for k in options:
            print(getattr(k, "dest"))
        exit(0)

    # Process arguments and check for errors
    try:
        identity_parameters = [float(x) for x in args.badread_identity.split(",")]
        assert (
            len(identity_parameters) == 3
        ), "Must specify 3 values for --badread-identity"
        args.badread_mean_identity = identity_parameters[0]
        args.badread_max_identity = identity_parameters[1]
        args.badread_identity_stdev = identity_parameters[2]
    except (ValueError, IndexError):
        sys.exit("Error: could not parse --identity values")
    if args.badread_mean_identity > 100.0:
        sys.exit("Error: mean read identity cannot be more than 100")
    if args.badread_max_identity > 100.0:
        sys.exit("Error: max read identity cannot be more than 100")
    if args.badread_mean_identity <= tksm_badread.SETTINGS_PY.MIN_MEAN_READ_IDENTITY:
        sys.exit(
            f"Error: mean read identity must be at least {tksm_badread.SETTINGS_PY.MIN_MEAN_READ_IDENTITY}"
        )
    if args.badread_max_identity <= tksm_badread.SETTINGS_PY.MIN_MEAN_READ_IDENTITY:
        sys.exit(
            f"Error: max read identity must be at least {tksm_badread.SETTINGS_PY.MIN_MEAN_READ_IDENTITY}"
        )
    if args.badread_mean_identity > args.badread_max_identity:
        sys.exit(
            f"Error: mean identity ({args.badread_mean_identity}) cannot be larger than max "
            f"identity ({args.badread_max_identity})"
        )
    if args.badread_identity_stdev < 0.0:
        sys.exit("Error: read identity stdev cannot be negative")
    if not args.badread and not args.perfect:
        parser.error("Must specify either --output or --perfect.")
    return args


def generate_fasta(fasta):
    name = ""
    seq = list()
    if fasta.endswith(".gz"):
        f = gzip.open(fasta, "rt")
    else:
        f = open(fasta, "r")
    for _, l in enumerate(f):
        l = l.rstrip("\n")
        if l[0] == ">":
            if len(seq) == 0:
                name = l[1:].split(" ")[0]
                continue
            yield (name, "".join(seq))
            seq = list()
            name = l[1:].split(" ")[0]
        else:
            seq.append(l)
    yield (name, "".join(seq))


def get_reference_seqs(reference):
    reference_seqs = dict()
    for ref in reference:
        print(f"Loading reference {ref}...")
        reference_seqs.update({name: seq for name, seq in generate_fasta(ref)})
    return reference_seqs


def mdf_generator(f):
    # MDF format:
    # Each entry starts with a header: +<mol_id>\t<depth>\t<comments>;
    #   where <mol_id> is the molecule ID, <depth> is the depth of the molecule, and each <comments> entry is key=value pair terminated by a semicolon.
    # Each line after the header, and until the next header or end of file, is an interval description: <ref_id>\t<start:int>\t<end:int>\t<strand>\t<modifications>
    read_id = None
    intervals = list()
    depth: int = 0
    for line in f:
        line = line.strip("\n").split("\t")
        if line[0][0] == "+":
            if read_id is not None:
                for _ in range(depth):
                    yield read_id, intervals
            read_id = line[0][1:]
            depth = int(line[1])
            intervals = list()
        else:
            chrom, start, end, strand, modifications = line

            start, end = int(start), int(end)
            intervals.append((chrom, start, end, strand, modifications))
    if read_id is not None:
        for _ in range(depth):
            yield read_id, intervals


def reverse_complement(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def apply_modifications(seq, modifications):
    if modifications == "":
        return seq
    else:
        modifications = modifications.split(",")
        seq = list(seq)
        for mod in modifications:
            char = mod[-1]
            pos = int(mod[:-1])
            seq[pos] = char
        return "".join(seq)


def badread(molecule_id, raw_seq, compute_qual=False):
    target_identity = identities.get_identity()
    (seq, quals, actual_identity, _) = tksm_badread.SIMULATE_PY.sequence_fragment(
        raw_seq,
        target_identity,
        error_model,
        qscore_model,
        tail_model,
        compute_qual,
    )
    info = [
        f"length={len(seq)}",
        f"error_free_length={len(raw_seq)}",
        f"read_identity={actual_identity * 100.0:.2f}%",
        f"molecule_id={molecule_id}",
    ]
    return (seq, quals, info)


def perfect(molecule_id, seq):
    quals = "K" * len(seq)
    info = f"length={len(seq)} molecule_id={molecule_id}"
    info = [
        f"length={len(seq)}",
        f"error_free_length={len(seq)}",
        f"read_identity={1 * 100.0:.2f}%",
        f"molecule_id={molecule_id}",
    ]
    return (seq, quals, info)


def fastq_formatter(read_id, seq, quals, info):
    result = [
        f"@{read_id} {' '.join(info)}",
        seq,
        "+",
        quals,
    ]
    return "\n".join(result)


def fasta_formatter(read_id, seq, _, info):
    result = [
        f">{read_id} {' '.join(info)}",
        seq,
    ]
    return "\n".join(result)


def get_output_file(outfile):
    if outfile.endswith(".gz"):
        f = gzip.open(outfile, "wt+")
        outfile = outfile[:-3]
    else:
        f = open(outfile, "w+")
    if outfile.endswith(".fastq") or outfile.endswith(".fq"):
        return f, fastq_formatter
    else:
        return f, fasta_formatter


def mdf_to_seq(mdf, targets=dict()):
    molecule_id, intervals = mdf
    seq = list()
    for chrom, start, end, strand, modifications in intervals:
        segment = reference_seqs.get(chrom, chrom)[start:end].upper()
        segment = apply_modifications(segment, modifications)
        if strand == "+":
            seq.append(segment)
        else:
            seq.append(reverse_complement(segment))
    seq = "".join(seq)

    results = dict()
    read_id = uuid.UUID(int=random.getrandbits(128))
    for k, (seq_producer, read_formatter) in targets.items():
        seq, quals, info = seq_producer(molecule_id, seq)
        results[k] = read_formatter(read_id, seq, quals, info)
    return results


if __name__ == "__main__":
    set_tksm_models_dicts()
    args = parse_args()
    reference_seqs = get_reference_seqs(args.references)

    targets = dict()
    target_outfiles = dict()
    if args.badread:
        identities = tksm_badread.IDENTITIES_PY.Identities(
            args.badread_mean_identity,
            args.badread_identity_stdev,
            args.badread_max_identity,
            sys.stderr,
        )
        error_model = tksm_badread.ERROR_MODEL_PY.ErrorModel(
            args.badread_error_model, sys.stderr
        )
        qscore_model = tksm_badread.QSCOREMODEL_PY.QScoreModel(
            args.badread_qscore_model, sys.stderr
        )
        tail_model = tksm_badread.TAIL_NOISE_MODEL_PY.KDE_noise_generator.load(
            args.badread_tail_model
        )

        output_file, output_file_formatter = get_output_file(args.badread)
        target_outfiles["badread"] = output_file
        partial_badread = functools.partial(
            badread,
            compute_qual=(not args.skip_qual_compute)
            and output_file_formatter == fastq_formatter,
        )
        targets["badread"] = [partial_badread, output_file_formatter]
    if args.perfect:
        output_file, output_file_formatter = get_output_file(args.perfect)
        target_outfiles["perfect"] = output_file
        targets["perfect"] = [perfect, output_file_formatter]

    mdf_to_seq_targeted = functools.partial(mdf_to_seq, targets=targets)
    with open(args.input, "r") as f:
        mdg = mdf_generator(f)
        if args.threads > 1:
            p = Pool(args.threads)
            mapper = functools.partial(p.imap_unordered, chunksize=10)
        else:
            mapper = map
        for read_dict in tqdm(mapper(mdf_to_seq_targeted, mdg), desc="Sequencing"):
            for k, v in read_dict.items():
                print(v, file=target_outfiles[k])

        if args.threads > 1:
            p.close()  # type: ignore

    for v in target_outfiles.values():
        v.close()
