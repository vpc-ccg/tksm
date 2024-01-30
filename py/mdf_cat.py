import argparse
from datetime import datetime
import queue
import sys
import threading

from tqdm import tqdm


### Multithreaded version of cat specifically for MDF files
### This scripts is part of enabling Snakemake to use pipes
### for files used by more than one Merge module.
def read_file(input_file, q):
    with open(input_file) as infile:
        mdf_lines = list()
        for line in infile:
            if line.startswith("+"):
                q.put("".join(mdf_lines))
                mdf_lines = list()
            mdf_lines.append(line)
        q.put("".join(mdf_lines))
    q.put(None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_files", nargs="+", help="input file names")
    parser.add_argument("output_file", help="output file name")
    args = parser.parse_args()

    input_queues = [queue.Queue() for _ in args.input_files]

    input_threads = [
        threading.Thread(target=read_file, args=(input_file, q))
        for input_file, q in zip(args.input_files, input_queues)
    ]
    for t in input_threads:
        t.start()

    done = [False for _ in input_queues]
    md_count = 0
    last_time = datetime.now()
    with open(args.output_file, "w") as outfile:
        while True:
            # if (datetime.now() - last_time).total_seconds() > 2:
            #     last_time = datetime.now()
            #     sys.stderr.write(
            #         f"[mdf_cat.py] {[int(x) for x in done]} files; {tuple(x.split('/')[-1] for x in args.input_files)}; {md_count} MDs to {args.output_file.split('/')[-2:]}.\n"
            #     )
            for idx, (d, q) in enumerate(zip(done, input_queues)):
                if d:
                    continue
                try:
                    mdf_lines = q.get(block=False)
                except queue.Empty:
                    pass
                else:
                    if mdf_lines == None:
                        done[idx] = True
                    else:
                        outfile.write(mdf_lines)
                        md_count += 1
            if all(done):
                break
    for t in input_threads:
        t.join()


if __name__ == "__main__":
    main()

