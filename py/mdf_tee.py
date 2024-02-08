import argparse
import queue
import threading

from tqdm import tqdm


### Multithreaded version of tee specifically for MDF files
### This scripts is part of enabling Snakemake to use pipes
### for files used by more than one Merge module.
def write_file(output_file, q):
    with open(output_file, "w") as f:
        while True:
            line = q.get()
            if line is None:
                break
            f.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="input file name")
    parser.add_argument("output_files", nargs="+", help="output file names")
    args = parser.parse_args()

    queues = [queue.Queue() for _ in args.output_files]
    threads = [
        threading.Thread(target=write_file, args=(output_file, q))
        for output_file, q in zip(args.output_files, queues)
    ]
    for t in threads:
        t.start()

    with open(args.input_file) as f, tqdm(
        f,
        desc=f"[mdf_tee.py] Tee-ing to {len(args.output_files)} file(s) from {args.input_file}",
    ) as tbar:
        mdf_lines = list()
        for line in f:
            if line.startswith("+"):
                tbar.update(1)
                for q in queues:
                    q.put("".join(mdf_lines))
                mdf_lines = list()
            mdf_lines.append(line)
        for q in queues:
            q.put("".join(mdf_lines))

    for q in queues:
        q.put(None)
    for t in threads:
        t.join()


if __name__ == "__main__":
    main()
