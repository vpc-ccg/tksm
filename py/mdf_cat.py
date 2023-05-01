import argparse
import queue
import threading

from tqdm import tqdm


### Multithreaded version of cat specifically for MDF files
### This scripts is part of enabling Snakemake to use pipes
### for files used by more than one Merge module.
def read_file(input_file, q):
    with open(input_file) as f:
        mdf_lines = list()
        for line in f:
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
    with open(args.output_file, "w") as f, tqdm(
        desc=f"[mdf_cat.py] Cat-ing from {len(args.input_files)} file(s) to {args.output_file}"
    ) as tbar:
        while True:
            for idx, q in enumerate(input_queues):
                try:
                    mdf_lines = q.get(block=False)
                except queue.Empty:
                    pass
                else:
                    if mdf_lines == None:
                        done[idx] = True
                    else:
                        f.write(mdf_lines)
                        tbar.update(1)
            if all(done):
                break
    for t in input_threads:
        t.join()


if __name__ == "__main__":
    main()
