import argparse
import queue
import threading


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
    with open(args.output_file, "w") as outfile:
        while True:
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
