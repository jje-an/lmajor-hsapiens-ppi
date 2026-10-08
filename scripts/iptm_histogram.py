import matplotlib.pyplot as plt
import numpy as np 
import sqlite3
import sys

def iptm_histogram(db_dir, output_dir="iptm_hist_9.24.26.png"):
    """
    Creates a log scale histogram of ipTM values from a input CSV file. 
    Input file can be created by esm_summaries_to_csv.py.

    Example usage:
    `python iptm_histogram.py esm_summary.csv`

    Arguments:
        db_dir - File path of directory to grab metrics from.
        output_dir - File path of output image file. matplotlib typically
            supports png, pfg, ps, eps, and svg file extensions.
    """
    
    con = sqlite3.connect(db_dir)
    cur = con.cursor()
    
    cur.execute("SELECT iptm FROM esm_summaries WHERE iptm IS NOT NULL")
    iptms = np.array([item[0] for item in cur.fetchall()])

    plt.xlabel("ipTM")
    plt.ylabel("Count")
    plt.title(f"ipTM of {len(iptms)} Structure Predictions")
    plt.hist(iptms, bins=50, log=True)
    plt.savefig(output_dir)
    plt.close()

    con.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(f"Usage: {sys.argv[0]} <db_dir> [output_dir]")
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "iptm_hist.png"
    iptm_histogram(sys.argv[1], output_csv)
