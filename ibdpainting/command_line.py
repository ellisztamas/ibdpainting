#!/usr/bin/env python

"""Console script for methlab."""

import argparse
import ibdpainting as ip
import argparse 
  
def main():
    parser = argparse.ArgumentParser(description='ibdpainting')

    parser.add_argument('-i', '--input',
        help='Path to an HDF5 file containing genotype data for one or more samples to check.  This should be the output of allel.vcf_to_hdf5().'
        )
    parser.add_argument('-n', '--sample_name',
        help ='Sample name for the individual to check. This must be present in the samples in the input file.'
    )
    parser.add_argument('-r', '--reference',
        help="Path to an HDF5 file containing genotype data for a panel of reference individuals to compare the input indivual against. This should be the output of allel.vcf_to_hdf5()."
    )
    parser.add_argument('-w', '--window_size',
        type=int, default=20000,
        help="Integer window size in base pairs."
    )
    parser.add_argument('--expected_match',
        help="Optional list of sample names in the reference panel that are expected to be ancestors of the test individual.",
        nargs = "+", required=False
    )
    parser.add_argument('--outdir',
        help="Directory to save the output."
    )
    parser.add_argument('--keep_ibd_table', 
        help="If set, write an intermediate text file giving genetic distance between the crossed individual and each candidate at each window in the genome. Defaults to False, because these can be quite large.",
        default=False,
        action=argparse.BooleanOptionalAction
    )
    parser.add_argument('--max_to_plot', 
        help="Optional number of the best matching candidates to plot so that the HTML files do not get too large and complicated. Ignored if this is more than the number of samples. Defaults to 20.",
        type=int, default = 10
    )
    parser.add_argument('--interactive',
        help="If set, save the output plot as an interactive HTML plot including information on candidates within the plot.",
        default=True,
        action=argparse.BooleanOptionalAction
        )
    parser.add_argument('--height',
        help="Height in centimetres of the output PNG file.",
        default=675
        )
    parser.add_argument('--width',
        help="Height in centimetres of the output PNG file.",
        default=900)
    args = parser.parse_args()

    # Data frame of IBD at all positions across the genome, and the plot of this
    itable = ip.ibd_table(args.input, args.reference, args.sample_name, args.window_size)
    scores = ip.ibd_scores(itable)
    fig = ip.plot_ibd_table(itable, args.sample_name, args.expected_match, args.max_to_plot)
    
    if args.keep_ibd_table:
        itable.to_csv(args.outdir + "/" + args.sample_name + "_ibd_table.csv", index=False)
    
    scores.to_csv( args.outdir + "/" + args.sample_name + "_ibd_scores.csv", index=False)

    fig.write_image(
        args.outdir + "/" + args.sample_name + "_plot_ibd.png",
        height = args.height, width = args.width
        )
    
    if args.interactive:
        fig.write_html(args.outdir + "/" + args.sample_name + "_plot_ibd.html")
    
        
    

if __name__ == '__main__':
    main()