#! /usr/bin/env python3
'''
PyTreeDate: A Python tool for dating a rooted phylogenetic tree using various molecular clock assumptions
'''
from datetime import datetime,timedelta
from gzip import open as gopen
from scipy.stats import linregress
from sys import stdin
from treeswift import read_tree_newick
import argparse

# convert a date (YYYY-MM-DD) to days since 0001-01-01
def date_to_days(sample_time):
    num_dashes = sample_time.count('-')
    if num_dashes == 2:   # YYYY-MM-DD
        tmp = datetime.strptime(sample_time, '%Y-%m-%d')
    elif num_dashes == 1: # YYYY-MM(-01)
        tmp = datetime.strptime('%s-01' % sample_time, '%Y-%m-%d')
    elif num_dashes == 0: # YYYY(-01-01)
        tmp = datetime.strptime('%s-01-01' % sample_time, '%Y-%m-%d')
    else:
        raise ValueError("Invalid sample date (should be YYYY-MM-DD): %s" % sample_time)
    return (tmp - datetime(1,1,1)).days # days since 0001-01-01

# convert days since 0001-01-01 to a date (YYYY-MM-DD)
def days_to_date(days):
    return (datetime(1,1,1) + timedelta(days=days)).strftime('%Y-%m-%d')

# parse a given sample dates file
def parse_dates(dates_filename):
    if dates_filename.lower().endswith('.gz'):
        lines = [l.strip() for l in gopen(dates_filename).read().decode().strip().splitlines()]
    else:
        lines = [l.strip() for l in open(dates_filename).read().strip().splitlines()]
    out = dict()
    for l in lines:
        if len(l) == 0:
            continue
        u,t = l.split('\t')
        out[u] = date_to_days(t)
    return out

# date the given tree using a strict molecular clock
def date_strict(tree, dates):
    rtt = dict(); x = list(); y = list() # x is time; y is root-to-tip
    for node in tree.traverse_preorder():
        if node.is_root():
            rtt[node] = 0
        else:
            rtt[node] = rtt[node.parent]
            if node.edge_length is not None:
                rtt[node] += node.edge_length
        if node.is_leaf():
            x.append(dates[node.label]); y.append(rtt[node])
    slope, y_intercept, r_value, p_value, std_err = linregress(x,y) # slope is mutations/time, x-intercept is tMRCA
    for node in tree.traverse_preorder():
        if node.edge_length is not None:
            node.edge_length /= slope
    return tree

MODES = {
    'strict':date_strict,
}

if __name__ == "__main__":
    # parse user args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_tree', required=False, type=str, default='stdin', help="Input Rooted Mutation Tree (Newick format)")
    parser.add_argument('-d', '--input_dates', required=True, type=str, help="Input Sample Dates (TSV, left column IDs, right column YYYY-MM-DD)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Rooted Time Tree (Newick format)")
    parser.add_argument('-m', '--mode', required=False, type=str, default='strict', help="Dating Mode (%s)" % ', '.join(MODES.keys()))
    args = parser.parse_args()
    if args.mode.lower() not in MODES:
        raise ValueError("Invalid mode: %s (valid options: %s" % (args.mode, ', '.join(MODES.keys())))
    if args.input_tree == 'stdin':
        tree_mut = read_tree_newick(stdin.read())
    else:
        tree_mut = read_tree_newick(args.input_tree)
    dates = parse_dates(args.input_dates)

    # date tree and output
    tree_time = MODES[args.mode.lower()](tree_mut, dates)
    if args.output == 'stdout':
        print(tree_time.newick())
    elif args.output.lower().endswith('.gz'):
        f = gopen(args.output,'wb'); f.write(tree_time.newick().encode()); f.close()
    else:
        f = open(args.output,'w'); f.write(tree_time.newick()); f.close()
