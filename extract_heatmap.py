import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from colorama import init, Fore, Style
from fuzzywuzzy import fuzz, process
import subprocess
import sys
import os
import math

def neighborhood(G, node, n):
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    return [node for node, length in path_lengths.items()
                    if length <= n]
def filter_graph(G, node, d, lr_threshold, p_threshold):
    edges = []
    for u,v,e in G.edges(data=True):
        if e['lr'] >= lr_threshold and e['p'] <= p_threshold:
            edges.append((u,v))
    H=G.edge_subgraph(edges)
    if node in H.nodes:
        return H.subgraph(neighborhood(H, node, d))
    return G.subgraph([node])

def sort_gene_names(gene_names):
    terminals = set( [gene[-1] for gene in gene_names] )
    sets = {k: [] for k in terminals}
    for gene in gene_names:
        sets[gene[-1]].append(gene)

    sorted_genes = []
    for key in sorted(sets.keys()):
        for gene in sorted(sets[key]):
            sorted_genes.append(gene)
    return sorted_genes

def escape_brackets(string):
    t = string.replace('(', '\(')
    t = t.replace(')', '\)')
    return t

## user input stuff
#allow the user to exit anytime
def check_exit(in_text):
    if in_text == 'exit':
        print("You have chosen to exit.")
        sys.exit()
    return in_text

# prompts
gene_prompt = "Please enter the name of the target gene:"
lr_prompt = "Please enter the lower bound of the likelihood ratio:"
p_prompt = "Please enter the upper bound of the p-value:"
d_prompt = "Please enter the neighborhood depth (as an integer):"

if __name__ == '__main__':
    cwd = os.getcwd()

    try:
        print("Loading main presence absence table. . .")
        pa = pd.read_table('data/all_categories_PA.csv', sep=',')
        print("Loading genome label informtation. . .")
        meta = pd.read_table('data/metadata_with_clades.csv',sep=',')
        meta.rename({'Run_accession': 'Isolate'}, axis=1,inplace=True)
        meta.set_index('Isolate', inplace=True)
        print("Loading graph. This may take a moment. . .")
        G = nx.graphml.read_graphml('./data/pagel_graph_relabelled.graphml')
        print("Data Loaded.")
    except:
        raise FileNotFoundError("Files must be located in a folder called ./data")


    #loop for user input
    print("Creating your data. Please answer the questions. You may type 'exit' any time to exit.")
    node_list = list(G.nodes)
    target = False
    while not target:
        gene = check_exit(input(Fore.GREEN + gene_prompt + Style.RESET_ALL))
        target = gene in G.nodes
        if not target:
            fuzzy_match = process.extract(gene, node_list, limit=3)
            print("{} not found in the network. Please double check spelling and capitalization.".format(gene))
            print("Did you mean:")
            for i in fuzzy_match:
                print(Fore.YELLOW + "\t{}".format(i[0]) + Style.RESET_ALL)


    p_valid = False
    while not p_valid:
        p = check_exit(input(Fore.GREEN + p_prompt + Style.RESET_ALL))
        try:
            p = float(p)
            assert p >= 0.0 and p <= 1.0
            p_valid = True
        except:
            print("p value must be a number between 0.0 and 1.0")

    lr_valid = False
    while not lr_valid:
        lr = check_exit(input(Fore.GREEN + lr_prompt + Style.RESET_ALL))
        try:
            lr = float(lr)
            lr_valid = True
        except:
            print("LR value must be numeric.")

    degree_valid = False
    while not degree_valid:
        d = check_exit(input(Fore.GREEN + d_prompt + Style.RESET_ALL))
        try:
            d = int(d)
            assert type(d) == int
            assert d > 0
            degree_valid = True
        except:
            print("Degree must be a positive integer.")

    print("Input done. Processing. . .")
    #make results directory
    results_dir = 'results/{}'.format(gene)
    if not os.path.exists('results'):
        os.makedirs('results')
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    outfile_csv = '{0}/{1}_neighborhood_PA.csv'.format(results_dir, gene)
    heatmapfile = '{0}/{1}_PA_Heatmap.png'.format(results_dir, gene)
    outfile_graphml = '{0}/{1}_neighborhood_graph.graphml'.format(results_dir, gene)
    outfile_habitats = '{0}/{1}_habitat_dist.png'.format(results_dir, gene)

    print("Computing subgraph. . .")
    S = filter_graph(G, gene, d, lr, p)
    print("Done. ")

    print(Fore.YELLOW + "Creating presence absence table.. This will save a file to your disk." + Style.RESET_ALL)
    included = sort_gene_names([node for node in S.nodes])
    pa.rename({'Isolate': "Genome_ID"}, axis=1, inplace=True)
    pa.set_index('Genome_ID')[included].to_csv(outfile_csv, sep=',')
    print(Fore.GREEN + "Saved csv file to", outfile_csv + Style.RESET_ALL)

    print(Fore.YELLOW + "Saving subgraph to graphml. . ." + Style.RESET_ALL)
    nx.readwrite.graphml.write_graphml(S, outfile_graphml)
    print(Fore.GREEN + "Saved graphml file to {}. This file may be opened in Cytoscape.".format(outfile_graphml) + Style.RESET_ALL)

    print(Fore.YELLOW + "Running R script to generate heatmap. This may take a moment." + Style.RESET_ALL)
    command = "Rscript {0}/PA_heatmaps_script.R {0}/data/core_gene_alignment.aln.treefile {0}/{1} {0}/data/isolate_color.csv {2}".format(cwd, escape_brackets(outfile_csv), escape_brackets(heatmapfile))
    status = subprocess.call(command, shell=True)
    if status == 0:
        print(Fore.GREEN + "Saved heatmap image file to {}".format(heatmapfile) + Style.RESET_ALL)
    else:
        print(Fore.RED + "There was a problem generating the heatmap." + Style.RESET_ALL)
        print("Command used: {0}".format(command))

    print(Fore.YELLOW + "Making distribution plot. . ." + Style.RESET_ALL)
    sns.set_context('talk')
    sns.set_palette('Paired')
    palette = ['C3', 'C5', 'C8', 'C1', 'C0']
    short_habs = {'Agriculture': 'AGRI',
                  'Clinical': "CLIN",
                  'Natural Water': "NW",
                  'Wastewater Agr.': 'AWW',
                  'Wastewater Human': 'MWW'}
    #list of genomes where query is present
    series = pa.set_index('Genome_ID').dropna()[gene]
    present = list(series[series != 0].index)
    #count presences per habitat
    t = meta.loc[present].reset_index()[['Habitat', 'Isolate']].groupby('Habitat',dropna=False).count()
    t.reset_index(inplace=True)
    for hab in meta['Habitat'].unique():
        if hab not in t['Habitat'].unique():
            t = t.concat({'Habitat': hab, 'Isolate':0}, ignore_index=True)
    t['Habitat'] = [short_habs[h] for h in t['Habitat']]

    padding = math.ceil(t['Isolate'].max() /100 * 10)
    plt.ylim((0, t['Isolate'].max() + padding))
    t.rename({'Isolate': '# Genomes Present'},axis=1, inplace=True)
    g = sns.barplot(data=t, x='Habitat', y='# Genomes Present', palette=palette )
    g.set_title('Frequency of {} by habitat'.format(gene))
    g.bar_label(g.containers[0], label_type='edge')
    plt.savefig(outfile_habitats, bbox_inches='tight')
    print(Fore.GREEN + "Saved distribution bar plot to {0}".format(outfile_habitats) + Style.RESET_ALL)
    print("Done!")
