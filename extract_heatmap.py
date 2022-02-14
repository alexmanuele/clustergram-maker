import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import networkx as nx
#import ete3
#import dendropy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
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
    node_list = []
    if d == 0:
        node_list.append(node)
    edges = []
    for u,v,e in G.edges(*node_list, data=True):
        if e['lr'] >= lr_threshold and e['p'] <= p_threshold:
            edges.append((u,v))
    H=G.edge_subgraph(edges)
    if node in H.nodes:
        if d==0:
            return H
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
    t = string.replace('(', '_')
    t = t.replace(')', '')
    t = t.replace('/', '-')
    return t

def red_string(string):
    return Fore.RED + string + Style.RESET_ALL

def yellow_string(string):
    return Fore.YELLOW + string + Style.RESET_ALL

def green_string(string):
     return Fore.GREEN + string + Style.RESET_ALL


## user input stuff
#allow the user to exit anytime
def check_exit(in_text):
    if in_text == 'exit':
        print("You have chosen to exit.")
        sys.exit()
    return in_text

#To validate from input sheet.
def validate_param_dict(params, G):
    valid = True
    gene = params['gene']
    p = params['p']
    lr = params['lr']
    d = params['d']
    target = gene in G.nodes
    node_list = list(G.nodes)
    if not target:
        fuzzy_match = process.extract(gene, node_list, limit=3)
        print(red_string("{} not found in the network. Please double check spelling and capitalization.".format(gene)))
        print(yellow_string("Did you mean:"))
        for i in fuzzy_match:
            print(yellow_string("\t{}".format(i[0])))
        valid = False
    if p < 0.0 or p > 1.0:
        print(red_string("p value must be between 0 and 1.0 inclusive. Found {}".format(p)))
        valid = False
    try:
        lr = float(lr)
    except:
        print(red_string("LR value must be numeric. Found {}".format(lr)))
        valid = False
    try:
        d = int(d)
        assert type(d) == int
        assert d > -1
    except:
        print(red_string("Degree must be an integer >= 0. Found {}".format(d)))
        valid = False
    if valid:
        return {k:v for k, v in zip(params.keys(), [gene, p, lr, d])}
    else:
        print(red_string("Errors in param form. No results will be computed for {0} {1} {2} {3}".format(gene, p, lr, d)))
        return None

# prompts
gene_prompt = "Please enter the name of the target gene:"
lr_prompt = "Please enter the lower bound of the likelihood ratio:"
p_prompt = "Please enter the upper bound of the p-value:"
d_prompt = "Please enter the neighborhood depth (as an integer). NOTE: You may type 0 to retrieve only the edges incident on the focal node:"

cmap = {
    #Clades
    'Type A': '#841e5a',
    'Type B': '#f6c19f',
    #Sampling Location
    'Clinical Ab': '#fb9a99',
    'Clinical UK': '#e31a1c',
    'Wastewater Mun. Ab': '#cab2d6',
    'Wastewater Mun. UK': '#6a3d9a',
    'Agricultural Ab': '#b2df8a',
    'Agricultural UK': '#33a02c',
    'Natural Water Ab': '#1f78b4',
    'Wastewater Agr. Ab': '#fdbf6f',
    #Geography
    'United Kingdom': '#9bbdff',
    'Canada/Alberta': '#e36951',
    #Habitat
    'Agricultural': '#33a02c',
    'Wastewater Mun.': '#cab2d6',
    'Clinical': '#e31a1c',
    'Natural Water': '#1f78b4',
    'Wastewater Agr.': '#fdc170',
    #Missing
    np.nan: '#FFFFFF',
}



if __name__ == '__main__':
    cwd = os.getcwd()

    try:
        print(yellow_string("Loading main presence absence table. . ."))
        pa = pd.read_table('data/all_categories_pa_new.csv',sep=',' ,index_col=0)
        print(yellow_string("Loading genome label informtation. . ."))
        meta = pd.read_table('data/genome_labels.csv', sep=',', index_col=0)
        # get colors for the genomes into a dataframe
        genome_colors = meta.copy()
        genome_colors.drop('Clade', axis=1, inplace=True)
        for col in genome_colors.columns:
            genome_colors[col] = [cmap[c] for c in genome_colors[col]]

        # format a legend key for plotting
        legend_keys = []
        for i, col in enumerate(meta.columns):
            legend_keys += sorted(list(meta[col].dropna().unique()))
        legend_cmap = {k:v for k, v in cmap.items() if k in legend_keys}

        # color map for presence/absence
        pa_cmap = sns.color_palette(['#f5f5f5', '#021657'])
        print(yellow_string("Loading graph. This may take a moment. . ."))
        G = nx.graphml.read_graphml('data/pagel_results_as_network_detailed.graphml')
        print(yellow_string("Loading tree data This may take a moment. . ."))

        dm = pd.read_table('data/phylo_distance_matrix.csv', sep=',', index_col=0)
        for genome in set(dm.index) - set(pa.index):
            pa.loc[genome] = 0
        # compressed distance matrix from the ultrametric tree
        um = squareform(dm[dm.index])
        # compute linkage for clustering
        ultrameric_link = linkage(um)
        print(green_string("Data Loaded."))
    except:
        raise FileNotFoundError("Files must be located in a folder called ./data")


    #Conditional: Pass input file OR start the input loop.
    infile = None
    params = []
    if len(sys.argv) > 1:
        input_loop = False
        from_input = False
        try:
            infile = pd.read_table(sys.argv[1], sep=',')
            assert list(infile.columns) == ['gene', 'p', 'lr', 'd']
        except:
            print(list(infile.columns))
            raise ValueError("Malformed input sheet. Consult README")
        for i, row in infile.iterrows():
            pdict = validate_param_dict({**row}, G)
            if pdict:
                params.append(pdict)
        if len(params) == 0:
            print(red_string("None of the inputs were valid! Exiting..."))
            sys.exit(1)
    #loop for user input
    else:
        input_loop = True
        from_input = True

    while input_loop:
        print(green_string("Creating your data. Please answer the questions. You may type 'exit' any time to exit."))
        node_list = list(G.nodes)
        target = False
        while not target:
            gene = check_exit(input(Fore.GREEN + gene_prompt + Style.RESET_ALL))
            target = gene in G.nodes
            if not target:
                fuzzy_match = process.extract(gene, node_list, limit=3)
                print(red_string("{} not found in the network. Please double check spelling and capitalization.".format(gene)))
                print(yellow_string("Did you mean:"))
                for i in fuzzy_match:
                    print(yellow_string("\t{}".format(i[0])))

        p_valid = False
        while not p_valid:
            p = check_exit(input(green_string(p_prompt)))
            try:
                p = float(p)
                assert p >= 0.0 and p <= 1.0
                p_valid = True
            except:
                print(yellow_string("p value must be a number between 0.0 and 1.0"))

        lr_valid = False
        while not lr_valid:
            lr = check_exit(input(green_string(lr_prompt)))
            try:
                lr = float(lr)
                lr_valid = True
            except:
                print(yellow_string("LR value must be numeric."))

        degree_valid = False
        while not degree_valid:
            d = check_exit(input(green_string(d_prompt)))
            try:
                d = int(d)
                assert type(d) == int
                assert d > -1
                degree_valid = True
            except:
                print(yellow_string("Degree must be a positive integer."))

        print("Input done. Processing. . .")

        print("Computing subgraph. . .")
        print("Done. ")
        S = filter_graph(G, gene, d, lr, p)
        if len(S.nodes) == 1:
            print(red_string("The filtered graph contains only one node!"))
            print(red_string("Cannot perform clustering on singleton. Please select new parameters, or type 'exit' to exit."))
        else:
            input_loop = False
        params.append({'gene': gene, 'p': p, 'lr': lr, 'd': d})

    for pdict in params:
        gene = pdict['gene']
        p = pdict['p']
        lr = pdict['lr']
        d = pdict['d']
        if not from_input:
            S = filter_graph(G, gene, d, lr, p)
            if len(S.nodes) == 1:
                print(red_string("The filtered graph contains only one node!"))
                print(red_string("No results can be computed for {0} {1} {2} {3}".format(gene, p, lr, d)))
                continue
        print("Computing for {0} p={1} lr={2} d={3}".format(gene, p, lr, d))
        #make results directory
        results_dir ='results/{0}/p{1}_lr{2}_d{3}'.format(escape_brackets(gene), p, lr, d)
        outfile_csv = '{0}/{1}_neighborhood_PA.csv'.format(results_dir, escape_brackets(gene))
        heatmapfile = '{0}/{1}_PA_Heatmap.png'.format(results_dir, escape_brackets(gene))
        outfile_graphml = '{0}/{1}_neighborhood_graph.graphml'.format(results_dir, escape_brackets(gene))
        outfile_habitats = '{0}/{1}_habitat_dist.png'.format(results_dir, escape_brackets(gene))
        outfile_habitats2 = '{0}/{1}_habitat_dist_proportion.png'.format(results_dir, escape_brackets(gene))
        if not os.path.exists('results'):
            os.makedirs('results')
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)

        print(yellow_string("Creating presence absence table..."))
        included = sort_gene_names([node for node in S.nodes])
        #pa.rename({'Isolate': "Genome_ID"}, axis=1, inplace=True)
        #pa.set_index('Genome_ID')[included].to_csv(outfile_csv, sep=',')
        subset_pa = pa.loc[dm.index][included]
        #print(Fore.GREEN + "Saved csv file to", outfile_csv + Style.RESET_ALL)

        print(yellow_string("Saving subgraph to graphml. . ."))
        target_p = {}
        target_lr = {}
        #feature_type = {}
        target_log2_lr = {}
        target_log10_lr = {}
        for node in S.nodes:
            if node == gene:
                #feature_type #Should this not be here?
                continue
            if S.has_edge(node, gene):
                edge = S.edges[(node, gene)]
                target_p[node] = edge['p']
                target_lr[node] = edge['lr']
                target_log2_lr[node] = max(0.000000001, np.log2(edge['lr']))
                target_log10_lr[node] = max(0.000000001, np.log10(edge['lr']))
            else:
                # TODO I think it makes more sense to leave this attribute empty for these nodes
                target_p[node] = p
                target_lr[node] = lr
                target_log2_lr[node] = max(0.000000001, np.log2(lr))
                target_log10_lr[node] = max(0.000000001, np.log10(lr))
        nx.set_node_attributes(S, name='target_p', values=target_p)
        nx.set_node_attributes(S, name='target_lr', values=target_lr)
        nx.set_node_attributes(S, name='target_log2_lr', values=target_log2_lr)
        nx.set_node_attributes(S, name='target_log10_lr', values=target_log10_lr)
        nx.readwrite.graphml.write_graphml(S, outfile_graphml)
        print(green_string("Saved graphml file to {}. This file may be opened in Cytoscape.".format(outfile_graphml)))
        print(yellow_string("Creating ClusterMap . . ."))
        title = "{} Neighborhood P/A".format(gene)
        handles = [Patch(facecolor=legend_cmap[name]) for name in legend_cmap]
        g = sns.clustermap(data=subset_pa.loc[dm.index],
                           row_linkage=ultrameric_link,
                           cbar_pos=None, row_colors=genome_colors,
                           cmap=pa_cmap,
                           method='complete', metric='cityblock',
                           xticklabels=1, yticklabels=False)
        g.ax_col_dendrogram.set_title(title)
        plt.legend(handles, legend_cmap, title='Genome Labels',
                   bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')
        #plt.xticks(fontsize=6)
        plt.savefig(heatmapfile, bbox_inches='tight', dpi=300)
        print(green_string("Saved heatmap image file to {}".format(heatmapfile)))

        print(yellow_string("Making distribution plot. . ."))
        sns.set_context('talk')
        sns.set_palette('Paired')
        palette = ['C3', 'C5', 'C8', 'C1', 'C0']
        short_habs = {'Agricultural': 'AGRI',
                      'Clinical': "CLIN",
                      'Natural Water': "NW",
                      'Wastewater Agr.': 'AWW',
                      'Wastewater Mun.': 'MWW'}
        #list of genomes where query is present
        series = pa[gene]
        present = list(series[series != 0].dropna().index)
        #count presences per habitat
        t = meta.loc[present].reset_index()[['Habitat', 'Isolate']].groupby('Habitat',dropna=False).count()
        t.reset_index(inplace=True)
        for hab in meta['Habitat'].unique():
            if hab not in t['Habitat'].unique():
                t = t.append({'Habitat': hab, 'Isolate':0}, ignore_index=True)
        t['Habitat'] = [short_habs[h] for h in t['Habitat']]

        padding = math.ceil(t['Isolate'].max() /100 * 10)
        plt.figure()
        plt.ylim((0, t['Isolate'].max() + padding))
        t.rename({'Isolate': '# Genomes Present'},axis=1, inplace=True)
        g = sns.barplot(data=t, x='Habitat', y='# Genomes Present', palette=palette, order=['AGRI', 'CLIN', 'MWW', 'AWW', 'NW'] )
        g.set_title('Frequency of {} by habitat'.format(gene))
        g.bar_label(g.containers[0], label_type='edge')
        plt.savefig(outfile_habitats, bbox_inches='tight')

        plt.figure()
        summ = t['# Genomes Present'].sum()
        t['Prop. Genomes Present'] = t['# Genomes Present'].map(lambda x: x/summ)
        padding=0.1
        plt.ylim((0, t['Prop. Genomes Present'].max() + padding))
        g2 = sns.barplot(data=t, x='Habitat', y='Prop. Genomes Present', palette=palette, order=['AGRI', 'CLIN', 'MWW', 'AWW', 'NW'] )
        g2.set_title('Proportion of genomes containing {0} by habitat\n({1} total genomes)'.format(gene, summ))
        g2.bar_label(g2.containers[0], labels=[int(i) for i in g.containers[0].datavalues], label_type='edge')#labels=[int(i) for i in g.containers[0].datavalues * summ], label_type='edge',)
        plt.savefig(outfile_habitats2, bbox_inches='tight')
        print(green_string("Saved distribution bar plot to {0}".format(outfile_habitats)))
        print("Done!")
