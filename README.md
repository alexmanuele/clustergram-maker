# clustergram-maker
Tailor-made solution for extracting some plots from a specific dataset. 

This tool contains a fully connected graph of genes in the data files, as well as presence absence data of those genes in an annotated set of genomes.
The tool will ask you for a focal gene around which to build a network, along with thresholds to filter graph edges. It will compute the local subgraph.

## Installation
Follow the instructions at `install_instructions.txt`

## Usage
### Parameter requirements
You will be asked to provide several parameters:
- `Gene name`. Must exactly match a gene that exists in the network. Using the interactive mode, suggestions will be made if your input doesn't match.
- `p` threshold. Edges have p-values between 0 and 1. The p threshold specifies the *upper* bound of p value and will remove all edges with p > this bound.
- `LR` threshold. Edges have numeric LR values. The LR theshold specifies the *lower* bound of LR values and will remove edges with LR < the threshold.
- `d`. Must be an integer. The number of neighbors away from the focal node. Providing a value of `0` will cause the program to return a graph containing only edges incident on the focal node.

### Option 1: Interactive
From the project directory, simply run `python extract_heatmaps.py`. An interactive terminal will guide you through the process.

### Option 2: CSV file
If you wish to quickly run several graphs, you may provide your parameters in a csv file. The file format is highly inflexible. It must use a comma as the delmiter. It must contain exactly four columns; `gene, p, lr, d`, in that order. 
Each row in the file will correspond to a new run. The parameters of each row will be validated individually. Console messages will indicate whether a row was successful or not.

To use csv mode, run `python extract_heatmap.py path/to/input.csv`. An example is available in the `data` folder.

## Output
For each gene and set of parameters, the program will create a directory with several results files.
The program will create a `results` directory if one doesn't exist. It will also create subdirectories named for the focal genes of each run.
Within the gene-named directories, the program will create a directory named for the parameters used. For example, if a user runs `gene=vanA, p=0.05, lr=100, d=1`, the program will create: 
`results/vanA/p0.05_lr100.0_d1/`

Within the result directory, the program will save 4 files:
- `gene_habitat_dist_proportion.png`. This file shows the number of genomes of each habitat containing the focal gene. On the y axis, the proportion of total genomes which the gene appears in is shown. e.g. `(n genomes of habitat Y / n genomes gene appears in)`.
- `gene_habitat_dist.png`. Similar to above, but y axis shows genome counts
- `gene_PA_heatmap.png`. A clustergram showing the presence/absence pattern of all genes in the computed subgraph
- `gene_neighborhood_graph.graphml`. A graphml file containing the computed subgraph. May be loaded into cytoscape.

## Extra
The program will generate a `.graphml` file which may be loaded into Cytoscape. We have also prepared a Cytoscape Style File, `data/NicheGraph.xml`. 
To use this, load the `.graphml` file into Cytoscape using `file -> import -> network from file`. Then, import the style file using `file -> import -> style from file`. You may then select "Niche Style" from the styles tab.
 
