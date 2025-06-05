import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def build_network_table(BERT_data: pd.DataFrame, synonyms: pd.DataFrame):
    '''
    Build a network table of edges with biomolecule IDs and their synonyms

    Parameters
    -----------
    BERT_data 
        The output table from run_bert() as a pandas DataFrame.
    
    synonyms
        The output table from map_synonyms() as a pandas DataFrame.
    
    
    Returns
    -------
        A network table of synonms, IDs, types (gene product, metabolite, lipid), and the source (literature or database)
    '''

    ############################
    ## CONSTRUCT A TERM TABLE ##
    ############################

    # Make a table for term1 and term2 and bind columns
    term1_table = pd.merge(BERT_data[["term_1"]].rename({"term_1":"Synonym"}, axis = 1), synonyms, how = "left")
    term1_table = term1_table.rename({"Synonym":"Synonym1", "ID":"ID1", "Type":"Type1"}, axis = 1)
    term2_table = pd.merge(BERT_data[["term_2"]].rename({"term_2":"Synonym"}, axis = 1), synonyms, how = "left")
    term2_table = term2_table.rename({"Synonym":"Synonym2", "ID":"ID2", "Type":"Type2"}, axis = 1)
    term_table = pd.concat([term1_table, term2_table], axis = 1)

    # Filter out any unknowns
    term_table = term_table.dropna()
    term_table = term_table[(term_table["ID1"] != "")]
    term_table = term_table[term_table["ID2"] != ""]

    # Filter out cases where the IDs are the same
    term_table = term_table[term_table["ID1"] != term_table["ID2"]].reset_index(drop = True)
    term_table["Source"] = "literature"

    # Return table
    return term_table

def visualize_network(network_table, 
                      gene_product_color = "#D55E00", 
                      metabolite_color = "#0072B2", 
                      lipid_color = "#E69F00", 
                      literature_color = "#56B4E9",
                      database_color = "#000000",
                      node_size = 30,
                      edge_weight = 4,
                      with_labels = False):
    '''
    Visualize a network_table    
    
    Args:
        network_table (Pandas DataFrame): Output of build_network_table, pull_protein_protein_interactions, 
           etc. Use pd.concat to concatenate multiple tables together. 
        gene_product_color (String): Hexadecimal for the gene product node color. Default is #D55E00 (vermillion).
        metabolite_color (String): Hexadecimal for the metabolite node color. Default is #0072B2 (blue).
        lipid_color (String): Hexadecimal for the lipid node color. Default is #E69F00 (orange).
        literature_color (String): Hexadecimal for the literature edge color. Default is #56B4E9 (skyblue).
        database_color (String): Hexadecimal for the database edge color. Default is #000000 (black).
        node_size (Integer): Size of the nodes. Default is 30. 
        edge_weight (Integer): Weight of the edges. Default is 4.
        with_labels (Logical): Whether labels should be included or not. Default is False. 
    Returns:
        A network object and the visualization of that object
    '''

    # Reset the index on the network table
    network_table = network_table.dropna().reset_index(drop = True)

    # Detect synonyms
    detected_synonyms = pd.concat([
        network_table[["ID1", "Type1"]].rename({"ID1":"ID", "Type1":"Type"}, axis = 1),
        network_table[["ID2", "Type2"]].rename({"ID2":"ID", "Type2":"Type"}, axis = 1)
    ]).drop_duplicates()

    ###################
    ## BUILD NETWORK ##
    ###################

    # List all nodes
    nodes = network_table["ID1"].to_list()
    nodes.extend(network_table["ID2"].to_list())
    nodes = list(set(nodes))

    # List all edges
    edges = []
    for row in range(len(network_table)):
        edges.append((network_table.loc[row, "ID1"], network_table.loc[row, "ID2"]))

    # Initiate network 
    network = nx.Graph()

    # Add nodes and edges - nodes stay in the same order, but edges do not
    network.add_nodes_from(nodes)
    network.add_edges_from(edges)

    #################
    ## COLOR NODES ##
    #################

    # Create node color list
    node_color_list = []

    for node in nodes: 
        type = detected_synonyms[detected_synonyms["ID"] == node]["Type"].values[0]
        if (type == "lipid"):
            node_color_list.append(lipid_color)
        elif (type == "metabolite"):
            node_color_list.append(metabolite_color)
        else:
            node_color_list.append(gene_product_color)

    #################
    ## EDGE COLORS ##
    #################

    # Create edge color list
    edge_color_list = []

    for u,v in network.edges():
        source = network_table[network_table["ID1"].isin([u, v]) & network_table["ID2"].isin([u, v])]["Source"].tolist()[0]
        if (source == "literature"):
            edge_color_list.append(literature_color)
        else:
            edge_color_list.append(database_color)

    ####################
    ## FINISH NETWORK ##
    ####################

    # Visualize network
    nx.draw(network, node_size = node_size, node_color = node_color_list, edge_color = edge_color_list, width = edge_weight, with_labels = with_labels, font_size = 8)
    plt.legend(handles = [
        Line2D([0], [0], marker='o', color='w', label = 'Gene Product', markerfacecolor = gene_product_color, markersize=10),
        Line2D([0], [0], marker='o', color='w', label = 'Metabolite', markerfacecolor = metabolite_color, markersize=10),
        Line2D([0], [0], marker='o', color='w', label = 'Lipid', markerfacecolor = lipid_color, markersize=10),
        Line2D([0], [0], color = literature_color, label = 'Literature', markersize=10),
        Line2D([0], [0], color = database_color, label = 'Database', markersize=10),
    ])
    return network

def calculate_network_metrics(network, metric = "all"):
    '''
    Calculate network metrics for the multi-omics network.    
    
    Args:
        network (networkx object): The output of visualize network
        metric (String): Either "number of components", "average component size", "degree centrality"
            "clustering coefficient", or "all". Default is "all". 
    Returns:
        Network summary metrics
    '''

    if metric == "number of components":
        return len(list(nx.connected_components(network)))
    
    elif metric == "average component size":
        lengths = [len(x) for x in list(nx.connected_components(network))]
        return np.round(np.mean(np.array(lengths)), 4)
    
    elif metric == "degree centrality":
        return pd.DataFrame({
            "Node": list(nx.degree_centrality(network).keys()),
            "Centrality": list(nx.degree_centrality(network).values())
        }).sort_values(by = "Centrality", ascending = False).reset_index(drop = True)
    
    elif metric == "clustering coefficient":
        return np.round(np.mean(np.array(list(nx.clustering(network).values()))), 4)
    
    elif metric == "all":
        return {
            "Number of Components": calculate_network_metrics(network, "number of components"),
            "Average Component Size": calculate_network_metrics(network, "average component size"),
            "Degree Centrality": calculate_network_metrics(network, "degree centrality"),
            "Clustering Coefficient": calculate_network_metrics(network, "clustering coefficient")
        }

    else: 
        print(metric + " is not a recognized metric")