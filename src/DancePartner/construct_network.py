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
        The output table from a synonym function as a pandas DataFrame. Can be built manually. Requires the columns DancePartnerID, ID, Synonym, Type.
    
    
    Returns
    -------
        A network table of synonyms, IDs, types (gene product, metabolite, lipid), and the source (literature or database)
    '''

    ############################
    ## CONSTRUCT A TERM TABLE ##
    ############################

    # Make a table for term1 and term2 and bind columns
    term1_table = pd.merge(BERT_data[["term_1"]].rename({"term_1":"Synonym"}, axis = 1), synonyms, how = "left")
    term1_table = term1_table.rename({"Synonym":"Synonym1", "DancePartnerID":"ID1", "Type":"Type1"}, axis = 1)
    term2_table = pd.merge(BERT_data[["term_2"]].rename({"term_2":"Synonym"}, axis = 1), synonyms, how = "left")
    term2_table = term2_table.rename({"Synonym":"Synonym2", "DancePartnerID":"ID2", "Type":"Type2"}, axis = 1)
    term_table = pd.concat([term1_table, term2_table], axis = 1)

    # Filter out any unknowns
    term_table = term_table.dropna()
    term_table = term_table[(term_table["ID1"] != "")]
    term_table = term_table[term_table["ID2"] != ""]
    term_table["ID1"] = term_table["ID1"].astype(int)
    term_table["ID2"] = term_table["ID2"].astype(int)

    # Filter out cases where the IDs are the same
    term_table = term_table[term_table["ID1"] != term_table["ID2"]].reset_index(drop = True)
    term_table["Source"] = "literature"

    # Return table
    return term_table

def visualize_network(network_table: pd.DataFrame, 
                      gene_product_color: str = "#D55E00", 
                      metabolite_color: str = "#0072B2", 
                      lipid_color: str = "#E69F00", 
                      literature_color:str = "#56B4E9",
                      database_color: str = "#000000",
                      node_size: int = 30,
                      edge_weight: int = 4,
                      with_labels: bool = False):
    '''
    Visualize a network table    
    
    Parameters
    ----------
    network_table
        Output of build_network_table, pull_protein_protein_interactions, etc. Use pd.concat to concatenate multiple tables together. 
    
    gene_product_color
        Hexadecimal for the gene product node color. Default is #D55E00 (vermillion).
    
    metabolite_color
        Hexadecimal for the metabolite node color. Default is #0072B2 (blue).
    
    lipid_color
        Hexadecimal for the lipid node color. Default is #E69F00 (orange).
    
    literature_color
        Hexadecimal for the literature edge color. Default is #56B4E9 (skyblue).
    
    database_color
        Hexadecimal for the database edge color. Default is #000000 (black).
    
    node_size
        Size of the nodes. Default is 30. 
    
    edge_weight
        Weight of the edges. Default is 4.
    
    with_labels
        Whether labels should be included or not. Default is False. 
    
    Returns
    -------
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

def calculate_network_metrics(network: nx.Graph, metric: str = "all"):
    '''
    Calculate network metrics for the multi-omics network.    
    
    Parameters
    ----------
    network
        The output of visualize_network
    
    metric
        Either "number of components", "average component size", "degree centrality", "clustering coefficient", or "all". Default is "all". 
    
    Returns
    -------
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

def run_cooccurrence(found_terms: pd.DataFrame, synonyms: pd.DataFrame, relational_term: bool = False):
    '''
    Use co-occurrence (whether two terms appear in the same sentence) or relational_term (co-occurrence with a check to ensure a relational
    term like "binds" is in the sentence) to determine relationships. The output is a network table.

    Parameters
    ----------
    found_terms
        The output table from find_terms_in_papers() as a pandas DataFrame.

    synonyms
        The output table from a synonym function as a pandas DataFrame. Can be built manually. Requires the columns DancePartnerID, ID, Synonym, Type.

    relational_term
        A boolean to indicate whether sentences should be filtered by a relational term. Default is False. 
    
    
    Returns
    -------
        A network table of synonyms, IDs, types (gene product, metabolite, lipid), and the source (literature or database)
    '''

    ################################
    ## FILTER TO RELATIONAL TERMS ##
    ################################

    if relational_term:

        # Show the list of terms
        rel_terms = ["abate", "abated","abating","abatement","abolish","abolished","abolishing","abolition","acetylate","acetylated","acetylating","acetylation","acrylate","acrylated","acrylating","acrylation","activate",
                 "activated","activating","activation","acylate","acylated","acylating","acylation","adhere","adhered","adhering","adhesion","affix","affixed","affixing","affixion","aggregate","aggregated","aggregating",
                 "aggregation","align","aligned","aligning","alignment","alkylate","alkylated","alkylating","alkylation","anabolize","anabolized","anabolizing","anabolism","annex","annexed","annexing","annexation","append",
                 "appended","appending","appendage","assemble", "assembled", "assembling", "assembly", "associate", "associated", "associating", "association", "attach", "attached", "attaching", "attachment", "attenuate", 
                 "attenuated", "attenuating", "attenuation", "bind", "bound", "binding", "binding", "block", "blocked", "blocking", "blockage", "bridge", "bridged", "bridging", "butylate", "butylated", "butylating", 
                 "butylation", "carboxylate", "carboxylated", "carboxylating", "carboxylation", "catabolize", "catabolized", "catabolizing", "catabolism", "catalyze", "catalyzed", "catalyzing", "catalysis", "change", 
                 "changed", "changing", "changeover", "chain", "chained", "chaining", "cleave", "cleaved", "cleaving", "cleavage", "cluster", "clustered", "clustering", "cohere", "cohered", "cohering", "cohesion", 
                 "combine", "combined", "combining", "combination", "complex", "complexed", "complexing", "complexation", "confine", "confined", "confining", "confinement", "connect", "connected", "connecting", 
                 "connection", "constrain", "constrained", "constraining", "constraint", "constrict", "constricted", "constricting", "constriction", "couple", "coupled", "coupling", "create", "created", "creating", 
                 "creation", "crylate", "crylated", "crylating", "crylation", "decrease", "decreased", "decreasing", "decrement", "detach", "detached", "detaching", "detachment", "detain", "detained", "detaining", 
                 "detention", "deter", "deterred", "deterring", "deterrence", "dilute", "diluted", "diluting", "dilution", "dimerize", "dimerized", "dimerizing", "dimerization", "diminish", "diminished", "diminishing", 
                 "diminishment", "disassemble", "disassembled", "disassembling", "disassembly", "dissipate", "dissipated", "dissipating", "dissipation", "elevate", "elevated", "elevating", "elevation", "eliminate", 
                 "eliminated", "eliminating", "elimination", "enhance", "enhanced", "enhancing", "enhancement", "ethylate", "ethylated", "ethylating", "ethylation", "extenuate", "extenuated", "extenuating", "extenuation", 
                 "facilitate", "facilitated", "facilitating", "facilitation", "fasten", "fastened", "fastening", "free", "freed", "freeing", "freedom", "fuse", "fused", "fusing", "fusion", "generate", "generated", 
                 "generating", "generation", "group", "grouped", "grouping", "grouping", "glycosylate", "glycosylated", "glycosylating", "glycosylation", "hemolyze", "hemolyzed", "hemolyzing", "hemolysis", "hinder", 
                 "hindered", "hindering", "hindrance", "hydrolyze", "hydrolyzed", "hydrolyzing", "hydrolysis", "impair", "impaired", "impairing", "impairment", "impede", "impeded", "impeding", "impedance", "increase", 
                 "increased", "increasing", "increment", "induce", "induced", "inducing", "induction", "inhibit", "inhibited", "inhibiting", "inhibition", "interact", "interacted", "interacting", "interaction", 
                 "intercept", "intercepted", "intercepting", "interception", "interfere", "interfered", "interfering", "interference", "intricate", "intricated", "intricating", "intrication", "join", "joined", 
                 "joining", "joining", "liberate", "liberated", "liberating", "liberation", "ligate", "ligated", "ligating", "ligation", "link", "linked", "linking", "linkage", "loosen", "loosened", "loosening", 
                 "metabolize", "metabolized", "metabolizing", "metabolism", "methylate", "methylated", "methylating", "methylation", "moderate", "moderated", "moderating", "moderation", "modulate", "modulated", 
                 "modulating", "modulation", "neutralize", "neutralized", "neutralizing", "neutralization", "obstruct", "obstructed", "obstructing", "obstruction", "occlude", "occluded", "occluding", "occlusion", 
                 "oligomerize", "oligomerized", "oligomerizing", "oligomerization", "organize", "organized", "organizing", "organization", "osmolyze", "osmolyzed", "osmolyzing", "osmolysis", "pair", "paired", 
                 "pairing", "phosphorylate", "phosphorylated", "phosphorylating", "phosphorylation", "prevent", "prevented", "preventing", "prevention", "prohibit", "prohibited", "prohibiting", "prohibition", 
                 "promote", "promoted", "promoting", "promotion", "produce", "produced", "producing", "production", "react", "reacted", "reacting", "reaction", "reduce", "reduced", "reducing", "reduction", 
                 "regulate", "regulated", "regulating", "regulation", "relate", "related", "relating", "relation", "release", "released", "releasing", "release", "repress", "repressed", "repressing", "repression", 
                 "restrict", "restricted", "restricting", "restriction", "silence", "silenced", "silencing", "slice", "sliced", "slicing", "slicing", "stimulate", "stimulated", "stimulating", "stimulation", "stop", 
                 "stopped", "stopping", "strap", "strapped", "strapping", "subdue", "subdued", "subduing", "subdual", "supplement", "supplemented", "supplementing", "supplementation", "suppress", "suppressed", 
                 "suppressing", "suppression", "tether", "tethered", "tethering", "trigger", "triggered", "triggering", "unite", "united", "uniting", "union", "ubiquitinylate", "ubiquitinylated", "ubiquitinylating", 
                 "ubiquitination", "weaken", "weakened", "weakening", "wrap", "wrapped", "wrapping", "xylosylate", "xylosylated", "xylosylating", "xylosylation", "zip", "zipped", "zipping"]
        
        # Indicate whether a row should be kept or tossed based on whether its sentence has one of the words
        found_terms = found_terms[[any([rel_term in segment for rel_term in rel_terms]) for segment in found_terms["segment"]]]

        # If there is nothing left, let the user know
        if (len(found_terms) == 0):
            print("No relational terms found. Returning None.")
            return None
    
    #####################
    ## CONSTRUCT TABLE ##
    #####################

    network_table = pd.concat([
        pd.DataFrame({"Synonym":found_terms["term_1"]}).merge(synonyms[["DancePartnerID", "Synonym", "Type"]]).rename(columns = {"Synonym": "Synonym1", "DancePartnerID": "ID1", "Type": "Type1"}),
        pd.DataFrame({"Synonym":found_terms["term_2"]}).merge(synonyms[["DancePartnerID", "Synonym", "Type"]]).rename(columns = {"Synonym": "Synonym2", "DancePartnerID": "ID2", "Type": "Type2"})
    ], axis = 1)
    network_table["Source"] = "literature"

    return(network_table)