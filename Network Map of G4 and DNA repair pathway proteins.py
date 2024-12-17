# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:20:47 2024

@author: Sudipta Lahiri
"""

import networkx as nx
import matplotlib.pyplot as plt
import community as community_louvain

# New dataset for DNA repair mechanisms and associated proteins
data = {
    'HR': {'ATM', 'ATR', 'PARP1', 'BLM', 'EXOI', 'WRN', 'FANCD', 'FANCJ', 'MRE11', 'RAD50', 'NBS1', 'RPA', 'BRCA1', 'BRCA2', 'RAD51', 'POLδ', 'LIGASE I'},
    'NHEJ': {'ARTEMIS', 'ATM', 'ATR', 'PARP1', 'XLF', 'XRCC4', 'KU70/80', 'DNAPKCs', 'POLµ', 'LIGASE IV'},
    'MMR': {'PCNA', 'RFC', 'MSH2', 'MSH3', 'MSH6', 'PMS1', 'PMS2', 'MLH1', 'MLH3', 'RPA', 'EXO1', 'POLδ', 'LIGASE I'},
    'BER': {'PCNA', 'PARP1', 'PARP2', 'OGG1', 'APE1', 'XRCC1', 'DNA2', 'FEN1', 'LIGASE I', 'LIGASE III', 'POLβ', 'NEIL1', 'NEIL2', 'NEIL3'},
    'NER': {'PCNA', 'XPA', 'XPB', 'XPC', 'XPD', 'XPE', 'XPG', 'XPF', 'XRCC1', 'ERCC1', 'POLδ', 'POLε', 'LIGASE I', 'RPA', 'TFIIH'},
    'G4': {'PARP1', 'BLM', 'WRN', 'PCNA', 'MRE11', 'RAD50', 'NBS1', 'RPA', 'FANCJ', 'BRCA1', 'BRCA2', 'RAD51', 'DNAPKCs', 'MSH2', 'MLH1', 'OGG1', 'APE1', 'NEIL1', 'NEIL2', 'NEIL3', 'XPB', 'XPD'}
}

# Step 1: Identify common proteins across all mechanisms
common_proteins = set.intersection(*data.values())

# Step 2: Create a bipartite graph where proteins and repair mechanisms are nodes
G = nx.Graph()

# Add nodes for proteins and repair mechanisms
for mechanism, proteins in data.items():
    for protein in proteins:
        G.add_node(protein, bipartite=0)  # Proteins as one set of nodes
        G.add_node(mechanism, bipartite=1)  # Repair mechanisms as another set
        G.add_edge(protein, mechanism)  # Create an edge between the protein and its mechanism

# Step 3: Apply Louvain community detection to the protein nodes
partition = community_louvain.best_partition(G)

# Step 4: Visualize the clusters of proteins
plt.figure(figsize=(15, 15))

# Create a layout for visualization, adjusting k and iterations for better spacing
pos = nx.spring_layout(G, k=0.5, iterations=200)  # Increased iterations and adjusted k

# Color nodes based on the community (cluster) they belong to, and identify common proteins
node_colors = []
node_sizes = []  # To store node sizes for each node
edges_color = []  # To store edge colors
edges_width = []  # To store edge widths
labels = {}  # Dictionary to store labels for specific nodes only

for node in G.nodes():
    if node in common_proteins:
        node_colors.append('blue')  # Common proteins are colored blue
        node_sizes.append(3000)  # Larger size for common proteins
        labels[node] = node  # Label common proteins
    elif node in data['HR']:
        node_colors.append('red')  # HR proteins are colored red
        node_sizes.append(2500)  # Larger size for HR proteins
        labels[node] = node  # Label HR proteins
    elif node in data['NHEJ']:
        node_colors.append('green')  # NHEJ proteins are colored green
        node_sizes.append(2500)  # Larger size for NHEJ proteins
        labels[node] = node  # Label NHEJ proteins
    elif node in data['MMR']:
        node_colors.append('pink')  # MMR proteins are colored pink
        node_sizes.append(2500)  # Larger size for MMR proteins
        labels[node] = node  # Label MMR proteins
    elif node in data['BER']:
        node_colors.append('lightblue')  # BER proteins are colored light blue
        node_sizes.append(2500)  # Larger size for BER proteins
        labels[node] = node  # Label BER proteins
    elif node in data['NER']:
        node_colors.append('orange')  # NER proteins are colored orange
        node_sizes.append(2500)  # Larger size for NER proteins
        labels[node] = node  # Label NER proteins
    elif node in data['G4']:
        node_colors.append('purple')  # G4 proteins are colored purple
        node_sizes.append(2500)  # Larger size for G4 proteins
        labels[node] = node  # Label G4 proteins
    else:
        node_colors.append('gray')  # Other proteins are colored gray (different communities)
        node_sizes.append(5000)  # Default size for other proteins

# Step 5: Assign edge color and width based on G4 interactions
for u, v in G.edges():
    if u in data['G4'] or v in data['G4']:
        edges_color.append('black')  # Black color for edges involving G4 proteins
        edges_width.append(4)  # Thicker edges
    else:
        edges_color.append('gray')  # Default color for other edges
        edges_width.append(1)  # Default width for other edges

# Step 6: Visualize the network with the updated font size for specific proteins
nx.draw(G, pos, with_labels=True, labels=labels, node_size=node_sizes, font_size=14, font_weight='bold', 
        node_color=node_colors, edge_color=edges_color, width=edges_width)

# Title and display
plt.title("Protein Clusters Across DNA Repair Mechanisms with G4 Interactions in Black and Thicker Edges")
plt.show()
