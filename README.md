# SV-MONA
SV-MONA (Single-Variable Multi-Omics Association and Network Analysis) is an R framework for integrating a single variable (e.g., environmental exposure, clinical outcome, biomarker, or phenotype) with one or more omics datasets. It identifies associated molecular features, constructs integrative correlation networks, detects network communities, evaluates network topology, and exports results for downstream analysis and Cytoscape visualization.

Features
Single-variable multi-omics association analysis. Unlike conventional univariate association analyses that produce independent feature lists or star-shaped networks, SV-MONA reconstructs an integrative multi-omics network among associated features, enabling community detection, network topology analysis, and systems-level biological interpretation.
Correlation network construction. Network-level statistics are also calculated, including network density, number of connected components, global clustering coefficient, average degree, average path length, network diameter, total nodes and edges.
Community detection and network topology analysis. Louvain community detection is applied to identify highly interconnected subnetworks among CPITN-associated omics features. Potential communities may represent microbial-metabolic interactions, protein-metabolite modules, coordinated periodontal severity signatures, candidate exposure-response pathways. These communities provide biologically interpretable systems-level information that is not obtainable from feature-wise association analyses alone.
Hub node identification. Hub nodes are identified based on network topology measures. Examples include High degree, High betweenness, High eigenvector centrality, High PageRank. Hub nodes may represent key microbial taxa, critical metabolites, important proteins, candidate biomarkers of periodontal severity. Importantly, hub identification is performed after integrating both CPITN-feature and omics-feature associations, thereby allowing clinically relevant molecular modules to emerge.
Computational optimizations. Several computational optimizations is incorporated to improve the analysis of high-dimensional multi-omics datasets. It performs blockwise pairwise correlation analysis to efficiently process large datasets while minimizing memory usage. Memory-efficient feature selection is applied prior to network construction to reduce computational overhead without compromising network integrity. The framework automatically calculates node- and network-level topology metrics and exports networks in GraphML and Cytoscape-compatible formats for downstream visualization and analysis. In addition, SV-MONA generates publication-ready network figures in PDF and PNG formats.
Applicable to metabolomics, proteomics, microbiome, transcriptomics, exposomics, and other high-dimensional datasets
Requirements
R (≥ 4.2 recommended)


If you use SV-MONA in your research, please cite this GitHub repository and the associated publication (when available).

Author

Xiaojia He, Ph.D.
UL Research Institutes' Chemical Insights
