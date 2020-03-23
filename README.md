# mobclustering
Benchmarking methods for plasmid clustering


#identifying smallest mash distance from report
python find_smallest_dist.py -i mash_dist_file.txt 


# converting fastANI results to distance matrix
python process_ani_results_to_matrix.py -i fastANI_distances.txt -o fastANI.matrix.txt

# Determine cluster metrics for mash based clusters
python replicon_relaxase_concordance_clusters_mash.py

-- accepts on of the 2019-12-10-mash-clustering-* files and outputs the following fields
threshold level
mean shannon entropy
stdev shannon entropy
mean cluster sizes 
stdev cluster sizes
mean number of types within a cluster
stdev number of types within a cluster

# Determine cluster metrics for ani based clusters
python replicon_relaxase_concordance_clusters_mash.py

-- accepts on of the 2019-12-10-ani-clustering-* files and outputs the following fields
threshold level
mean shannon entropy
stdev shannon entropy
mean cluster sizes 
stdev cluster sizes
mean number of types within a cluster
stdev number of types within a cluster

# taxonomy convergence analysis
python taxonomy_convergence.py

-- Uses the Taxonomy_NCBI_Plasmids.txt file to output the following fields based on the taxa associated with each feature
accession
overall_convergence_rank
overall_convergence_name
convergence_rank_replicon
convergence_name_replicon
convergence_rank_relaxase
convergence_name_relaxase
convergence_rank_mobcluster
convergence_name_mobcluster