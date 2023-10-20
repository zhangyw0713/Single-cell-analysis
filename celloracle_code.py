import celloracle as co
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns



co.check_python_requirements()
#anndt = co.data_conversion.seurat_object_to_anndata("/cluster/home/zyw_jh/projects/mef2d/cellchat/cellchat_TCF3.rds",delete_tmp_file=True)

adata = sc.read_h5ad("/cluster/home/zyw_jh/projects/mef2d/celloracle/recluster.hdnew.h5ad") #change to newone
print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}")

n_cells_downsample = 30000
if adata.shape[0] > n_cells_downsample:
    # Let's dowmsample into 30K cells
    sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)

print(f"Cell number is :{adata.shape[0]}")
base_GRN = co.data.load_human_promoter_base_GRN(version="hg38_gimmemotifsv5_fpr2")

oracle = co.Oracle()
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

sc.pl.pca(adata, color='CST3')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.draw_graph(adata)
sc.tl.tsne(adata)
adata.layers['counts'] = adata.raw.X.copy()

adata.X = adata.layers["counts"].copy()
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="cell_types1",
                                   embedding_name="X_draw_graph_fa")

oracle.import_TF_data(TF_info_matrix=base_GRN)


oracle.perform_PCA()
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,b_maxl=k*4, n_jobs=4)

oracle.to_hdf5("/cluster/home/zyw_jh/projects/mef2d/celloracle/hdnew.celloracle.oracle")
#oracle = co.load_hdf5("/cluster/home/zyw_jh/projects/mef2d/celloracle/mef2d.celloracle.oracle")

#sc.pl.draw_graph(oracle.adata, color="cell_types_broad",save="fa2_hd_crop.pdf")
links = oracle.get_links(cluster_name_for_GRN_unit="cell_types1", alpha=10, verbose_level=10)
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)


links.get_network_score()

links.to_hdf5(file_path="/cluster/home/zyw_jh/projects/mef2d/celloracle/links.hd_crop1.celloracle.links")
#links = co.load_hdf5(file_path="/cluster/home/zyw_jh/projects/mef2d/celloracle/links.hd_crop1.celloracle.links")

links.plot_score_comparison_2D(value="degree_centrality_all",
                               cluster1="proB", cluster2="preB_I",
                               percentile=98, save="/cluster/home/zyw_jh/figures/hd_pro_preI_score_comparison")
links.plot_scores_as_rank(cluster="preB_I", n_gene=30, save="/cluster/home/zyw_jh/figures/hd_ranked_score")



oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)
goi = "EGR1"
sc.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis",save="EGR1.HD_crop1.exp.pdf")

oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)
                      
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

oracle.calculate_embedding_shift(sigma_corr=0.05)                     
                      
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

#tune scale to make vector more visible
scale = 20
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")


plt.savefig("/cluster/home/zyw_jh/figures/hd_quiver20.pdf")
#plt.show()
                      
n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)                   
                      
min_mass = 0.01
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)                    
                      
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale_simulation = 0.5
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
plt.savefig("/cluster/home/zyw_jh/figures/FOS.simulationflow.pdf")
#plt.show()


fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=10)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)



###pseudotime

pt = Pseudotime_calculator(adata=adata,
                           obsm_key="X_draw_graph_fa", # Dimensional reduction data name
                           cluster_column_name="cell_types_broad" # Clustering data name
                           )
print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)
pt.plot_cluster(fontsize=8)

clusters_in_B_lineage = ['CLP', 'pre_proB', 'proB', 'preB_I', 'preB_II', 'immatureB',
                          'matureB']


# Make a dictionary
lineage_dictionary = {"Lineage_B": clusters_in_B_lineage}

# Input lineage information into pseudotime object
pt.set_lineage(lineage_dictionary=lineage_dictionary)

# Visualize lineage information
pt.plot_lineages()



fig, ax = plt.subplots(figsize=[6,6])

sc.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name, ax=ax, cmap="rainbow",
                color=["Pseudotime"],save="pseudotime.pdf")


