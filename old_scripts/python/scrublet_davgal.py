###scrublet

#in R
#Make Seurat object, assume that doublets havent been filtered for umi counts....
load("r_objects/190107_Seurat_allCells_cellphoneDB.rdata")
writeMM(TA@data,file="~/Desktop/davDat.mtx")
write.table(rownames(TA@data), col.names=F,row.names=F,quote=F,file="~/Desktop/davGenes.tsv")

source activate ml_course
spyder

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

input_dir = '/Users/briglo/Desktop'
counts_matrix = scipy.io.mmread(input_dir + '/davDat.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/Davgenes.tsv', delimiter='\t', column=0))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
#still bitches about div zero                                                          n_prin_comps=30)

scrub.call_doublets(threshold=0.25)

scrub.plot_histogram();

np.savetxt(fname="/Users/briglo/Desktop/doubletScores.tsv",X=doublet_scores)
#added as meta.data

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

# # Uncomment to run tSNE - slow
# print('Running tSNE...')
# scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

# # Uncomment to run force layout - slow
# print('Running ForceAtlas2...')
# scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))
    
print('Done.')

scrub.plot_embedding('UMAP', order_points=True);

# scrub.plot_embedding('tSNE', order_points=True);
# scrub.plot_embedding('FA', order_points=True);