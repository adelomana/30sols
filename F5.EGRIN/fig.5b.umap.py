import scanpy
scanpy.settings.verbosity=5

# 1. read data
adata=scanpy.read_csv('ribo.expression.data.csv')
print(adata)

# 2. associate corem labels
coremLabels=['green', 'green', 'green', 'green', 'green', 'green', 'green', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue']
adata.obs['coremLabels']=coremLabels

# 2. visualization

# 2.1. PCA
print('PCA...')
scanpy.tl.pca(adata, svd_solver='arpack')
scanpy.pl.pca(adata,save='scanpy.pca.pdf',show=False,color='coremLabels')

# 2.2. tSNE
print('tSNE...')
scanpy.tl.tsne(adata)
scanpy.pl.tsne(adata,color='coremLabels',alpha=0.5,save='scanpy.tSNE.pdf',show=False)

# 2.2. UMAP
print('UMAP...')
scanpy.pp.neighbors(adata, n_neighbors=9, n_pcs=50)
scanpy.tl.umap(adata)
scanpy.pl.umap(adata, color='coremLabels',save='scanpy.UMAP.pdf',show=False,palette=['blue','green','magenta','red'],size=1000,alpha=2/3)
