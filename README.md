# monocle_tools

Convert a Seurat object to CDS for use in monocle V3 processing.

Preserves UMAP and PCA embedding and clustering

Commented out lines plot Seurat and CDS versions to inspect - this could be a option, but now requires manual editing.

At some point, Monocle and/or Seurat may support this directly... Until then.


Some operational notes:
* This is not robust to changes in Monocle, which was not out of beta when I wrote it.
* If you subset in Monocle, it's fast, but it breaks. I don't know why. Subset the SO, then reconvert.
