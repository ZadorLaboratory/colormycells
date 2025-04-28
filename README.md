# ColorMyCells

A Python package for creating perceptually uniform colormaps for cell types in single-cell data analysis.

## Description

`colormycells` provides a way to create colormaps where the perceptual distance between colors is equal to the biological similarity between cell types. This makes visualizations more intuitive, as similar cell types appear in similar colors.

The package is designed to work with AnnData objects, commonly used in single-cell analysis with packages like Scanpy, Seurat, or AnnData.

## Installation

For local development, clone the repository and install the package in editable mode:
```bash
git clone https://github.com/ZadorLaboratory/colormycells.git
cd colormycells
pip install -e .
```

## Dependencies

- numpy
- pandas
- matplotlib
- seaborn
- scikit-learn
- scipy
- colour-science

## Usage

```python
import scanpy as sc
from colormycells import get_colormap

# Load your data
adata = sc.read_h5ad("your_data.h5ad")

# Create a colormap based on cell type similarities
colors = get_colormap(adata, key="cell_type")

# Use the colormap for plotting
sc.pl.umap(adata, color="cell_type", palette=colors)
```

## Parameters

- **adata**: AnnData object with observations/cells as rows and variables/genes as columns
- **key**: Key in adata.obs encoding cell type information (default: "cell_type")
- **plot_colorspace**: Whether to visualize the colorspace (default: False)
- **include_unknown**: Whether to include "Unknown" category in the colormap (default: False)
- **unknown_color**: Color to use for "Unknown" category if not included (default: 'w')
- **deficiency**: Type of color vision deficiency to simulate (options: None, "Deuteranomaly", "Protanomaly", "Tritanomaly", default: None)
- **severity**: Severity of color vision deficiency (0-100, default: 0)

## How It Works

The function creates a similarity matrix between cell types based on their average gene expression profiles. It then uses Multidimensional Scaling (MDS) to embed this similarity in 3D space, which is mapped to the LUV color space. This ensures that cell types with similar expression profiles receive perceptually similar colors.

The colormap changes each time you run the function, so if you need consistent colors across multiple visualizations, save the colormap dictionary.

## Note

Color vision deficiency simulation is currently not fully implemented.

## License

GNU General Public License v3.0 (GPL-3.0)