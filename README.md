# ColorMyCells

## A Biological Approach to Cell Type Visualization

`colormycells` is a specialized Python package that solves a common problem in single-cell analysis: **creating colormaps where the perceptual distance between colors meaningfully represents the biological similarity between cell types**.

## The Problem: Color Selection in Single-Cell Visualization

Standard colormaps like `tab10` or `rainbow` aren't great when applied to single-cell data:

- **Biological meaning is lost**: Default colormaps assign colors arbitrarily, with no relation to cell type similarity
- **Perceptual imbalance**: Some colors jump out while others blend together, drawing attention to cell types for no biological reason
- **Limited palette**: Most standard colormaps support 10-20 colors, but you may have hundreds of cell types

<table>
  <tr>
    <td><img src="imgs/tab10.png" width="200"/></td>
    <td><img src="imgs/color_my_cells.png" width="200"/></td>
  </tr>
  <tr>
    <td><center>Standard Colormap</center></td>
    <td><center>ColorMyCells</center></td>
  </tr>
</table>


## Installation

```bash
pip install colormycells
```

For 3D visualization support (optional):

```bash
pip install colormycells[full]
```

## Dependencies

- numpy
- pandas
- matplotlib
- seaborn
- scikit-learn
- scipy
- colour-science

### Optional Dependencies (for 3D visualization)
- pillow
- ipython

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

# Visualize the color space with 2D and 3D plots
colors = get_colormap(adata, key="cell_type", plot_colorspace=True)
# A 2D scatter plot will be displayed, and in Jupyter notebooks
# an interactive 3D rotating visualization will also be shown

# Create reproducible colormap with seed
colors1 = get_colormap(adata, key="cell_type", seed=123)
colors2 = get_colormap(adata, key="cell_type", seed=123)  # colors1 and colors2 will be identical
```

You can also pass a file path directly:

```python
# Load directly from file
colors = get_colormap("your_data.h5ad", key="cell_type")

# Works with various file formats
colors = get_colormap("expression_matrix.csv", key="cell_type")
```

## Parameters

- **adata**: AnnData object with observations/cells as rows and variables/genes as columns
- **key**: Key in adata.obs encoding cell type information (default: "cell_type")
- **plot_colorspace**: Whether to visualize the colorspace (default: False)
- **include_unknown**: Whether to include "Unknown" category in the colormap (default: False)
- **unknown_color**: Color to use for "Unknown" category if not included (default: 'w')
- **deficiency**: Type of color vision deficiency to simulate (options: None, "Deuteranomaly", "Protanomaly", "Tritanomaly", default: None)
- **severity**: Severity of color vision deficiency (0-100, default: 0)
- **seed**: Random seed for reproducible colormaps (default: 42, set to None for stochastic behavior)


## Our Approach: Biology-Driven Color Assignment

`colormycells` takes a fundamentally different approach:

1. **Biological similarity drives color selection**: Similar cell types receive similar colors
2. **Gene expression determines color**: We use the average expression profile of each cell type (pseudobulk) to measure cell type similarity
3. **Perceptually uniform color space**: We map cell type relationships to the LUV color space, where perceptual distances are uniform
4. **Intuitive visualization**: The result is a colormap where visual intuition aligns with biological reality

The result is a colormap where:
- Similar cell types appear in similar colors
- Color distances reflect biological relationships
- Visualizations become more intuitive to interpret

![Description](imgs/cell_types_3d.gif)


## Note

Color vision deficiency simulation is currently not fully implemented.

## License

MIT