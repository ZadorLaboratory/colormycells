"""
Main module containing the get_colormap function for creating perceptually uniform colormaps.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Literal, Dict, Union, Optional, Any
from sklearn import manifold
from scipy.stats import special_ortho_group
import colour
import warnings


def get_colormap(
    adata: Any,
    key: str = "cell_type",
    plot_colorspace: bool = False,
    include_unknown: bool = False,
    unknown_color: str = 'w',
    deficiency: Optional[Literal["Deuteranomaly", "Protanomaly", "Tritanomaly"]] = None,
    severity: int = 0
) -> Dict[str, np.ndarray]:
    """
    Returns a dictionary of colors where perceptual distance equals cell type dissimilarity.
    
    Creates a colormap that changes each time the function is run. The function uses
    pseudobulk expression profiles to calculate similarities between cell types, then 
    maps these similarities to perceptually uniform color distances.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with observations/cells as rows and variables/genes as columns.
        Can be an actual AnnData object or a path to a file (.h5ad, .loom, etc.)
    key : str, default="cell_type"
        Key in adata.obs encoding cell type information
    plot_colorspace : bool, default=False
        Whether to visualize the colorspace
    include_unknown : bool, default=False
        Whether to include "Unknown" category in the colormap
    unknown_color : str, default='w'
        Color to use for "Unknown" category if not included
    deficiency : {None, "Deuteranomaly", "Protanomaly", "Tritanomaly"}, default=None
        Type of color vision deficiency to simulate
    severity : int, default=0
        Severity of color vision deficiency (0-100)
        
    Returns
    -------
    dict
        Dictionary mapping cell types to RGB colors
        
    Notes
    -----
    Similarity is calculated using the pseudobulk expression of cells in each cell type.
    The function:
    1. Calculates average gene expression for each cell type
    2. Computes correlations between these pseudobulk profiles
    3. Maps the similarity matrix to 3D space using MDS
    4. Transforms these coordinates to LUV color space for perceptual uniformity
    
    Uses the LUV color space with adjustments to make colors brighter/more vivid.
    
    Raises
    ------
    TypeError
        If adata is not an AnnData object or a path to a valid file
    KeyError
        If the specified key is not in adata.obs
    ValueError
        If severity is not between 0 and 100 or for other invalid parameters
    """
    # Handle different input types
    try:
        if isinstance(adata, str):
            # If adata is a path to a file
            import os
            ext = os.path.splitext(adata)[1].lower()
            if ext == '.h5ad':
                import anndata as ad
                adata = ad.read_h5ad(adata)
            elif ext == '.loom':
                import anndata as ad
                adata = ad.read_loom(adata)
            elif ext in ['.csv', '.txt', '.tsv']:
                import anndata as ad
                sep = '\t' if ext == '.tsv' else ','
                adata = ad.read_csv(adata, sep=sep).T
            else:
                raise ValueError(f"Unsupported file format: {ext}. Supported formats: .h5ad, .loom, .csv, .txt, .tsv")
    except ImportError as e:
        raise ImportError(f"Failed to import required module: {e}. Make sure anndata is installed.") from e
    except Exception as e:
        raise ValueError(f"Failed to read AnnData from file: {e}") from e
        
    # Validate inputs
    if key not in adata.obs.columns:
        raise KeyError(f"Key '{key}' not found in adata.obs. Available keys: {adata.obs.columns.tolist()}")
    
    if severity < 0 or severity > 100:
        raise ValueError("Severity must be between 0 and 100")
    
    if deficiency is not None and deficiency not in ["Deuteranomaly", "Protanomaly", "Tritanomaly"]:
        raise ValueError("Deficiency must be one of: None, 'Deuteranomaly', 'Protanomaly', 'Tritanomaly'")
    
    # Extract cell type labels
    labels = adata.obs[key].unique()
    if not include_unknown:
        labels = labels[labels != "Unknown"]
    
    if len(labels) < 2:
        raise ValueError(f"Found fewer than 2 unique labels in '{key}'. Need at least 2 cell types to create a colormap.")
    
    # Calculate pseudobulk profiles
    bulks = []
    for label in labels:
        try:
            cells = adata[adata.obs[key] == label]
            if cells.shape[0] == 0:
                warnings.warn(f"No cells found for type '{label}'. Skipping.")
                continue
                
            # Handle different AnnData matrix types
            if isinstance(adata.X, np.ndarray):
                pseudobulk = cells.X.mean(axis=0)
            else:
                # For sparse matrices
                pseudobulk = cells.X.mean(axis=0).A1
                
            bulks.append(pseudobulk)
        except Exception as e:
            raise RuntimeError(f"Error calculating pseudobulk for cell type '{label}': {e}") from e
    
    if len(bulks) < 2:
        raise ValueError("Could not calculate pseudobulk for at least 2 cell types. Check your data.")
    
    # Create similarity matrix
    bulks = np.array(np.stack(bulks))
    
    # Center the data for better correlation
    bulks_centered = bulks - bulks.mean(axis=0, keepdims=True)
    
    # Calculate correlation coefficients
    try:
        similarities = np.corrcoef(bulks_centered)
    except np.linalg.LinAlgError as e:
        raise RuntimeError(f"Failed to compute correlation matrix: {e}") from e
    
    # Plot similarity matrix if requested
    if plot_colorspace:
        plt.figure(figsize=(10, 8))
        sns.heatmap(pd.DataFrame(similarities, index=labels, columns=labels))
        plt.title(f"Cell Type Similarity Matrix (based on {key})")
        plt.tight_layout()
        plt.show()
    
    # Create dissimilarity matrix for MDS
    dissimilarities = 1 - similarities
    
    # Handle any numerical issues in the dissimilarity matrix
    # Ensure the matrix is symmetric
    dissimilarities = (dissimilarities + dissimilarities.T) / 2
    
    # Ensure values are in valid range [0, 1]
    dissimilarities = np.clip(dissimilarities, 0, 1)
    
    # Fix any NaN values
    np.fill_diagonal(dissimilarities, 0)
    
    # Create 3D embedding with MDS
    try:
        embed3 = manifold.MDS(
            n_components=3, 
            dissimilarity="precomputed",
            random_state=42,  # For reproducibility
            n_init=10  # Multiple initializations for better results
        )
        colors3 = embed3.fit_transform(dissimilarities)
    except Exception as e:
        raise RuntimeError(f"MDS embedding failed: {e}") from e
    
    # Apply random rotation for variety in colors
    random_3d_rotation = special_ortho_group.rvs(3)
    colors3 = np.matmul(colors3, random_3d_rotation)
    
    # Scale to fit in LUV color space
    # Normalize to [-1, 1] range
    min_vals = colors3.min(axis=0, keepdims=True)
    max_vals = colors3.max(axis=0, keepdims=True)
    range_vals = max_vals - min_vals
    # Avoid division by zero
    range_vals[range_vals == 0] = 1
    colors3_norm = 2 * (colors3 - min_vals) / range_vals - 1
    
    # Convert to LUV color space
    luv = colors3_norm.copy()
    luv[:, 0] = luv[:, 0] * 0.5 + 0.5  # Adjust lightness (L) to be 0-1
    luv[:, 1:] *= 2  # Increase chromaticity (u,v) for more vivid colors
    
    # Convert LUV to XYZ to RGB
    try:
        xyz = colour.Luv_to_XYZ(luv * 100)  # Scale to 0-100 for colour package
        colors_rgb = np.maximum(np.minimum(colour.XYZ_to_sRGB(xyz), 1), 0)
    except Exception as e:
        raise RuntimeError(f"Color space conversion failed: {e}") from e
    
    # Handle color vision deficiency
    if deficiency is not None:
        warnings.warn(
            f"Color vision deficiency simulation ('{deficiency}') is not fully implemented. "
            "Returning regular colormap instead."
        )
    
    # Visualize the 2D color space if requested
    if plot_colorspace:
        try:
            # Create 2D embedding with MDS for visualization
            embed2 = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=42)
            colors2d = embed2.fit_transform(dissimilarities)
            
            plt.figure(figsize=(10, 8))
            scatter = plt.scatter(
                colors2d[:, 0], 
                colors2d[:, 1], 
                c=colors_rgb, 
                s=100, 
                edgecolor='black'
            )
            
            # Add labels
            for i, label in enumerate(labels):
                plt.annotate(
                    label, 
                    (colors2d[i, 0], colors2d[i, 1]),
                    textcoords="offset points",
                    xytext=(0, 10),
                    ha='center'
                )
                
            plt.gca().set_facecolor('lightgray')
            plt.title("2D Embedding of Cell Types with Assigned Colors")
            plt.tight_layout()
            plt.show()
        except Exception as e:
            warnings.warn(f"Could not create 2D visualization: {e}")
    
    # Create the final colormap dictionary
    valid_labels = labels[:len(colors_rgb)]  # In case we skipped some labels
    colormap = {cat: c for cat, c in zip(valid_labels, colors_rgb)}
    
    # Add 'Unknown' category if requested
    if include_unknown and 'Unknown' in adata.obs[key].unique():
        if 'Unknown' not in colormap:
            colormap['Unknown'] = np.array([1.0, 1.0, 1.0])  # white
    elif not include_unknown:
        colormap['Unknown'] = unknown_color if isinstance(unknown_color, np.ndarray) else np.array([1.0, 1.0, 1.0])
    
    return colormap