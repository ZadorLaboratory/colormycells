"""
Module for generating perceptually uniform colormaps based on cell type similarity.

This module provides tools for creating colormaps where the perceptual distance
between colors matches the biological similarity between cell types.
"""

from __future__ import annotations

import io
import warnings
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any, Callable, Dict, Literal, Optional, Tuple, Union, cast

import colour
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.animation import FuncAnimation
from sklearn import manifold
from scipy.stats import special_ortho_group


# Type aliases for clarity
AnnDataLike = Any  # Could be anndata.AnnData or a string path
ColorDict = Dict[str, np.ndarray]
ColorDeficiency = Literal["Deuteranomaly", "Protanomaly", "Tritanomaly"]
DataMatrix = np.ndarray
LabelVector = np.ndarray


@dataclass
class ColorParameters:
    """Parameters for color generation and visualization."""
    
    key: str
    plot_colorspace: bool
    include_unknown: bool
    unknown_color: Union[str, np.ndarray]
    deficiency: Optional[ColorDeficiency]
    severity: int


def get_colormap(
    adata: AnnDataLike,
    key: str = "cell_type",
    plot_colorspace: bool = False,
    include_unknown: bool = False,
    unknown_color: str = 'w',
    deficiency: Optional[ColorDeficiency] = None,
    severity: int = 0
) -> ColorDict:
    """
    Generate a colormap where perceptual distance equals cell type dissimilarity.
    
    Creates a colormap that changes each time the function is run. The function uses
    pseudobulk expression profiles to calculate similarities between cell types, then 
    maps these similarities to perceptually uniform color distances.
    
    Parameters
    ----------
    adata : AnnData or str
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
    deficiency : {"Deuteranomaly", "Protanomaly", "Tritanomaly"} or None, default=None
        Type of color vision deficiency to simulate
    severity : int, default=0
        Severity of color vision deficiency (0-100)
        
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping cell types to RGB colors
        
    Notes
    -----
    Similarity is calculated using the pseudobulk expression of cells in each cell type.
    The function:
    1. Calculates average gene expression for each cell type
    2. Computes correlations between these pseudobulk profiles
    3. Maps the similarity matrix to 3D space using MDS
    4. Transforms these coordinates to LUV color space for perceptual uniformity
    
    Examples
    --------
    >>> import scanpy as sc
    >>> from colormycells import get_colormap
    >>> adata = sc.read_h5ad("my_data.h5ad")
    >>> colors = get_colormap(adata, key="cell_type")
    >>> sc.pl.umap(adata, color="cell_type", palette=colors)
    
    >>> # With visualization
    >>> colors = get_colormap(adata, key="cell_type", plot_colorspace=True)
    
    >>> # Loading directly from file
    >>> colors = get_colormap("my_data.h5ad", key="cell_type")
    """
    # Validate parameters
    _validate_parameters(severity, deficiency)
    
    # Create a parameter object to pass around
    params = ColorParameters(
        key=key,
        plot_colorspace=plot_colorspace,
        include_unknown=include_unknown,
        unknown_color=unknown_color,
        deficiency=deficiency,
        severity=severity
    )
    
    # Process data
    adata_obj = _load_adata(adata)
    labels, bulks = _calculate_pseudobulks(adata_obj, params)
    similarities, valid_labels = _compute_similarity_matrix(bulks, labels, params)
    colors_rgb, colors3d = _generate_colors(similarities)
    
    # Optional visualization
    if plot_colorspace:
        _visualize_colorspace(similarities, colors_rgb, valid_labels, colors3d)
    
    # Create and return the colormap
    return _create_colormap(valid_labels, colors_rgb, params)


def _validate_parameters(severity: int, deficiency: Optional[ColorDeficiency]) -> None:
    """
    Validate input parameters for colormap generation.
    
    Parameters
    ----------
    severity : int
        Severity of color vision deficiency
    deficiency : Optional[ColorDeficiency]
        Type of color vision deficiency
        
    Raises
    ------
    ValueError
        If parameters are invalid
    """
    if severity < 0 or severity > 100:
        raise ValueError("Severity must be between 0 and 100")
    
    if deficiency is not None and deficiency not in ["Deuteranomaly", "Protanomaly", "Tritanomaly"]:
        raise ValueError("Deficiency must be one of: None, 'Deuteranomaly', 'Protanomaly', 'Tritanomaly'")


def _load_adata(adata: AnnDataLike) -> Any:
    """
    Load AnnData from various sources.
    
    Parameters
    ----------
    adata : AnnDataLike
        AnnData object or path to file
        
    Returns
    -------
    AnnData
        Loaded AnnData object
        
    Raises
    ------
    ValueError
        If the file format is unsupported or loading fails
    ImportError
        If required dependencies are missing
    """
    if not isinstance(adata, str):
        return adata
    
    try:
        # Convert to Path for cleaner handling
        path = Path(adata)
        suffix = path.suffix.lower()
        
        if suffix == '.h5ad':
            import anndata as ad
            return ad.read_h5ad(path)
        elif suffix == '.loom':
            import anndata as ad
            return ad.read_loom(path)
        elif suffix in ['.csv', '.txt', '.tsv']:
            import anndata as ad
            sep = '\t' if suffix == '.tsv' else ','
            return ad.read_csv(path, sep=sep).T
        else:
            raise ValueError(
                f"Unsupported file format: {suffix}. "
                f"Supported formats: .h5ad, .loom, .csv, .txt, .tsv"
            )
    except ImportError as e:
        raise ImportError(
            f"Failed to import required module: {e}. "
            f"Make sure anndata is installed."
        ) from e
    except Exception as e:
        raise ValueError(f"Failed to read AnnData from file: {e}") from e


def _calculate_pseudobulks(
    adata: Any, params: ColorParameters
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate pseudobulk profiles for each cell type.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    params : ColorParameters
        Parameters for processing
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Cell type labels and their corresponding pseudobulk profiles
        
    Raises
    ------
    KeyError
        If the specified key is not found in adata.obs
    ValueError
        If fewer than 2 cell types are found or pseudobulk calculation fails
    """
    # Check if key exists
    if params.key not in adata.obs.columns:
        available_keys = adata.obs.columns.tolist()
        raise KeyError(
            f"Key '{params.key}' not found in adata.obs. "
            f"Available keys: {available_keys}"
        )
    
    # Extract unique labels
    labels = adata.obs[params.key].unique()
    if not params.include_unknown:
        labels = labels[labels != "Unknown"]
    
    if len(labels) < 2:
        raise ValueError(
            f"Found fewer than 2 unique labels in '{params.key}'. "
            f"Need at least 2 cell types to create a colormap."
        )
    
    # Calculate pseudobulk profiles
    bulks = []
    valid_labels = []
    
    for label in labels:
        try:
            cells = adata[adata.obs[params.key] == label]
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
            valid_labels.append(label)
        except Exception as e:
            raise ValueError(f"Error calculating pseudobulk for cell type '{label}': {e}") from e
    
    if len(bulks) < 2:
        raise ValueError("Could not calculate pseudobulk for at least 2 cell types. Check your data.")
    
    return np.array(valid_labels), np.array(np.stack(bulks))


def _compute_similarity_matrix(
    bulks: np.ndarray, labels: np.ndarray, params: ColorParameters
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute similarity matrix between cell types based on pseudobulk profiles.
    
    Parameters
    ----------
    bulks : np.ndarray
        Pseudobulk profiles for each cell type
    labels : np.ndarray
        Cell type labels
    params : ColorParameters
        Parameters for processing
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Similarity matrix and corresponding valid labels
        
    Raises
    ------
    RuntimeError
        If similarity calculation fails
    """
    try:
        # Center the data for better correlation
        bulks_centered = bulks - bulks.mean(axis=0, keepdims=True)
        
        # Calculate correlation coefficients
        similarities = np.corrcoef(bulks_centered)
        
        # Plot similarity matrix if requested
        if params.plot_colorspace:
            plt.figure(figsize=(10, 8))
            sns.heatmap(pd.DataFrame(similarities, index=labels, columns=labels))
            plt.title(f"Cell Type Similarity Matrix (based on {params.key})")
            plt.tight_layout()
            plt.show()
            
        return similarities, labels
    
    except np.linalg.LinAlgError as e:
        raise RuntimeError(f"Failed to compute correlation matrix: {e}") from e


def _generate_colors(similarities: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate perceptually uniform colors from similarity matrix.
    
    Parameters
    ----------
    similarities : np.ndarray
        Similarity matrix between cell types
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        RGB colors and their 3D coordinates
        
    Raises
    ------
    RuntimeError
        If color generation fails
    """
    try:
        # Create dissimilarity matrix
        dissimilarities = 1 - similarities
        
        # Ensure the matrix is symmetric and in valid range
        dissimilarities = (dissimilarities + dissimilarities.T) / 2
        dissimilarities = np.clip(dissimilarities, 0, 1)
        np.fill_diagonal(dissimilarities, 0)
        
        # Create 3D embedding with MDS
        embed3 = manifold.MDS(
            n_components=3, 
            dissimilarity="precomputed",
            random_state=42,
            n_init=10
        )
        colors3d = embed3.fit_transform(dissimilarities)
        
        # Apply random rotation for variety in colors
        random_3d_rotation = special_ortho_group.rvs(3)
        colors3d_rotated = np.matmul(colors3d, random_3d_rotation)
        
        # Scale coordinates to fit in LUV color space
        colors3d_norm = _normalize_coordinates(colors3d_rotated)
        
        # Convert to LUV color space
        luv = colors3d_norm.copy()
        luv[:, 0] = luv[:, 0] * 0.5 + 0.5  # Adjust lightness (L) to be 0-1
        luv[:, 1:] *= 2  # Increase chromaticity (u,v) for more vivid colors
        
        # Convert LUV to XYZ to RGB
        xyz = colour.Luv_to_XYZ(luv * 100)  # Scale to 0-100 for colour package
        colors_rgb = np.maximum(np.minimum(colour.XYZ_to_sRGB(xyz), 1), 0)
        
        return colors_rgb, colors3d
    
    except Exception as e:
        raise RuntimeError(f"Failed to generate colors: {e}") from e


def _normalize_coordinates(coords: np.ndarray) -> np.ndarray:
    """
    Normalize 3D coordinates to [-1,1] range.
    
    Parameters
    ----------
    coords : np.ndarray
        3D coordinates
        
    Returns
    -------
    np.ndarray
        Normalized coordinates
    """
    min_vals = coords.min(axis=0, keepdims=True)
    max_vals = coords.max(axis=0, keepdims=True)
    range_vals = max_vals - min_vals
    
    # Avoid division by zero
    range_vals[range_vals == 0] = 1
    
    # Scale to [-1, 1]
    return 2 * (coords - min_vals) / range_vals - 1


def _visualize_colorspace(
    similarities: np.ndarray,
    colors_rgb: np.ndarray,
    labels: np.ndarray,
    colors3d: np.ndarray
) -> None:
    """
    Visualize the colorspace in 2D and 3D.
    
    Parameters
    ----------
    similarities : np.ndarray
        Similarity matrix
    colors_rgb : np.ndarray
        RGB colors
    labels : np.ndarray
        Cell type labels
    colors3d : np.ndarray
        3D coordinates
    """
    try:
        # Create 2D embedding with MDS for visualization
        dissimilarities = 1 - similarities
        embed2 = manifold.MDS(
            n_components=2, 
            dissimilarity="precomputed", 
            random_state=42
        )
        colors2d = embed2.fit_transform(dissimilarities)
        
        # Plot 2D visualization
        _plot_2d_colorspace(colors2d, colors_rgb, labels)
        
        # Plot 3D visualization
        _plot_3d_colorspace(colors3d, colors_rgb, labels)
        
    except Exception as e:
        warnings.warn(f"Could not create colorspace visualization: {e}")


def _plot_2d_colorspace(
    coords: np.ndarray, colors: np.ndarray, labels: np.ndarray
) -> None:
    """
    Plot 2D visualization of the colorspace.
    
    Parameters
    ----------
    coords : np.ndarray
        2D coordinates
    colors : np.ndarray
        RGB colors
    labels : np.ndarray
        Cell type labels
    """
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(
        coords[:, 0], 
        coords[:, 1], 
        c=colors, 
        s=100, 
        edgecolor='black'
    )
    
    # Add labels
    for i, label in enumerate(labels):
        plt.annotate(
            label, 
            (coords[i, 0], coords[i, 1]),
            textcoords="offset points",
            xytext=(0, 10),
            ha='center'
        )
        
    plt.gca().set_facecolor('lightgray')
    plt.title("2D Embedding of Cell Types with Assigned Colors")
    plt.tight_layout()
    plt.show()


def _plot_3d_colorspace(
    coords: np.ndarray, colors: np.ndarray, labels: np.ndarray
) -> None:
    """
    Plot 3D visualization of the colorspace with rotation.
    
    Parameters
    ----------
    coords : np.ndarray
        3D coordinates
    colors : np.ndarray
        RGB colors
    labels : np.ndarray
        Cell type labels
    """
    try:
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib
        from IPython.display import HTML, display
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot the 3D scatter
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            coords[:, 2],
            c=colors,
            s=100,
            edgecolor='black'
        )
        
        # Add labels
        for i, label in enumerate(labels):
            ax.text(
                coords[i, 0],
                coords[i, 1],
                coords[i, 2],
                label,
                size=8,
                zorder=1,
                color='black'
            )
        
        ax.set_title("3D Embedding of Cell Types with Assigned Colors")
        ax.set_facecolor('lightgray')
        
        # Function to update the plot for animation
        def rotate(angle):
            ax.view_init(elev=30, azim=angle)
            return [scatter]
        
        # Create animation
        ani = FuncAnimation(fig, rotate, frames=range(0, 360, 5), interval=100, blit=True)
        
        # Detect environment to handle the animation display
        try:
            # Check if running in a Jupyter notebook
            get_ipython
            is_notebook = True
        except NameError:
            is_notebook = False
        
        if is_notebook:
            # Display in notebook
            plt.close()  # Prevent the static plot from displaying
            display(HTML(ani.to_jshtml()))
        else:
            # Save animation to file or show static view
            try:
                ani.save('cell_types_3d.gif', writer='pillow', fps=10)
                print("3D animation saved as 'cell_types_3d.gif'")
            except Exception as save_error:
                print(f"Could not save animation: {save_error}")
                plt.show()  # Fall back to static display
    
    except ImportError as ie:
        warnings.warn(
            f"Could not create 3D visualization: {ie}. "
            f"You may need to install additional packages."
        )
    except Exception as e:
        warnings.warn(f"Could not create 3D visualization: {e}")


def _create_colormap(
    valid_labels: np.ndarray, colors_rgb: np.ndarray, params: ColorParameters
) -> ColorDict:
    """
    Create the final colormap dictionary.
    
    Parameters
    ----------
    valid_labels : np.ndarray
        Cell type labels
    colors_rgb : np.ndarray
        RGB colors
    params : ColorParameters
        Parameters for processing
        
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping cell types to RGB colors
    """
    # Create the colormap from valid labels and colors
    colormap = {str(cat): c for cat, c in zip(valid_labels, colors_rgb)}
    
    # Add 'Unknown' category if requested
    if params.include_unknown and 'Unknown' in valid_labels:
        pass  # Already included
    elif not params.include_unknown:
        if isinstance(params.unknown_color, str):
            # Convert string to RGB if needed
            if params.unknown_color == 'w':
                colormap['Unknown'] = np.array([1.0, 1.0, 1.0])
            else:
                # Simple conversion for common color strings
                color_map = {
                    'k': np.array([0.0, 0.0, 0.0]),  # black
                    'r': np.array([1.0, 0.0, 0.0]),  # red
                    'g': np.array([0.0, 1.0, 0.0]),  # green
                    'b': np.array([0.0, 0.0, 1.0]),  # blue
                    'y': np.array([1.0, 1.0, 0.0]),  # yellow
                    'c': np.array([0.0, 1.0, 1.0]),  # cyan
                    'm': np.array([1.0, 0.0, 1.0]),  # magenta
                }
                colormap['Unknown'] = color_map.get(params.unknown_color, np.array([1.0, 1.0, 1.0]))
        else:
            colormap['Unknown'] = params.unknown_color
    
    return colormap