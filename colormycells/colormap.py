"""
Main module containing the get_colormap function for creating perceptually uniform colormaps.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Literal, Dict, Any
from sklearn import manifold
from scipy.stats import special_ortho_group
import colour


def get_colormap(adata, key="cell_type", plot_colorspace=False, include_unknown=False, unknown_color='w',
                deficiency: Literal[None, "Deuteranomaly", "Protanomaly", "Tritanomaly"] = None,
                severity=0):
    """ Returns a dictionary of colors in which the perceptual distance is equal to the type/type dissimilarity.
    The colormap changes each time this is run. 
    Optionally, you can specify a color deficiency (e.g. "Deuteranomaly") and severity (0-100) to create colors
    that are approximately perceptually uniform for a certain form of colorblindness.
    This uses the CVD simulator from Machado et al. (2009).
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with observations/cells as rows and variables/genes as columns
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
    Similarity is specifically 3d MDS embedding of the "pseudobulk" expression of cells in each cell type.
    What is the average gene expression across types? What is the similarity between those?
    We'll use these to define the cell-cell similarity and then select a colormap in which
    perceptual distance is equal to the type/type dissimilarity.
            
    Uses the LUV color space. Check the code to make this more brighter/vivid
    """
    labels = adata.obs[key].unique()
    if not include_unknown:
        labels = labels[labels!="Unknown"]
    bulks = []
    for label in labels:
        pseudobulk = adata[adata.obs[key]==label].X.mean(0)
        bulks.append(pseudobulk)
    bulks = np.array(np.stack(bulks))
    similarities = np.corrcoef(bulks- bulks.mean(axis=0,keepdims=True))
    if plot_colorspace:
        sns.heatmap(pd.DataFrame(similarities, index=labels, columns=labels))
        plt.show()
    embed3 = manifold.MDS(n_components=3, dissimilarity="precomputed")
    colors3 = embed3.fit_transform(1-similarities)
    random_3d_rotation = special_ortho_group.rvs(3)
    colors3 = np.matmul(colors3,random_3d_rotation)
    luv=colors3.copy()
    luv[:,0]=luv[:,0]*0.5 + .5  # squish the lightness and make it lighter
    luv[:,1:]*=2  # more vivid
    xyz = colour.Luv_to_XYZ(luv*100)
    colors_rgb = np.maximum(np.minimum(colour.XYZ_to_sRGB(xyz),1),0)
    if deficiency is not None:
        matrix = colour.blindness.matrix_cvd_Machado2009(deficiency, severity)
        # this maps normal rgb -> simulated rgb. how can we choose colors in this space?
        raise NotImplementedError("Color vision deficiency simulation not yet implemented")
    
    if plot_colorspace:
        embed = manifold.MDS(n_components=2, dissimilarity="precomputed")
        colors = embed.fit_transform(1-similarities)
        plt.scatter(colors[:,0], colors[:,1], c=colors_rgb, s=100)
        plt.gca().set_facecolor('gray')
    
    d = {cat: c for cat, c in zip(labels, colors_rgb)}
    if not include_unknown:
        d['Unknown'] = unknown_color
    return d