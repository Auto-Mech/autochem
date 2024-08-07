""" Draw 2D images
"""

from ..extern import rdkit_
from ._conv import smiles


def draw(chi, save_path=None):
    """ Draw an image
    """
    smi = smiles(chi, res_stereo=False)
    rdm = rdkit_.from_smiles(smi)
    img = rdkit_.draw(rdm)
    if save_path is not None:
        img.save(save_path)

    return img


def draw_grid(chis, legends=None,
              img_per_row=3, sub_img_size=(200, 200),
              save_path=None):
    """ Draw an image
    """
    smis = [smiles(chi, res_stereo=False) for chi in chis]
    rdms = tuple(rdkit_.from_smiles(smi) for smi in smis)
    img = rdkit_.draw_grid(
        rdms,
        img_per_row=img_per_row,
        sub_img_size=sub_img_size,
        legends=legends)
    if save_path is not None:
        img.save(save_path)

    return img
