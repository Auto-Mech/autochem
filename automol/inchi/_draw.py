""" Draw 2D images
"""

from automol.extern import rdkit_


def draw(ich, save_path=None):
    """ Draw an image
    """
    rdm = rdkit_.from_inchi(ich)
    img = rdkit_.draw(rdm)
    if save_path is not None:
        img.save(save_path)


def draw_grid(ichs, names, save_path=None):
    """ Draw an image
    """
    rdms = tuple(rdkit_.from_inchi(ich) for ich in ichs)
    img = rdkit_.draw_grid(rdms, names=names)
    if save_path is not None:
        img.save(save_path)
