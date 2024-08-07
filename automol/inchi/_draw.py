""" Draw 2D images
"""

from ..extern import rdkit_


def draw(ich, save_path=None, highlight_radicals=False):
    """ Draw an image
    """
    rdm = rdkit_.from_inchi(ich)
    img = rdkit_.draw(rdm, save_path, highlight_radicals=highlight_radicals)

    return img


def draw_grid(ichs, legends=None,
              img_per_row=3, sub_img_size=(200, 200),
              save_path=None):
    """ Draw an image
    """
    rdms = tuple(rdkit_.from_inchi(ich) for ich in ichs)
    img = rdkit_.draw_grid(
        rdms,
        img_per_row=img_per_row,
        sub_img_size=sub_img_size,
        legends=legends)
    if save_path is not None:
        img.save(save_path)

    return img
