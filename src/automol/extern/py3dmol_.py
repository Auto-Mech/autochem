""" py3Dmol interface
"""

from typing import Tuple

import numpy
import py3Dmol


def create_view(image_size: int = 400) -> py3Dmol.view:
    """Create a new 3D viewer

    :param image_size: The image size, defaults to 400
    :type image_size: int, optional
    :return: A 3D viewer
    :rtype: py3Dmol.view
    """
    return py3Dmol.view(width=image_size, height=image_size)


def view_molecule_from_xyz(
    xyz_str: str, view: py3Dmol.view = None, image_size: int = 400, vib: bool = False
) -> py3Dmol.view:
    """Visualize a molecule in a 3D view, using an XYZ file.

    :param xyz_str: MOLFile block string
    :param view: An existing 3D view to append to, defaults to None
    :param image_size: The image size, if creating a new view, defaults to 400
    :param vib: Animate vibrations for the molecule?
    :return: A 3D view containing the molecule
    :rtype: py3Dmol.view
    """
    if view is None:
        view = create_view(image_size=image_size)

    options = {"vibrate": {"frames": 10, "amplitude": 1}} if vib else {}

    view.addModel(xyz_str, "xyz", options)
    view.setStyle({"stick": {}, "sphere": {"radius": 0.3}})
    view.setBackgroundColor("0xeeeeee")
    if vib:
        view.animate({"loop": "backAndForth"})

    view.zoomTo()
    return view


def view_molecule_from_xyz_trajectory(
    xyz_str: str, view: py3Dmol.view = None, image_size: int = 400
) -> py3Dmol.view:
    """Visualize a molecule in a 3D view, using an XYZ trajectory file.

    :param xyz_str: MOLFile block string
    :param view: An existing 3D view to append to, defaults to None
    :param image_size: The image size, if creating a new view, defaults to 400
    :return: A 3D view containing the molecule
    :rtype: py3Dmol.view
    """
    if view is None:
        view = create_view(image_size=image_size)

    # options = {"vibrate": {"frames": 10, "amplitude": 1}} if vib else {}

    view.addModelsAsFrames(xyz_str, "xyz")
    view.setStyle({"stick": {}, "sphere": {"radius": 0.3}})
    view.setBackgroundColor("0xeeeeee")
    view.animate({"loop": "forward"})
    view.zoomTo()
    return view


def view_molecule_from_molfile(
    mfl: str, view: py3Dmol.view = None, image_size: int = 400
) -> py3Dmol.view:
    """Visualize a molecule in a 3D view, using a MOLFile string

    :param mfl: MOLFile block string
    :type mfl: str
    :param image_size: The image size, if creating a new view, defaults to 400
    :type image_size: int, optional
    :param view: An existing 3D view to append to, defaults to None
    :type view: py3Dmol.view, optional
    :return: A 3D view containing the molecule
    :rtype: py3Dmol.view
    """
    if view is None:
        view = create_view(image_size=image_size)

    view.addModel(mfl, "sdf")
    view.setStyle({"stick": {}, "sphere": {"radius": 0.3}})
    view.setBackgroundColor("0xeeeeee")
    view.zoomTo()
    return view


def view_vector(
    vec_xyz: Tuple[float, float, float],
    orig_xyz: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    color: str = "black",
    view: py3Dmol.view = None,
    image_size: int = 400,
) -> py3Dmol.view:
    """Visualize a vector in a 3D view

    :param vec_xyz: The vector coordinates
    :type vec_xyz: Tuple[float, float, float]
    :param orig_xyz: The vector origin, defaults to (0., 0., 0.)
    :type orig_xyz: Tuple[float, float, float], optional
    :param color: The color of the vector (hex or name), default "black"
    :type color: str, optional
    :param view: An existing 3D view to append to, defaults to None
    :type view: py3Dmol.view, optional
    :param image_size: The image size, if creating a new view, defaults to 400
    :type image_size: int, optional
    :return: A 3D view containing the molecule
    :rtype: py3Dmol.view
    """
    if view is None:
        view = create_view(image_size=image_size)

    view.addArrow(
        {
            "start": dict(zip("xyz", orig_xyz)),
            "end": dict(zip("xyz", numpy.add(vec_xyz, orig_xyz))),
            "color": color,
        }
    )
    view.zoomTo()
    return view
