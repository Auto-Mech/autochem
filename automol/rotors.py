import dataclasses
from typing import List, Tuple

from automol import graph, zmat


@dataclasses.dataclass
class Torsion:
    """Encodes information for a single torsion, which is one component of a rotor

    :param name: The z-matrix coordinate name
    :type name: str
    :param axis: The pair of atom keys defining the rotational axis
    :type axis: Tuple[int, int]
    :param groups: The sets of atoms keys defining the rotational groups
    :type groups: Tuple[List[int], List[int]]
    :type symmetry: The rotational symmetry number of the torsion
    :type symmetry: int
    """

    name: str
    axis: Tuple[int, int]
    groups: Tuple[List[int], List[int]]
    symmetry: int


Rotor = List[Torsion]


@dataclasses.dataclass
class Rotors:
    """Encodes information for all rotors of an equilibrium or transition state
    structure

    :param zmatrix: The z-matrix
    :type zmatrix: automol zmat data structure
    :param rotors: The list of rotors for this structure
    :type rotors: List[Rotor]
    """

    zmatrix: tuple  # automol zmat data structure
    rotors: List[Rotor]


# Constructors
def from_zmatrix(zma, gra=None, multi: bool = False) -> Rotors:
    """Construct rotors for a z-matrix

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param gra: Optionally, specify the z-matrix connectivity with a graph
    :type gra: automol graph data structure
    :param multi: Construct multi-dimensional rotors? defaults to False
    :type multi: bool, optional
    :returns: The rotors for this structure
    :rtype: Rotors
    """
    gra = zmat.graph(zma, stereo=True, dummy=True) if gra is None else gra
    lin_keys = graph.linear_atom_keys(gra, dummy=True)

    rot_bkeys = graph.rotational_bond_keys(gra, lin_keys=lin_keys)

    # Read in the torsion coordinate names and sort them in z-matrix order
    tors_names = [zmat.torsion_coordinate_name(zma, *bk) for bk in rot_bkeys]
    tors_names = [n for n in zmat.dihedral_angle_names(zma) if n in tors_names]

    tors_lst = []
    for name in tors_names:
        axis = zmat.dihedral_axis(zma, name)
        tors = Torsion(
            name=name,
            axis=axis,
            groups=graph.rotational_groups(gra, *axis),
            symmetry=graph.rotational_symmetry_number(gra, *axis, lin_keys=lin_keys),
        )

        tors_lst.append(tors)

    if multi:
        rotor_lst = group_torsions_into_multidimensional_rotors(tors_lst)
    else:
        rotor_lst = [[t] for t in tors_lst]

    return Rotors(zmatrix=zma, rotors=rotor_lst)


def group_torsions_into_multidimensional_rotors(
    tors_lst: List[Torsion], gra=None
) -> List[Rotor]:
    """Group torsions into multi-dimensional rotors, reducing dimensionality where needed

    :param tors_lst: A list of torsions
    :type tors_lst: List[Torsion]
    :param gra: A graph describing connectivity
    :type gra: automol graph data structure
    :returns: The multi-dimensional rotors
    :rtype: List[Rotor]
    """
    # If there are <= 4 torsions, group them into a single multi-dimensional rotor
    if len(tors_lst) <= 4:
        return tors_lst

    # Otherwise, reduce the dimensionality by dropping methyl rotors
    tors_lst = [t for t in tors_lst if not graph.is_methyl_rotor(gra, *t.axis)]

    # Now, if there are <= 4 torsions, group them into a single multi-dimensional rotor
    if len(tors_lst) <= 4:
        return tors_lst

    # If removing methyl rotors didn't reduce the dimensionality down to <= 4, flatten
    # this down to a list of one-dimensional rotors
    return [[t] for t in tors_lst]
