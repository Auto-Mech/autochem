""" transition state graph data structure
"""
import operator
import functools
import yaml
import automol.par
from automol.util import dict_
from automol.graph._graph_base import string
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import bond_keys
from automol.graph._graph_base import bond_orders
from automol.graph._graph_base import set_bond_orders
from automol.graph._graph_base import yaml_dictionary
from automol.graph._graph_base import from_yaml_dictionary
from automol.graph._graph import subgraph
from automol.graph._graph import add_bonds
from automol.graph._graph import remove_bonds
from automol.graph._graph import union_from_sequence
from automol.graph._graph import connected_components


class TSGraph():
    """ a transition-state graph, along with its associated information

    :param graph: a graph representing the transition state; forming bonds are
        included as 0-order bonds; breaking bonds are regular bonds
    :param reaction_class: the name of this graph's reaction class
    :param form_bonds: bonds that are forming in this transition state
    :param break_bonds: bonds that are breaking in this transition state
    :param reactant_keys: a tuple of tuples containing the atom keys for each
        reactant, in order
    """

    # constructors
    @classmethod
    def from_data(cls, rct_gras, frm_bnd_keys, brk_bnd_keys, rxn_class):
        """ constructor

        :param rct_gras: graphs for the reactant; in the desired order
            (attacking atom second for bimolecular reactions) with
            non-overlapping keys
        :param frm_bnd_keys: bond keys formed in the reaction
        :param brk_bnd_keys: bond keys broken in the reaction
        :param rxn_class: the reaction class name, as a string
        """
        frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
        brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))

        rcts_keys = tuple(map(atom_keys, rct_gras))

        # Now, form the graph itself, checking that the breaking bonds are
        # present and adding in the forming ones as 0-order bonds
        gra = union_from_sequence(rct_gras)

        gra = add_bonds(gra, frm_bnd_keys, {k: 0 for k in frm_bnd_keys})
        return cls(graph=gra, reaction_class=rxn_class,
                   form_bonds=frm_bnd_keys, break_bonds=brk_bnd_keys,
                   reactants_keys=rcts_keys)

    # operations
    def reactant_graphs(self):
        """ reactant graphs for this transition state
        """
        gra = remove_bonds(self.graph, self.form_bonds)
        rct_gras = tuple(subgraph(gra, keys) for keys in self.reactants_keys)
        return rct_gras

    def product_graphs(self):
        """ product graphs for this transition state
        """
        return self.reverse().reactant_graphs()

    def reverse(self):
        """ transition-state graph for the reverse reaction

        (this will need to be case-by-case for each reaction class)
        """
        rct_gras = self.product_graphs()
        frm_bnd_keys = self.break_bonds
        brk_bnd_keys = self.form_bonds
        rxn_class = automol.par.reverse_reaction_class(self.reaction_class)
        if rxn_class is None:
            raise NotImplementedError("Reversal for {} is not implemented."
                                      .format(self.reaction_class))
        return TSGraph.from_data(
            rct_gras=rct_gras, frm_bnd_keys=frm_bnd_keys,
            brk_bnd_keys=brk_bnd_keys, rxn_class=rxn_class)

    # I/O
    def string(self, one_indexed=True):
        """ write the object to a string
        """
        def _encode_bond(bnd_key):
            atm1_key, atm2_key = bnd_key
            if one_indexed:
                bnd_str = '{}-{}'.format(atm1_key+1, atm2_key+1)
            else:
                bnd_str = '{}-{}'.format(atm1_key, atm2_key)
            return bnd_str

        # reformat data for YAML
        frm_bnd_strs = list(map(_encode_bond, self.form_bonds))
        brk_bnd_strs = list(map(_encode_bond, self.break_bonds))
        rcts_keys = list(map(list, self.reactants_keys))
        if one_indexed:
            rcts_keys = [[k+1 for k in ks] for ks in rcts_keys]
        yaml_gra_dct = yaml_dictionary(self.graph, one_indexed=one_indexed)

        yaml_dct = {}
        yaml_dct['reaction class'] = self.reaction_class
        yaml_dct['forming bonds'] = frm_bnd_strs
        yaml_dct['breaking bonds'] = brk_bnd_strs
        yaml_dct['reactants keys'] = rcts_keys
        yaml_dct.update(**yaml_gra_dct)
        tsg_str = yaml.dump(
            yaml_dct, default_flow_style=None, sort_keys=False)
        return tsg_str

    @classmethod
    def from_string(cls, tsg_str, one_indexed=True):
        """ construct a TSGraph object from a string
        """
        def _decode_bond(bnd_str):
            atm1_key, atm2_key = map(int, bnd_str.split('-'))
            if one_indexed:
                bnd_key = frozenset({atm1_key-1, atm2_key-1})
            else:
                bnd_key = frozenset({atm1_key, atm2_key})
            return bnd_key

        yaml_dct = yaml.load(tsg_str, Loader=yaml.FullLoader)
        rxn_class = yaml_dct['reaction class']
        frm_bnd_keys = list(map(_decode_bond, yaml_dct['forming bonds']))
        brk_bnd_keys = list(map(_decode_bond, yaml_dct['breaking bonds']))
        rcts_keys = yaml_dct['reactants keys']
        if one_indexed:
            rcts_keys = [[k-1 for k in ks] for ks in rcts_keys]
        gra = from_yaml_dictionary(yaml_dct, one_indexed=one_indexed)
        return cls(graph=gra, reaction_class=rxn_class,
                   form_bonds=frm_bnd_keys, break_bonds=brk_bnd_keys,
                   reactants_keys=rcts_keys)

    def __init__(self, graph, reaction_class, form_bonds, break_bonds,
                 reactants_keys):
        """ raw constructor
        """
        form_bonds = frozenset(map(frozenset, form_bonds))
        break_bonds = frozenset(map(frozenset, break_bonds))
        reactants_keys = tuple(map(frozenset, reactants_keys))

        # check that this is a valid reaction class
        assert automol.par.is_reaction_class(reaction_class)

        # check the forming/breaking bond keys
        bnd_keys = bond_keys(graph)
        assert form_bonds | break_bonds <= bnd_keys, (
            "Forming {} or breaking bonds {} don't match graph bonds {}"
            .format(*map(str, (form_bonds, break_bonds, bnd_keys))))

        # check the forming/breaking bond orders
        frm_ords = dict_.values_by_key(bond_orders(graph), form_bonds)
        assert all(o == 0 for o in frm_ords), (
            "Bond orders for forming bonds {} should be 0"
            .format(str(form_bonds)))

        brk_ords = dict_.values_by_key(bond_orders(graph), break_bonds)
        assert all(o == 1 for o in brk_ords), (
            "Bond orders for breaking bonds {} should be 1"
            .format(str(break_bonds)))

        # check that the reactant keys are valid
        all_rct_keys = functools.reduce(operator.or_, reactants_keys)
        atm_keys = atom_keys(graph)
        assert all_rct_keys == atm_keys, (
            "Reactant keys {} don't match graph keys:\n{}"
            .format(str(all_rct_keys), string(graph, one_indexed=False)))

        # now, set the attributes
        self.graph = graph
        self.reaction_class = reaction_class
        self.form_bonds = form_bonds
        self.break_bonds = break_bonds
        self.reactants_keys = reactants_keys
