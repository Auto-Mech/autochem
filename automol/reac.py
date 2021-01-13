""" reaction classifiers and reaction-class-specific functions

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.
"""
import itertools
import more_itertools as mit
import numpy
from qcelemental import constants as qcc
import automol.convert.graph
import automol.geom.ts
from automol import par
from automol.graph import ts
from automol.graph import atom_keys
from automol.graph import bond_keys
from automol.graph import string
from automol.graph import atom_count
from automol.graph import heavy_atom_count
from automol.graph import electron_count
from automol.graph import union
from automol.graph import explicit
from automol.graph import add_bonds
from automol.graph import remove_bonds
from automol.graph import full_isomorphism
from automol.graph import union_from_sequence
from automol.graph import unsaturated_atom_keys
from automol.graph import without_stereo_parities
from automol.graph import sorted_atom_neighbor_keys
from automol.graph import add_atom_explicit_hydrogen_keys
from automol.graph import rings_bond_keys
from automol.graph import rings_atom_keys
from automol.graph import cycle_ring_atom_key_to_front

ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')


class GraphReaction:
    """ Describes a specific reaction

    Methods with trailing underscores have a side-effect

    :param class_: the name of the reaction class
    :type class_: str
    :param forward_ts_graph: a graph representing the transition state in the
        forward direction; keys must match `reactants_keys`
    :param backward_ts_graph: a graph representing the transition state in the
        backward direction; keys must match `products_keys`
    :param reactants_keys: a sequence of keys, one for each reactant, for
        extracting the reactants in order from the `forward_ts_graph`
    :type reactants_keys: tuple[tuple[int]]
    :param products_keys: a sequence of keys, one for each product, for
        extracting the products in order from the `product_ts_graph`
    :type products_keys: tuple[tuple[int]]
    """

    def __init__(self, rxn_cls, forw_tsg, back_tsg, rcts_keys, prds_keys):
        """ constructor
        """
        rcts_keys = tuple(map(tuple, map(sorted, rcts_keys)))
        prds_keys = tuple(map(tuple, map(sorted, prds_keys)))

        # Check the reaction class
        assert par.is_reaction_class(rxn_cls), (
            "{} is not a reaction class".format(rxn_cls))

        # Check the reactant keys and the forward transition state graph
        all_rcts_keys = set(itertools.chain(*rcts_keys))
        assert all_rcts_keys == atom_keys(forw_tsg), (
            "{} != {}".format(str(all_rcts_keys), str(atom_keys(forw_tsg))))

        # Check the product keys and the backward transition state graph
        all_prds_keys = set(itertools.chain(*prds_keys))
        assert all_prds_keys == atom_keys(back_tsg), (
            "{} != {}".format(str(all_prds_keys), str(atom_keys(back_tsg))))

        # Check that the reactants and products are consistent
        assert full_isomorphism(ts.reverse(forw_tsg), back_tsg)

        # Set attributes
        self.class_ = rxn_cls
        self.reactants_keys = rcts_keys
        self.products_keys = prds_keys
        self.forward_ts_graph = forw_tsg
        self.backward_ts_graph = back_tsg

    def sort_order(self):
        """ determine the appropriate sort order for reactants and products,
        based on their keys
        """
        if len(self.reactants_keys) == 1:
            rct_idxs = [0]
        else:
            rct_keys = list(map(tuple, map(sorted, self.reactants_keys)))
            rct_idxs = numpy.argsort(numpy.array(rct_keys, dtype=object))

        if len(self.products_keys) == 1:
            prd_idxs = [0]
        else:
            prd_keys = list(map(tuple, map(sorted, self.products_keys)))
            prd_idxs = numpy.argsort(numpy.array(prd_keys, dtype=object))

        rct_idxs, prd_idxs = map(tuple, (rct_idxs, prd_idxs))
        return rct_idxs, prd_idxs

    def key_map(self, reverse=False):
        """ get the key map taking atoms from the reactant into atoms from the
        product
        """
        iso_dct = full_isomorphism(ts.reverse(self.forward_ts_graph),
                                   self.backward_ts_graph)
        if reverse:
            iso_dct = dict(map(reversed, iso_dct.items()))

        return iso_dct

    def reverse_(self):
        """ get the reaction class for the reverse reaction
        """
        self.class_ = par.reverse_reaction_class(self.class_)
        self.forward_ts_graph, self.backward_ts_graph = (
            self.backward_ts_graph, self.forward_ts_graph)
        self.reactants_keys, self.products_keys = (
            self.products_keys, self.reactants_keys)

    def standardize_keys_(self):
        """ standardize keys and, optionally, sort reactant and product
        geometries in the standard order
        """
        rct_keys = list(map(sorted, self.reactants_keys))
        prd_keys = list(map(sorted, self.products_keys))
        rct_key_dct = {k: i for i, k in enumerate(itertools.chain(*rct_keys))}
        prd_key_dct = {k: i for i, k in enumerate(itertools.chain(*prd_keys))}
        self.reactants_keys = tuple(tuple(map(rct_key_dct.__getitem__, keys))
                                    for keys in self.reactants_keys)
        self.products_keys = tuple(tuple(map(prd_key_dct.__getitem__, keys))
                                   for keys in self.products_keys)
        self.forward_ts_graph = automol.graph.relabel(self.forward_ts_graph,
                                                      rct_key_dct)
        self.backward_ts_graph = automol.graph.relabel(self.backward_ts_graph,
                                                       prd_key_dct)

    def standardize_keys_and_sort_geometries_(self, rct_geos, prd_geos):
        """ standardize keys and line up geometries to match
        """
        rct_idxs, prd_idxs = self.sort_order()
        rct_geos = tuple(map(rct_geos.__getitem__, rct_idxs))
        prd_geos = tuple(map(prd_geos.__getitem__, prd_idxs))
        self.standardize_keys_()
        return rct_geos, prd_geos

    def is_standardized(self):
        """ has this GraphReaction been standardized?
        """
        other = self.copy()
        other.standardize_keys_()
        return self == other

    def copy(self):
        """ return a copy of this GraphReaction
        """
        return GraphReaction(
            self.class_, self.forward_ts_graph, self.backward_ts_graph,
            self.reactants_keys, self.products_keys)

    def __eq__(self, other):
        """ equality operator
        """
        if not isinstance(other, GraphReaction):
            ret = False
        else:
            ret = (self.class_ == other.class_ and
                   self.forward_ts_graph == other.forward_ts_graph and
                   self.backward_ts_graph == other.backward_ts_graph and
                   self.reactants_keys == other.reactants_keys and
                   self.products_keys == other.products_keys)
        return ret


def trivial(rct_gras, prd_gras):
    """ find a trivial reaction, with the same reactants and products
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == len(prd_gras):
        prd_gras = list(prd_gras)

        rct_idxs = []
        prd_idxs = []

        # One at a time, find matches for each reactant; track the positions to
        # get the right sort order
        for rct_idx, rct_gra in enumerate(rct_gras):
            prd_idx = next((idx for idx, prd_gra in enumerate(prd_gras)
                            if full_isomorphism(rct_gra, prd_gra)), None)

            if prd_idx is not None:
                rct_idxs.append(rct_idx)
                prd_idxs.append(prd_idx)
                prd_gras.pop(prd_idx)
            else:
                break

        if rct_idxs and prd_idxs:
            # reorder the reactants and products
            rct_gras = list(map(rct_gras.__getitem__, rct_idxs))
            prd_gras = list(map(prd_gras.__getitem__, prd_idxs))

            rcts_gra = union_from_sequence(rct_gras)
            prds_gra = union_from_sequence(prd_gras)

            rxns.append(GraphReaction(
                rxn_cls=par.ReactionClass.TRIVIAL,
                forw_tsg=ts.graph(rcts_gra, [], []),
                back_tsg=ts.graph(prds_gra, [], []),
                rcts_keys=list(map(atom_keys, rct_gras)),
                prds_keys=list(map(atom_keys, prd_gras)),
            ))

    return tuple(rxns)


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migrations(rct_gras, prd_gras):
    """ find hydrogen migrations consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Hydrogen migrations are identified by adding a hydrogen to an unsaturated
    site of the reactant and adding a hydrogen to an unsaturated site of the
    product and seeing if they match up. If so, we have a hydrogen migration
    between these two sites.
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 1:
        gra1, = rct_gras
        gra2, = prd_gras
        h_atm_key1 = max(atom_keys(gra1)) + 1
        h_atm_key2 = max(atom_keys(gra2)) + 1

        atm_keys1 = unsaturated_atom_keys(gra1)
        atm_keys2 = unsaturated_atom_keys(gra2)
        for atm_key1, atm_key2 in itertools.product(atm_keys1, atm_keys2):
            gra1_h = add_atom_explicit_hydrogen_keys(
                gra1, {atm_key1: [h_atm_key1]})
            gra2_h = add_atom_explicit_hydrogen_keys(
                gra2, {atm_key2: [h_atm_key2]})

            iso_dct = full_isomorphism(gra1_h, gra2_h)
            inv_dct = dict(map(reversed, iso_dct.items()))
            if inv_dct:
                f_frm_bnd_key = (atm_key1, inv_dct[h_atm_key2])
                f_brk_bnd_key = (inv_dct[atm_key2], inv_dct[h_atm_key2])
                b_frm_bnd_key = (atm_key2, iso_dct[h_atm_key1])
                b_brk_bnd_key = (iso_dct[atm_key1], iso_dct[h_atm_key1])
                forw_tsg = ts.graph(gra1,
                                    frm_bnd_keys=[f_frm_bnd_key],
                                    brk_bnd_keys=[f_brk_bnd_key])
                back_tsg = ts.graph(gra2,
                                    frm_bnd_keys=[b_frm_bnd_key],
                                    brk_bnd_keys=[b_brk_bnd_key])

                # Create the reaction object
                rxns.append(GraphReaction(
                    rxn_cls=par.ReactionClass.HYDROGEN_MIGRATION,
                    forw_tsg=forw_tsg,
                    back_tsg=back_tsg,
                    rcts_keys=[atom_keys(gra1)],
                    prds_keys=[atom_keys(gra2)],
                ))

    return tuple(rxns)


def hydrogen_migration_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for a hydrogen migration transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.HYDROGEN_MIGRATION
    assert rxn.is_standardized()
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    frm_bnd_dist = 1.7 if angstrom else 1.7 * ANG2BOHR

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key] = frm_bnd_dist

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, angstrom=angstrom)

    geo_init, = rct_geos
    relax_ang = True
    relax_tors = True

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# 2. Beta scissions
def beta_scissions(rct_gras, prd_gras):
    """ find beta scission reactions

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Implemented as the reverse of additions.
    """
    rxns = additions(prd_gras, rct_gras)
    for rxn in rxns:
        rxn.reverse_()

    return tuple(rxns)


def beta_scission_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for a beta scission transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.BETA_SCISSION
    assert rxn.is_standardized()
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    brk_bnd_dist = 1.5 if angstrom else 1.5 * ANG2BOHR

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[brk_bnd_key] = brk_bnd_dist

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, angstrom=angstrom)

    geo_init, = rct_geos
    relax_ang = False
    relax_tors = False

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# 3. Ring-forming scissions
def ring_forming_scissions(rct_gras, prd_gras):
    """ find ring-forming scissions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Ring-forming scissions are found by breaking ring-bonds on one product and
    joining the ends to unsaturated sites on the other product
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 2:
        rgra, = rct_gras
        pgra = union_from_sequence(prd_gras)
        for pgra1, pgra2 in itertools.permutations(prd_gras):
            bnd_keys = list(itertools.chain(*rings_bond_keys(pgra1)))
            atm_keys = unsaturated_atom_keys(pgra2)

            for bnd_key, atm_key in itertools.product(bnd_keys, atm_keys):
                # Break a ring bond
                gra = remove_bonds(pgra, [bnd_key])

                for end_key in bnd_key:
                    # Add to one end of the broken ring
                    fgra = add_bonds(gra, [(atm_key, end_key)])
                    inv_dct = full_isomorphism(fgra, rgra)
                    if inv_dct:
                        other_end_key, = bnd_key - {end_key}
                        f_frm_bnd_key = (inv_dct[end_key],
                                         inv_dct[other_end_key])
                        f_brk_bnd_key = (inv_dct[end_key], inv_dct[atm_key])
                        b_frm_bnd_key = (end_key, atm_key)
                        b_brk_bnd_key = (end_key, other_end_key)
                        forw_tsg = ts.graph(rgra,
                                            frm_bnd_keys=[f_frm_bnd_key],
                                            brk_bnd_keys=[f_brk_bnd_key])
                        back_tsg = ts.graph(pgra,
                                            frm_bnd_keys=[b_frm_bnd_key],
                                            brk_bnd_keys=[b_brk_bnd_key])

                        # Create the reaction object
                        rxns.append(GraphReaction(
                            rxn_cls=par.ReactionClass.RING_FORM_SCISSION,
                            forw_tsg=forw_tsg,
                            back_tsg=back_tsg,
                            rcts_keys=[atom_keys(rgra)],
                            prds_keys=[atom_keys(pgra1), atom_keys(pgra2)],
                        ))

    return tuple(rxns)


def ring_forming_scission_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for a ring-forming scission transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.RING_FORM_SCISSION
    assert rxn.is_standardized()
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    frm_bnd_dist = 1.7 if angstrom else 1.7 * ANG2BOHR
    brk_bnd_dist = 1.5 if angstrom else 1.5 * ANG2BOHR
    a234 = 85.
    d1234 = 170.

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key] = frm_bnd_dist
    dist_dct[brk_bnd_key] = brk_bnd_dist

    key2, = frm_bnd_key & brk_bnd_key
    key1, = frm_bnd_key - {key2}
    key3, = brk_bnd_key - {key2}
    ang_dct = {(key2, key3, None): a234}
    dih_dct = {(key1, key2, key3, None): d1234}

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, ang_dct=ang_dct, dih_dct=dih_dct, degree=True,
        angstrom=angstrom)

    geo_init, = rct_geos
    relax_ang = True
    relax_tors = True

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# 4. Eliminations
def eliminations(rct_gras, prd_gras):
    """ find eliminations consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Eliminations are identified by forming a bond between an attacking heavy
    atom and another atom not initially bonded to it, forming a ring. The bond
    adjacent to the attacked atom is then broken, along with a second bond in
    the ring, downstream of the attacking heavy atom, away from the attacked
    atom.
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 2:
        rgra, = rct_gras
        pgra = union_from_sequence(prd_gras)

        rngb_keys = sorted_atom_neighbor_keys(rgra)

        frm1_keys = atom_keys(rgra, excl_syms=('H',))
        frm2_keys = atom_keys(rgra)
        bnd_keys = bond_keys(rgra)

        frm_bnd_keys = [(frm1_key, frm2_key) for frm1_key, frm2_key
                        in itertools.product(frm1_keys, frm2_keys)
                        if frm1_key != frm2_key and
                        not frozenset({frm1_key, frm2_key}) in bnd_keys]

        for frm1_key, frm2_key in frm_bnd_keys:
            # Bond the radical atom to the hydrogen atom
            rgra_ = add_bonds(rgra, [(frm2_key, frm1_key)])

            # Get keys to the ring formed by this extra bond
            rng_keys = next((ks for ks in rings_atom_keys(rgra_)
                             if frm2_key in ks and frm1_key in ks), None)
            if rng_keys is not None:
                for nfrm2_key in rngb_keys[frm2_key]:
                    # Break the bond between the attacked atom and its neighbor
                    rgra_ = remove_bonds(rgra_, [(frm2_key, nfrm2_key)])

                    # Sort the ring keys so that they start with the radical
                    # atom and end with the hydrogen atom
                    keys = cycle_ring_atom_key_to_front(rng_keys, frm1_key,
                                                        end_key=frm2_key)

                    # Break one ring bond at a time, starting from the rind,
                    # and see what we get
                    for brk_key1, brk_key2 in mit.windowed(keys[:-1], 2):
                        gra = remove_bonds(rgra_, [(brk_key1, brk_key2)])

                        inv_dct = full_isomorphism(gra, pgra)
                        if inv_dct:
                            f_frm_bnd_key = (frm2_key, frm1_key)
                            f_brk_bnd_key1 = (frm2_key, nfrm2_key)
                            f_brk_bnd_key2 = (brk_key1, brk_key2)
                            b_frm_bnd_key1 = (inv_dct[frm2_key],
                                              inv_dct[nfrm2_key])
                            b_frm_bnd_key2 = (inv_dct[brk_key1],
                                              inv_dct[brk_key2])
                            b_brk_bnd_key = (inv_dct[frm2_key],
                                             inv_dct[frm1_key])
                            forw_tsg = ts.graph(rgra,
                                                frm_bnd_keys=[f_frm_bnd_key],
                                                brk_bnd_keys=[f_brk_bnd_key1,
                                                              f_brk_bnd_key2])
                            back_tsg = ts.graph(pgra,
                                                frm_bnd_keys=[b_frm_bnd_key1,
                                                              b_frm_bnd_key2],
                                                brk_bnd_keys=[b_brk_bnd_key])

                            rcts_atm_keys = list(map(atom_keys, rct_gras))
                            prds_atm_keys = list(map(atom_keys, prd_gras))

                            if inv_dct[frm2_key] not in prds_atm_keys[1]:
                                prds_atm_keys = list(reversed(prds_atm_keys))

                            # Create the reaction object
                            rxns.append(GraphReaction(
                                rxn_cls=par.ReactionClass.ELIMINATION,
                                forw_tsg=forw_tsg,
                                back_tsg=back_tsg,
                                rcts_keys=rcts_atm_keys,
                                prds_keys=prds_atm_keys,
                            ))

    return tuple(rxns)


def elimination_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for an elimination transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.ELIMINATION
    assert rxn.is_standardized()
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    frm_bnd_dist = 1.6 if angstrom else 1.6 * ANG2BOHR
    frm_rng_keys, = ts.forming_rings_atom_keys(rxn.forward_ts_graph)

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key] = frm_bnd_dist

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, rings_keys=[frm_rng_keys], degree=True,
        angstrom=angstrom)

    geo_init, = rct_geos
    relax_ang = True
    relax_tors = True

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstractions(rct_gras, prd_gras):
    """ find hydrogen abstractions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Hydrogen abstractions are identified first by checking whether the
    molecular formulas are consistent with a reaction of the form R1H + R2 =>
    R2H + R1. If they do, we identify the abstraction sites by adding hydrogens
    to unsaturated sites of the R1 product to see if we get the R1H reactant.
    We then do the same for the R2 reactant and the R2H product.
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_fmls = list(map(automol.convert.graph.formula, rct_gras))
        prd_fmls = list(map(automol.convert.graph.formula, prd_gras))

        ret = automol.formula.reac.argsort_hydrogen_abstraction(
            rct_fmls, prd_fmls)
        if ret:
            rct_idxs_, prd_idxs_ = ret
            rct_gras = list(map(rct_gras.__getitem__, rct_idxs_))
            prd_gras = list(map(prd_gras.__getitem__, prd_idxs_))

            q1h_gra, q2_gra = rct_gras
            q2h_gra, q1_gra = prd_gras

            rets1 = _partial_hydrogen_abstraction(q1h_gra, q1_gra)
            rets2 = _partial_hydrogen_abstraction(q2h_gra, q2_gra)
            for ret1, ret2 in itertools.product(rets1, rets2):
                f_q1h_q_atm_key, f_q1h_h_atm_key, b_q2_q_atm_key = ret1
                b_q1h_q_atm_key, b_q1h_h_atm_key, f_q2_q_atm_key = ret2

                # Create the forward/backward ts graphs
                rcts_gra = union_from_sequence(rct_gras)
                prds_gra = union_from_sequence(prd_gras)
                f_frm_bnd_key = (f_q2_q_atm_key, f_q1h_h_atm_key)
                f_brk_bnd_key = (f_q1h_q_atm_key, f_q1h_h_atm_key)
                b_frm_bnd_key = (b_q2_q_atm_key, b_q1h_h_atm_key)
                b_brk_bnd_key = (b_q1h_q_atm_key, b_q1h_h_atm_key)
                forw_tsg = ts.graph(rcts_gra,
                                    frm_bnd_keys=[f_frm_bnd_key],
                                    brk_bnd_keys=[f_brk_bnd_key])
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[b_frm_bnd_key],
                                    brk_bnd_keys=[b_brk_bnd_key])

                # Create the reaction object
                rxns.append(GraphReaction(
                    rxn_cls=par.ReactionClass.HYDROGEN_ABSTRACTION,
                    forw_tsg=forw_tsg,
                    back_tsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                ))

    return tuple(rxns)


def hydrogen_abstraction_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for a hydrogen abstraction transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.HYDROGEN_ABSTRACTION
    assert rxn.is_standardized()
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    frm_bnd_dist = 1.6 if angstrom else 1.6 * ANG2BOHR
    a123 = 170.
    a234 = 85.

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key] = frm_bnd_dist

    # key2 the hydrogen atom and key3 is the attacking atom; this is guaranteed
    # to be true since the GraphReaction has been standardized
    key2, key3 = sorted(frm_bnd_key)
    ang_dct = {(None, key2, key3): a123, (key2, key3, None): a234}

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, ang_dct=ang_dct, degree=True, angstrom=angstrom)

    geo1, geo2 = rct_geos
    geo_init = automol.geom.ts.join(geo1, geo2, key2=key2, key3=key3,
                                    r23=frm_bnd_dist, a123=a123, a234=a234)
    relax_ang = False
    relax_tors = False

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# 2. Additions
def additions(rct_gras, prd_gras):
    """ find additions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Additions are identified by joining an unsaturated site on one reactant to
    an unsaturated site on the other. If the result matches the products, this
    is an addition reaction.
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        x_gra, y_gra = rct_gras
        prd_gra, = prd_gras
        x_atm_keys = unsaturated_atom_keys(x_gra)
        y_atm_keys = unsaturated_atom_keys(y_gra)

        for x_atm_key, y_atm_key in itertools.product(x_atm_keys, y_atm_keys):
            xy_gra = add_bonds(
                union(x_gra, y_gra), [{x_atm_key, y_atm_key}])

            iso_dct = full_isomorphism(xy_gra, prd_gra)
            if iso_dct:
                rcts_gra = union_from_sequence(rct_gras)
                prds_gra = prd_gra
                f_frm_bnd_key = (x_atm_key, y_atm_key)
                b_brk_bnd_key = (iso_dct[x_atm_key], iso_dct[y_atm_key])
                forw_tsg = ts.graph(rcts_gra,
                                    frm_bnd_keys=[f_frm_bnd_key],
                                    brk_bnd_keys=[])
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[],
                                    brk_bnd_keys=[b_brk_bnd_key])

                # sort the reactants so that the largest species is first
                rct_idxs = _argsort_reactants(rct_gras)
                rct_gras = list(map(rct_gras.__getitem__, rct_idxs))

                # Create the reaction object
                rxns.append(GraphReaction(
                    rxn_cls=par.ReactionClass.ADDITION,
                    forw_tsg=forw_tsg,
                    back_tsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                ))

    return tuple(rxns)


def addition_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for an addition transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.ADDITION
    assert rxn.is_standardized()
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    frm_bnd_dist = 1.9 if angstrom else 1.9 * ANG2BOHR
    a123 = 85.
    a234 = 85.
    d1234 = 85.

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key] = frm_bnd_dist

    # key2 the hydrogen atom and key3 is the attacking atom; this is guaranteed
    # to be true since the GraphReaction has been standardized
    key2, key3 = sorted(frm_bnd_key)
    ang_dct = {(None, key2, key3): a123, (key2, key3, None): a234}
    dih_dct = {(None, key2, key3, None): d1234}

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, ang_dct=ang_dct, dih_dct=dih_dct, degree=True,
        angstrom=angstrom)

    geo1, geo2 = rct_geos
    geo_init = automol.geom.ts.join(geo1, geo2, key2=key2, key3=key3,
                                    r23=frm_bnd_dist, a123=a123, a234=a234,
                                    d1234=d1234)
    relax_ang = False
    relax_tors = False

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# 3. Insertions
def insertions(rct_gras, prd_gras):
    """ find insertions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Implemented as the reverse of an addition reaction.
    """
    rxns = eliminations(prd_gras, rct_gras)
    for rxn in rxns:
        rxn.reverse_()

    return tuple(rxns)


def insertion_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for an insertion transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.INSERTION
    assert rxn.is_standardized()

    # set the formed bond distance based on the number of hydrogens
    frm_bnd_dist_dct = {
        0: 1.9,
        1: 1.6,
        2: 1.2,
    }

    tsg = rxn.forward_ts_graph

    frm_bnd_key1, frm_bnd_key2 = ts.forming_bond_keys(tsg)
    hcnt1 = automol.graph.atom_count_by_type(tsg, 'H', keys=frm_bnd_key1)
    hcnt2 = automol.graph.atom_count_by_type(tsg, 'H', keys=frm_bnd_key2)
    frm_bnd_dist1 = frm_bnd_dist_dct[hcnt1]
    frm_bnd_dist2 = frm_bnd_dist_dct[hcnt2]
    a123 = 170.
    frm_rng_keys, = ts.forming_rings_atom_keys(tsg)

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key1] = frm_bnd_dist1
    dist_dct[frm_bnd_key2] = frm_bnd_dist2

    gra = ts.reactants_graph(tsg)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, rings_keys=[frm_rng_keys], angstrom=angstrom)

    # Join on the end without hydrogen, if there is one
    if hcnt1 <= hcnt2:
        key2, key3 = sorted(frm_bnd_key1)
    else:
        key2, key3 = sorted(frm_bnd_key2)

    geo1, geo2 = rct_geos
    geo_init = automol.geom.ts.join(geo1, geo2, key2=key2, key3=key3,
                                    r23=frm_bnd_dist1, a123=a123)
    relax_ang = True
    relax_tors = True

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


# 4. Substitutions
def substitutions(rct_gras, prd_gras):
    """ find substitutions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Substitutions are identified by breaking one bond in the reactants and one
    bond from the products and checking for isomorphism.
    """
    _assert_is_valid_reagent_graph_list(rct_gras)
    _assert_is_valid_reagent_graph_list(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_gra = union_from_sequence(rct_gras)
        prd_gra = union_from_sequence(prd_gras)

        for rgra1, rgra2 in itertools.permutations(rct_gras):
            bnd_keys = bond_keys(rgra1)
            rad_keys = unsaturated_atom_keys(rgra2)

            for bnd_key, rad_key in itertools.product(bnd_keys, rad_keys):
                gra = remove_bonds(rct_gra, [bnd_key])

                for brk_key1 in bnd_key:
                    gra = add_bonds(gra, [(brk_key1, rad_key)])

                    inv_dct = full_isomorphism(gra, prd_gra)
                    if inv_dct:
                        brk_key2, = bnd_key - {brk_key1}
                        f_frm_bnd_key = (brk_key1, rad_key)
                        f_brk_bnd_key = (brk_key1, brk_key2)
                        b_frm_bnd_key = (inv_dct[brk_key1], inv_dct[brk_key2])
                        b_brk_bnd_key = (inv_dct[brk_key1], inv_dct[rad_key])

                        forw_tsg = ts.graph(rct_gra,
                                            frm_bnd_keys=[f_frm_bnd_key],
                                            brk_bnd_keys=[f_brk_bnd_key])
                        back_tsg = ts.graph(prd_gra,
                                            frm_bnd_keys=[b_frm_bnd_key],
                                            brk_bnd_keys=[b_brk_bnd_key])

                        rcts_atm_keys = [atom_keys(rgra1), atom_keys(rgra2)]

                        prds_atm_keys = list(map(atom_keys, prd_gras))
                        if inv_dct[rad_key] not in prds_atm_keys[0]:
                            prds_atm_keys = list(reversed(prds_atm_keys))

                        # Create the reaction object
                        rxns.append(GraphReaction(
                            rxn_cls=par.ReactionClass.SUBSTITUTION,
                            forw_tsg=forw_tsg,
                            back_tsg=back_tsg,
                            rcts_keys=rcts_atm_keys,
                            prds_keys=prds_atm_keys,
                        ))

    return tuple(rxns)


def substitution_ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ atom distance ranges for a substitution transition state

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    """
    assert rxn.class_ == par.ReactionClass.SUBSTITUTION
    assert rxn.is_standardized()

    # set the formed bond distance based on the number of hydrogens
    frm_bnd_dist_dct = {
        0: 1.9,
        1: 1.6,
        2: 1.2,
    }

    tsg = rxn.forward_ts_graph
    a123 = 170.
    a234 = 95.

    frm_bnd_key, = ts.forming_bond_keys(tsg)
    hcnt = automol.graph.atom_count_by_type(tsg, 'H', keys=frm_bnd_key)
    frm_bnd_dist = frm_bnd_dist_dct[hcnt]

    dist_dct = automol.geom.ts.distances(rct_geos, angstrom=angstrom)
    dist_dct[frm_bnd_key] = frm_bnd_dist

    # key2 the hydrogen atom and key3 is the attacking atom; this is guaranteed
    # to be true since the GraphReaction has been standardized
    key2, key3 = sorted(frm_bnd_key)
    ang_dct = {(None, key2, key3): a123, (key2, key3, None): a234}

    gra = ts.reactants_graph(tsg)
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        gra, dist_dct, ang_dct=ang_dct, degree=True, angstrom=angstrom)

    geo1, geo2 = rct_geos
    geo_init = automol.geom.ts.join(geo1, geo2, key2=key2, key3=key3,
                                    r23=frm_bnd_dist, a123=a123, a234=a234)
    relax_ang = True
    relax_tors = True

    keys = sorted(atom_keys(gra))
    xmat = automol.geom.coordinates(geo_init, angstrom=angstrom)
    lmat, umat = automol.graph.embed.join_distance_bounds_matrices(
        gra, keys, dist_range_dct, geos=rct_geos, relax_angles=relax_ang,
        relax_torsions=relax_tors)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    return xmat, lmat, umat, chi_dct, pla_dct


def classify(rct_gras, prd_gras):
    """ classify a reaction
    """
    # check whether this is a valid reaction
    rct_fmls = list(map(automol.convert.graph.formula, rct_gras))
    prd_fmls = list(map(automol.convert.graph.formula, prd_gras))
    rct_strs = list(map(automol.formula.string, rct_fmls))
    prd_strs = list(map(automol.formula.string, prd_fmls))
    assert automol.formula.reac.is_valid_reaction(rct_fmls, prd_fmls), (
        "Invalid reaction: {:s} -> {:s}".format(str(rct_strs), str(prd_strs)))

    # Cycle through the different finders and gather all possible reactions
    finders_ = [
        trivial,
        # unimolecular reactions
        hydrogen_migrations,
        beta_scissions,
        ring_forming_scissions,
        eliminations,
        # bimolecular reactions
        hydrogen_abstractions,
        additions,
        insertions,
        substitutions,
    ]

    rxns = tuple(itertools.chain(*(f_(rct_gras, prd_gras) for f_ in finders_)))

    return rxns


def ts_embedding_info(rxn, rct_geos, angstrom=True):
    """ reaction-class-specific embedding info

    :param rxn: a hydrogen migration GraphReaction object
    :param rct_geos: the reactant geometries
    :returns: an initial geometry, a dictionary of distance ranges, and two
        booleans specifying whether or not to relax angles and torsions,
        respectively
    """
    function_dct = {
        # unimolecular
        par.ReactionClass.HYDROGEN_MIGRATION:
        hydrogen_migration_ts_embedding_info,
        par.ReactionClass.BETA_SCISSION:
        beta_scission_ts_embedding_info,
        par.ReactionClass.RING_FORM_SCISSION:
        ring_forming_scission_ts_embedding_info,
        par.ReactionClass.ELIMINATION:
        elimination_ts_embedding_info,
        # bimolecular
        par.ReactionClass.HYDROGEN_ABSTRACTION:
        hydrogen_abstraction_ts_embedding_info,
        par.ReactionClass.ADDITION:
        addition_ts_embedding_info,
        par.ReactionClass.INSERTION:
        insertion_ts_embedding_info,
        par.ReactionClass.SUBSTITUTION:
        substitution_ts_embedding_info,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, rct_geos, angstrom=angstrom)
    return ret


# helpers
def _partial_hydrogen_abstraction(qh_gra, q_gra):
    rets = []

    h_atm_key = max(atom_keys(q_gra)) + 1
    uns_atm_keys = unsaturated_atom_keys(q_gra)
    for atm_key in uns_atm_keys:
        q_gra_h = add_atom_explicit_hydrogen_keys(
            q_gra, {atm_key: [h_atm_key]})
        inv_atm_key_dct = full_isomorphism(q_gra_h, qh_gra)
        if inv_atm_key_dct:
            qh_q_atm_key = inv_atm_key_dct[atm_key]
            qh_h_atm_key = inv_atm_key_dct[h_atm_key]
            q_q_atm_key = atm_key
            rets.append((qh_q_atm_key, qh_h_atm_key, q_q_atm_key))

    return rets


def _assert_is_valid_reagent_graph_list(gras):
    gras_str = '\n---\n'.join(map(string, gras))
    assert _are_all_explicit(gras), (
        "Implicit hydrogens are not allowed here!\nGraphs:\n{}"
        .format(gras_str))
    assert _have_no_stereo_assignments(gras), (
        "Stereo assignments are not allowed here!\nGraphs:\n{}"
        .format(gras_str))
    assert _have_no_common_atom_keys(gras), (
        "Overlapping atom keys are not allowed here!\nGraphs:\n{}"
        .format(gras_str))


def _are_all_explicit(gras):
    return all(gra == explicit(gra) for gra in gras)


def _have_no_stereo_assignments(gras):
    return all(gra == without_stereo_parities(gra) for gra in gras)


def _have_no_common_atom_keys(gras):
    atm_keys = list(itertools.chain(*map(atom_keys, gras)))
    return len(atm_keys) == len(set(atm_keys))


def _argsort_reactants(gras):

    def __sort_value(args):
        _, gra = args
        val = (-heavy_atom_count(gra),
               -atom_count(gra),
               -electron_count(gra))
        return val

    idxs = tuple(idx for idx, gra in sorted(enumerate(gras), key=__sort_value))
    return idxs
