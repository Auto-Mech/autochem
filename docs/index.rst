AutoChem
=========

*A Package of the AutoMech Suite*

Andreas V. Copan, Kevin B. Moore III, Sarah N. Elliott, and Stephen J. Klippenstein

--------------------------------------------------------------------------------------

.. toctree::
    :glob:
    :maxdepth: 3

    index
--------------------------------------------------------------------------------------

Overview
~~~~~~~~

Autochem is a library of molecular toolkit packages, including automol, phydat, and transformations.

AutoMol
^^^^^^^

Automol provides the framework for any molecule object used across the AutoMech suite and
facilitates the conversion between objects.  Some of these objects are internal representations of
well known objects in chemistry. Automol, however, also has its own molecular objects like our
molecuar graph, which stores connectivity and stereo information. Perhaps the most novel features in
automol are its transformation functions on any of these objects.  Such tools enable us to search reaction profiles,
set up scans along torsional angles, and identify configurations for alternative conformers. Our tutorials
and code documentation demonstrate the extensive scope of functions within this toolkit.

|

.. centered::
    **SMILES string: 2D representation**

Simplified molecular-input line-entry system (SMILES) are a standard convention in large scale chemistry.
SMILES strings encode molecular connectivity.

.. list-table::
   :widths: 10
   :header-rows: 1

   * - Vinyl Radical (CH2CH)
   * - C=[CH]

|

|

.. centered::
    **InChI string: 2D representation**

IUPAC International Chemical Identifiers are also widely used in large scale chemistry.  The string encodes
atoms and their bond connectivity, tautomeric information, isotope information, stereochemistry, and electronic charge information in
various, and optional, layers.

.. list-table::
   :widths: 10
   :header-rows: 1

   * - Vinyl Radical (CH2CH)
   * - InChI=1S/C2H3/c1-2/h1H,2H2

|

|

.. centered::
    **Graph object: 2D representation**

The automol graph, designed for the AutoMech suite, gives atoms, their hydrogen valences, and their stereo parities as vertices with
the atom connectivies as edges.  The graph provides a 2D representation that is transformable and manipulable.  For this reason,
the graph object is essential to reaction identification.
There are two implementations of the graph object, with trivial conversion between the two, with implicit and explicit hydrogen atoms

.. table::
    :class: tight-table

    +-----------------------------------------------------------+
    | Vinyl Radical (CH2CH) Internal                            |
    +===========================================================+
    | *Implicit*                                                |
    |                                                           |
    |      ({                                                   |
    |      0: ('C', 1, None),                                   |
    |      1: ('C', 2, None)},                                  |
    |      {frozenset({0, 1}): (1, None)})                      |
    |                                                           |
    |                                                           |
    | *Explicit*                                                |
    |                                                           |
    |  ({                                                       |
    |  0: ('C', 0, None),                                       |
    |  1: ('C', 0, None),                                       |
    |  2: ('H', 0, None),                                       |
    |  3: ('H', 0, None),                                       |
    |  4: ('H', 0, None)},                                      |
    |  {                                                        |
    |  frozenset({0, 1}): (1, None),                            |
    |  frozenset({0, 2}): (1, None),                            |
    |  frozenset({1, 3}): (1, None),                            |
    |  frozenset({1, 4}): (1, None)})                           |
    +-----------------------------------------------------------+

+-------------------------------------------------------------------------+
| Vinyl Radical (CH2CH)  External                                         |
+=========================================================================+
| *Implicit*                                                              |
|                                                                         |
| atoms:                                                                  |
|                                                                         |
|   1: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}     |
|                                                                         |
|   2: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}     |
|                                                                         |
| bonds:                                                                  |
|                                                                         |
|    1-2: {order: 1, stereo_parity: null}                                 |
|                                                                         |
|                                                                         |
| *Explicit*                                                              |
|                                                                         |
| atoms:                                                                  |
|                                                                         |
|   1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}     |
|                                                                         |
|   2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}     |
|                                                                         |
|   3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}     |
|                                                                         |
|   4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}     |
|                                                                         |
|   5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}     |
|                                                                         |
| bonds:                                                                  |
|                                                                         |
|    1-2: {order: 1, stereo_parity: null}                                 |
|                                                                         |
|    1-3: {order: 1, stereo_parity: null}                                 |
|                                                                         |
|    2-4: {order: 1, stereo_parity: null}                                 |
|                                                                         |
|    2-5: {order: 1, stereo_parity: null}                                 |
+-------------------------------------------------------------------------+

|

|

.. centered::
    **Geometry Object: 3D representation**

Our most simple 3D representation of the molecule gives the location
of its atoms in Cartesian coordinates. Externally, coordinates are in
Angstrom and, internally, coordinates are in Bohr.  Automol has simple
conversions to/from the 2D descriptors from/to a geometry object using
calls to RDkit and OpenBabel.

.. table::
    :class: tight-table

    +---------------------------------------------------------------------------+
    | Vinyl Radical (CH2CH)  Internal                                           |
    +===========================================================================+
    | (('C', (-1.7385058975595666, 0.2715876398092092, 0.7998932295794201)),    |
    | ('C', (0.9772590604406942, 0.0027272925514021072, 0.028477437733657032)), |
    | ('H', (-3.250682519670039, -0.08885591643834877, -0.5801409688754181)),   |
    | ('H', (2.160193047410571, 1.6990073384324165, -0.20556820450994961)),     |
    | ('H', (1.8517363093783363, -1.8844663543546802, -0.04266149392770653)))   |
    +---------------------------------------------------------------------------+

.. table::
    :class: tight-table

    +----+------------+--------------+--------------------------------+
    | Vinyl Radical (CH2CH)  External                                 |
    +====+============+==============+================================+
    |  C |  -0.919978 |     0.143718 |   0.423285                     |
    +----+------------+--------------+--------------------------------+
    |  C |  0.517143  |    0.001443  |  0.015070                      |
    +----+------------+--------------+--------------------------------+
    |  H | -1.720187  |   -0.047021  | -0.306997                      |
    +----+------------+--------------+--------------------------------+
    |  H | 1.143125   |  0.899076    |  -0.108782                     |
    +----+------------+--------------+--------------------------------+
    |  H | 0.979897   |  -0.997217   |  -0.022575                     |
    +----+------------+--------------+--------------------------------+

|

|

.. centered::
    **Z-matrix object: 3D representation**

The Z-matrix is a standard representation for a molecule in computational
chemistry.  It uses internal coordinates between atoms, specifying
an atom *a*, its distance from atom *b*, the central angle *a*-*b*-*c*
for a third atom c, and the dihedral angle *a*-*b*-*c*-*d* for a fourth
atom *d*. Many possible Z-matices can represent a single molecule, but
automol will choose a set of internal coordinates that captures the connections
and torisonal angles along the backbone of the molecule. Automol
provides functions to seamlessly convert between geometry and Z-matrix
object and the 2D representations.

.. table::
    :class: tight-table

    +------------------------------------------------------------------------------+
    | Vinyl Radical (CH2CH)  Internal                                              |
    +==============================================================================+
    | (('C', (None, None, None), (None, None, None), (None, None, None)),          |
    | ('C', (0, None, None), ('R1', None, None), (2.4533342041496025, None, None)),|
    | ('H',                                                                        |
    | (0, 1, None),                                                                |
    | ('R2', 'A2', None),                                                          |
    | (2.04887881761982, 2.106806084047999, None)),                                |
    | ('H',                                                                        |
    | (0, 1, 2),                                                                   |
    | ('R3', 'A3', 'D3'),                                                          |
    | (2.048878892080454, 2.106806109048073, 3.1415931294947175)),                 |
    | ('X',                                                                        |
    | (1, 0, 2),                                                                   |
    | ('R4', 'A4', 'D4'),                                                          |
    | (1.8897261254578281, 1.5707963267948968, 1.4901161193847656e-08)),           |
    | ('H',                                                                        |
    | (1, 4, 0),                                                                   |
    | (('R5', 'A5', 'D5'),                                                         |
    | ((2.0154889471042523, 1.5707962916270888, 3.141592595877844)))               |
    +------------------------------------------------------------------------------+

.. table::
    :class: tight-table
    :widths: 10 10 10 10 10 10 10

    +-------------------+
    | Vinyl Radical     |
    | (CH2CH)           |
    | Internal          |
    +===================+
    |C                  |
    +-+-+---------------+
    |C|1|R1             |
    +-+-+--+-+----------+
    |H|1|R2|2|A2        |
    +-+-+--+-+--+-+-----+
    |H|1|R3|2|A3|3|D3   |
    +-+-+--+-+--+-+-----+
    |X|2|R4|1|A4|3|D4   |
    +-+-+--+-+--+-+-----+
    |H|2|R5|5|A5|1|D5   |
    +-+-+--+-+--+-+-----+
    |                   |
    |                   |
    +--+-+--------------+
    |R1|=|  1.298249    |
    +--+-+--------------+
    |R1|=|  1.298249    |
    +--+-+--------------+
    |R2|=|  1.084220    |
    +--+-+--------------+
    |R3|=|  1.084220    |
    +--+-+--------------+
    |R4|=|  1.000000    |
    +--+-+--------------+
    |R5|=|  1.066551    |
    +--+-+--------------+
    |A2|=|120.711097    |
    +--+-+--------------+
    |A3|=|120.711098    |
    +--+-+--------------+
    |A4|=| 90.000000    |
    +--+-+--------------+
    |A5|=| 89.999998    |
    +--+-+--------------+
    |D3|=|180.000027    |
    +--+-+--------------+
    |D4|=|  0.000001    |
    +--+-+--------------+
    |D5|=|179.999997    |
    +--+-+--------------+

|

|

.. centered::
    **Reaction object: 2D representation**

The reaction object in automol is made up of the reaction class string,
(*e.g.*, hydrogen abstraction, beta-scission, hydrogen migration, addition, insertion,
substitution), two graphs -- one describing  the forward transformation from reactant to product
molecules and one describing the backward transformation from product to reactant molecules -- and
lists of the keys that identify either the atom vertices from reactant or product molecules.
Forming bonds, in these graphs, are given a value of 0.1, while breaking bonds are given
a value of 0.9.  The transition state itself can be represented as a Z-matrix or Geometry.

.. table::
    :class: tight-table

    +----------------------------------------------------------------------+
    | CH2CH + CH4 = CH2CH2 + CH3                                           |
    +======================================================================+
    | reaction class: hydrogen abstraction                                 |
    +----------------------------------------------------------------------+
    | forward TS atoms:                                                    |
    +----------------------------------------------------------------------+
    |   1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   2: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   6: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   8: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null} |
    +----------------------------------------------------------------------+
    |   11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null} |
    +----------------------------------------------------------------------+
    |   12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null} |
    +----------------------------------------------------------------------+
    | forward TS bonds:                                                    |
    +----------------------------------------------------------------------+
    |   1-2: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   1-3: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   1-4: {order: 0.9, stereo_parity: null}                             |
    +----------------------------------------------------------------------+
    |   1-5: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   4-6: {order: 0, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   4-7: {order: 0.1, stereo_parity: null}                             |
    +----------------------------------------------------------------------+
    |   7-8: {order: 0, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   7-9: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   7-10: {order: 1, stereo_parity: null}                              |
    +----------------------------------------------------------------------+
    |   9-11: {order: 1, stereo_parity: null}                              |
    +----------------------------------------------------------------------+
    |   9-12: {order: 1, stereo_parity: null}                              |
    +----------------------------------------------------------------------+
    | reactants keys:                                                      |
    +----------------------------------------------------------------------+
    | - [1, 2, 3, 4, 5, 6]                                                 |
    +----------------------------------------------------------------------+
    | - [7, 8, 9, 10, 11, 12]                                              |
    +----------------------------------------------------------------------+
    | backward TS atoms:                                                   |
    +----------------------------------------------------------------------+
    |   1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}  |
    +----------------------------------------------------------------------+
    |   10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null} |
    +----------------------------------------------------------------------+
    | backward TS bonds:                                                   |
    +----------------------------------------------------------------------+
    |   1-2: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   1-3: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   1-4: {order: 0.9, stereo_parity: null}                             |
    +----------------------------------------------------------------------+
    |   2-5: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   2-6: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   4-7: {order: 0.1, stereo_parity: null}                             |
    +----------------------------------------------------------------------+
    |   7-8: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   7-9: {order: 1, stereo_parity: null}                               |
    +----------------------------------------------------------------------+
    |   7-10: {order: 1, stereo_parity: null}                              |
    +----------------------------------------------------------------------+
    | products keys:                                                       |
    +----------------------------------------------------------------------+
    | - [1, 2, 3, 4, 5, 6]                                                 |
    +----------------------------------------------------------------------+
    | - [7, 8, 9, 10]                                                      |
    +----------------------------------------------------------------------+




PhyDat
^^^^^^^

PhyDat stores physical data and constants used across AutoMech suite.


Transformations
^^^^^^^^^^^^^^^

author: Christoph Gohlke

site: https://www.lfd.uci.edu/~gohlke/

Calculate homogeneous transformation matrices for translating, rotating,
reflecting, scaling, shearing, projecting, orthogonalizing, and superimposing
3D homogeneous coordinates, convert between rotation matrices, Euler angles,
and quaternions, decompose transformation matrices, and provide an Arcball
control object.



Getting Started
~~~~~~~~~~~~~~~
Installation
^^^^^^^^^^^^^
.. code-block:: python

    >>> conda install autochem -c auto-mech

Tutorial
^^^^^^^^
The first step is to make sure the installation was successful by importing some of the modules in each package

.. code-block:: python

    >>> import automol
    >>> import phydat
    >>> import transformations

Then we can move on to using the autochem modules:

* Automol Tutorial\: :ref:`tutorial-hub`
    * :ref:`ioformat-tutorial-doc`
    * :ref:`autoparse-tutorial-doc`
    * :ref:`autoread-tutorial-doc`
* PhyDat tutorial
* Transformations tutorial


Documentation
~~~~~~~~~~~~~
    .. toctree::
        :maxdepth: 3

        submodule_automol
        submodule_phydat
        submodule_transformations
