.. _run_glossary::

Run Glossary
=============

.. raw:: html

    <link rel="stylesheet" href="_static/css/glossary.css">




--------------------------------------------------------------------------------------

The run.dat file contains the input for workflow.  This includes specifications for

 - general `input`_ (e.g., run/save locations, type of mechanism)
 - species and reactions `numbers`_ of that mechanisms the computations should be done on
 - the electronic structure `tasks`_ to be run and with what options
 - the transport `tasks`_  to be run and with what options
 - the thermo or k(T, P)  `tasks`_ to be run and with what options
 - the information we want to run processing `tasks`_ on to output csv and txt files

--------------------------------------------------------------------------------------

|

.. _input:

Glossary of Input Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The input block has the following format

.. code-block:: python

    input
        run_prefix = My_AutoMech_Database/RUN
        save_prefix = My_AutoMech_Database/SAVE
    end input
|
|

.. list-table:: Input Keywords
   :widths: 10 20
   :header-rows: 1

   * - Input Keyword
     - Description
   * - run_prefix\*
     - path to run directory
   * - save_prefix\*
     - path to save directory
   * - inp_mech
     - the type of mechanism in mechanism.dat (chemkin, rmg, etc)
   * - out_mech
     - the type of mechanism to be outputted (chemkin, rmg, etc)
   * - inp_spc
     - the format of the species dictionary (csv)
   * - out_spc
     - the format of outputted species dictionary (csv)

*compulsory keywords

|

.. _numbers:

Glossary of SPC/PES Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spc and pes blocks have no keywords.  The spc block require only that a species number
(or comma seperated list of species numbers) is given that corresponds to the species in
the species.csv that the tasks should be performed on. The pes block, similarly, selects
for the potential energy surface number in the sorted mechanism to run on, and can be
broken down further by selecting specific channels.

.. code-block:: python

    pes
        1: 1, 2, 3
    end pes

.. code-block:: python

    spc
        1,2,3,4,5
    end spc


|
|

.. _tasks:

Glossary of Task Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The remaining blocks (els, trans, thermo, and ktp) take on a similar format to one another.
For instance, an example electronic structure block (els):

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
    end els

Here els is the block name (compulsory), init_geom is the task keyword (compulsory),
spc is the type of molecule (type is compulsory), and runlvl is one of the options of this
task keyword (in this case the run level is compulsory, but many options are optional).

The glossary for these tasks is below, where the Task Keyword can be clicked on for more `detail`_.
Compulsory information is identified with an asterisk.

.. list-table:: Task Keywords
   :widths: 8 10 9 12
   :header-rows: 1

   * - Task Keyword
     - Description
     - Type to Run on\*
     - Options
   * -
     -
     -
     -

.. list-table:: els
   :widths: 10 20 10 20
   :header-rows: 0

   * - `init_geom`_\* *\*compulsory for new species only*
     - initialize a geometry for a new species
     - spc
     - runlvl\*, retryfail, overwrite
   * - `conf_samp`_
     - search for additional conformers
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, cnf_range
   * - `conf_opt`_
     - runs an optimization job on any number of conformers
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `conf_energy`_
     - runs a single point energy for any number of conformers
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `conf_grad`_
     - runs a gradient computation on any number of conformers
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `conf_hess`_
     - runs a hessian computation on any number of conformers
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `conf_vpt2`_
     - runs an vpt2 anharmonic analysis on any number of conformers
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `conf_prop`_
     -
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `hr_scan`_
     - runs a hindered rotor scan on the lowest energy conformer
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite, tors_model
   * - `hr_reopt`_
     - runs a geometry optimization for each step of a hindered rotor
       scan using the geometry optimized at an inplvl of theory
       for that dihedral angle
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite, tors_model
   * - `hr_grad`_
     - runs gradient computations along the steps of a hindered rotor
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite, tors_model
   * - `hr_hess`_
     - runs hessian computations along the steps of a hindered rotor
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite, tors_model
   * - `hr_energy`_
     - runs single point energy  computations along the steps of a hindered rotor
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite, tors_model
   * - `hr_vpt2`_
     - runs vpt2  computations along the steps of a hindered rotor
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite, tors_model
   * - `find_ts`_
     - search for a transition state
     - ts
     - runlvl\*, inplvl\*, retryfail, overwrite, nobarrier,
       var_splvl1, var_splvl2, var_scnlvl
   * - `conf_pucker`_ (dev)
     - search for any ring puckering conformations (in development)
     - spc, ts
     - runlvl\*, inplvl\*, retryfail, overwrite, conf_range
   * - `tau_samp`_
     - sample addition configurations that dont need to be local energy minima
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite
   * - `tau_energy`_
     - run single point energy for tau geometries
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite
   * - `tau_grad`_
     - run gradient computations for tau geometries
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite
   * - `tau_hess`_
     - run hessian computations for tau geometries
     - spc, ts, all
     - runlvl\*, inplvl\*, retryfail, overwrite

.. list-table:: trans
   :widths: 10 20 10 20
   :header-rows: 0

   * - `onedmin`_
     -
     - spc, ts, all
     - inplvl\*, runlvl\*, retryfail, overwrite,
       bath, nsamp, njobs, smin, smax,
       conf, pot
   * -
     -
     -
     -

.. list-table:: thermo
   :widths: 10 20 10 20
   :header-rows: 0

   * - `write_mess_thermo`_
     - write the MESS partition function input file for each species
     - *no type prefix for this section*
     - kin_model, spc_model, overwrite
   * - `run_mess_thermo`_
     - run MESS for each species
     - *no type prefix for this section*
     - kin_model, spc_model, overwrite, inpname
   * - `run_fits_thermo`_
     - produce NASA polynomials and CHEMKIN style inputs for each speices
     - *no type prefix for this section*
     - kin_model

.. list-table:: ktp
   :widths: 10 20 10 20
   :header-rows: 0

   * - `write_mess`_
     - write the MESS rate constant input file for each connected PES
     - *no type prefix for this section*
     - kin_model, spc_model, overwrite
   * - `run_mess`_
     - run MESS for each connected PES
     - *no type prefix for this section*
     - kin_model, spc_model, overwrite, inpname
   * - `run_fits`_
     - produce Arhennius fits and CHEMKIN style input for the rate constants
     - *no type prefix for this section*
     - kin_model

.. list-table:: process
   :widths: 10 20 10 20
   :header-rows: 0

   * - `freqs`_
     - produce a csv file of frequencies
     - spc, ts, vdw, all
     - geolvl, proplvl, nconfs, econfs, scale
   * - `energy`_
     - produce a csv file of energies
     - spc, ts, vdw, all
     - geolvl, proplvl, nconfs, econfs
   * - `geo`_
     - produce a txt file of geometries
     - spc, ts, vdw, all
     - geolvl, proplvl, nconfs, econfs
   * - `zmatrix`_
     - produce a txt file of zmatrices
     - spc, ts, vdw, all
     - geolvl, proplvl, nconfs, econfs
   * - `enthalpy`_
     - produce a csv file of 0 K heats of formation
     - spc, ts, vdw, all
     - geolvl, proplvl, nconfs, econfs
   * - `coeffs`_
     - produce a csv file of the reference molecules required
       for a heat of fomration calculation
     - spc, ts, vdw, all
     - *None*

|
|
|

.. _detail:

Detailed Description of Task Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _init_geom:

init_geom
^^^^^^^^^^

This tasks initializes a geometry for your species based on a geometry built from the default force fields from, on first attempt,
RDKit, and on second attempt, OpenBabel.  It then optimizes it with the level of theory specified by runlvl.  This is a key that
connects it to a theory in the theory.dat file. The geometry is then saved in the filesystem


.. list-table::
   :class: options
   :widths: 10 55 25 10
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
     - Default
   * - **runlvl**\*
     - the task will be performed at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **retryfail**
     - if the electronic structure job fails, reads the failed output, processes why it may have failed,  and resubmits it.
     - **True** or **False**
     - **False**
   * - **overwrite**
     - overwrite any existing or running data saved for this molecule and level of theory in the filesystem
     - **True** or **False**
     - **False**

**Example**:

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
    end els


Check out our (some sort of theory manual) to see how to set up wb97_sto3g as a theory key in theory.dat.

|

.. _conf_samp:

conf_samp
^^^^^^^^^^^^
This tasks runs a monte carlo sampling over the torsional coordinates of a molecule to produce geometry samples that
are then optimized at the runlvl of theory.  The starting geometry to generate new conformers is a stored geometry
that was optimized at inplvl of theory already.  If no such geometry is saved in the saved filesystems,
the user needs to have run init_geom for this quantum chemistry method and basis set. The number of samples is
set with the cnf_range keyword.

.. list-table::
   :class: options
   :widths: 10 55 25 10
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
     - Default
   * - **inplvl**\*
     - the starting geometry is from prior optimization at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **runlvl**\*
     - the task will be performed at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **retryfail**
     - if the electronic structure job fails, reads the failed output, processes why it may have failed,  and resubmits it.
     - **True** or **False**
     - **False**
   * - **overwrite**
     - overwrite any existing or running data saved for this molecule and level of theory in the filesystem
     - **True** or **False**
     - **False**
   * - **cnf_range**
     - the number of geometries to generated by the stochastic sampling routine. The value can be either an integer,
       and that exact number of samples will be optimized or an array of [a, b, c, d] which will take N (the number
       of torsions) and use the formula nsamp = min(d, a + ( b * c^N))
     - <int> or

       <[int, int, int, int]>
     - [3, 1, 3, 100]

**Example**:

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
        spc  conf_samp     runlvl=m062x_ccpvdz inplvl=wb97_sto3g  cnf_range=[6,1,3,200]
    end els

Check out our (some sort of theory manual) to see how to set up wb97_sto3g as a theory key in theory.dat.

|

.. _conf_opt:

conf_opt
^^^^^^^^^

This task runs geometry optimizations on a set of conformers for a molecule, specified by conf_range, at the
runlvl of theory, where the starting geometries and sorting is done based on the inplvl of theory.

.. list-table::
   :class: options
   :widths: 10 50 30 10
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
     - Default
   * - **inplvl**\*
     - the starting geometry is from prior optimization at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **runlvl**\*
     - the task will be performed at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **retryfail**
     - if the electronic structure job fails, reads the failed output, processes why it may have failed,  and resubmits it.
     - **True** or **False**
     - **False**
   * - **overwrite**
     - overwrite any existing or running data saved for this molecule and level of theory in the filesystem
     - **True** or **False**
     - **False**
   * - **cnf_range**
     - the conformers to run the task on, sorted using the inplvl of theory
     - **min**: the lowest energy conformer

       **n** <int>: the lowest <int> number of conformers

       **e** <float>: any conformers within <float> kcal/mol of the lowest energy conformer

     - **min**

**Example**:

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
        spc  conf_samp     runlvl=wb97_sto3g   inplvl=wb97_sto3g
        spc  conf_opt       runlvl=m062x_ccpvdz inplvl=wb97_sto3g  cnf_range=n5
    end els


|

.. _conf_energy:

conf_energy
^^^^^^^^^^^

This task runs single point energies on a set of conformers for a molecule, specified by conf_range, at the
runlvl of theory, where the starting geometries and sorting is done based on the inplvl of theory.

.. list-table::
   :class: options
   :widths: 10 50 30 10
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
     - Default
   * - **inplvl**\*
     - the starting geometry is from prior optimization at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **runlvl**\*
     - the task will be performed at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **retryfail**
     - if the electronic structure job fails, reads the failed output, processes why it may have failed,  and resubmits it.
     - **True** or **False**
     - **False**
   * - **overwrite**
     - overwrite any existing or running data saved for this molecule and level of theory in the filesystem
     - **True** or **False**
     - **False**
   * - **cnf_range**
     - the conformers to run the task on, sorted using the inplvl of theory
     - **min**: the lowest energy conformer

       **n** <int>: the lowest <int> number of conformers

       **e** <float>: any conformers within <float> kcal/mol of the lowest energy conformer

     - **min**

**Example**:

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
        spc  conf_samp     runlvl=wb97_sto3g   inplvl=wb97_sto3g
        spc  conf_energy   runlvl=m062x_ccpvdz inplvl=wb97_sto3g  cnf_range=n5
    end els

|

.. _conf_grad:

conf_grad
^^^^^^^^^

This task runs gradients on a set of conformers for a molecule, specified by conf_range, at the
runlvl of theory, where the starting geometries and sorting is done based on the inplvl of theory.

.. list-table::
   :class: options
   :widths: 10 50 30 10
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
     - Default
   * - **inplvl**\*
     - the starting geometry is from prior optimization at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **runlvl**\*
     - the task will be performed at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **retryfail**
     - if the electronic structure job fails, reads the failed output, processes why it may have failed,  and resubmits it.
     - **True** or **False**
     - **False**
   * - **overwrite**
     - overwrite any existing or running data saved for this molecule and level of theory in the filesystem
     - **True** or **False**
     - **False**
   * - **cnf_range**
     - the conformers to run the task on, sorted using the inplvl of theory
     - **min**: the lowest energy conformer

       **n** <int>: the lowest <int> number of conformers

       **e** <float>: any conformers within <float> kcal/mol of the lowest energy conformer

     - **min**

**Example**:

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
        spc  conf_samp     runlvl=wb97_sto3g   inplvl=wb97_sto3g
        spc  conf_grad     runlvl=wb97_sto3g   inplvl=wb97_sto3g  cnf_range=n5
    end els


|

.. _conf_hess:

conf_hess
^^^^^^^^^

This task runs hessians on a set of conformers for a molecule, specified by conf_range, at the
runlvl of theory, where the starting geometries and sorting is done based on the inplvl of theory.

.. list-table::
   :class: options
   :widths: 10 50 30 10
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
     - Default
   * - **inplvl**\*
     - the starting geometry is from prior optimization at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **runlvl**\*
     - the task will be performed at this level of theory
     - <theory level key> *as defined in theory.dat*
     - No Default
   * - **retryfail**
     - if the electronic structure job fails, reads the failed output, processes why it may have failed,  and resubmits it.
     - **True** or **False**
     - **False**
   * - **overwrite**
     - overwrite any existing or running data saved for this molecule and level of theory in the filesystem
     - **True** or **False**
     - **False**
   * - **cnf_range**
     - the conformers to run the task on, sorted using the inplvl of theory
     - **min**: the lowest energy conformer

       **n** <int>: the lowest <int> number of conformers

       **e** <float>: any conformers within <float> kcal/mol of the lowest energy conformer

     - **min**

**Example**:

.. code-block:: python

    els
        spc  init_geom     runlvl=wb97_sto3g
        spc  conf_samp     runlvl=wb97_sto3g   inplvl=wb97_sto3g
        spc  conf_hess     runlvl=wb97_sto3g   inplvl=wb97_sto3g  cnf_range=n5
    end els

|

.. _conf_vpt2:

conf_vpt2
^^^^^^^^^^^^

lalala

.. _conf_prop:

conf_prop
^^^^^^^^^^^^

lalala

.. _conf_pucker:

conf_pucker
^^^^^^^^^^^^

lalala


.. _hr_scan:

hr_scan
^^^^^^^^^

Runs a hindered rotor scan on either a species (at an energy minima) or a transition state (at a saddle point).
The hindered rotor scans can be run with multiple models, which are described in the options table below.


.. list-table::
   :class: options
   :widths: 10 25 65
   :header-rows: 1

   * - Option Keyword
     - Description
     - Values
   * - **inplvl**
     - the starting geometry comes from this level of theory
     - <theory level> *as defined in theory.dat*
   * - **runlvl**
     - the task will be performed at this level of theory
     - <theory level> *as defined in theory.dat*
   * - **tors_model**
     - the type of optimization run at each step along the torsional profile
     - **1dhr**: Scans are along one torsional coordinate at a time, that coordinate is frozen, and all other coordinates are optimized

       **1dhrfa**: Scans are along one torsional coordinate at a time, all coordinates are frozen

       **1dhrf**: ?????????????????

       **mdhr**: 2 or 3 torsional coordinates are scanned together to define a rotor, those coordinates are frozen, and all other coordinates are optimized

       **mdhrv**: ?????????????????


.. _hr_grad:

hr_grad
^^^^^^^^^


.. _hr_hess:

hr_hess
^^^^^^^^^

.. _hr_energy:

hr_energy
^^^^^^^^^


.. _hr_vpt2:

hr_vpt2
^^^^^^^^^


.. _hr_reopt:

hr_opt
^^^^^^^^^


.. _find_ts:

find_ts
^^^^^^^^^
lalala

.. _tau_samp:

tau_samp
^^^^^^^^^


.. _tau_grad:

tau_grad
^^^^^^^^^


.. _tau_energy:

tau_energy
^^^^^^^^^^

.. _tau_hess:

tau_hess
^^^^^^^^^
