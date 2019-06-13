""" test autorxn.trans
"""
import automol
import autorxn
from autorxn import trans


def test__migration():
    """ test trans.migration
    """
    rcts_ich = 'InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H'
    prds_ich = 'InChI=1S/C5H5O/c1-2-3-4-5-6/h2-4H,1H2'

    rcts_cgr = automol.inchi.graph(rcts_ich)
    prds_cgr = automol.inchi.graph(prds_ich)

    tra = trans.migration(rcts_cgr, prds_cgr)
    inv_tra = autorxn.trans.invert(tra)

    trans_rcts_cgr = autorxn.trans.apply(tra, rcts_cgr)
    assert trans_rcts_cgr == prds_cgr

    trans_prds_cgr = autorxn.trans.apply(inv_tra, prds_cgr)
    assert trans_prds_cgr == rcts_cgr


def test__addition():
    """ test trans.addition
    """
    rcts_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H'
    prds_ich = 'InChI=1S/C4H5F2O/c5-3-1-2-4(6)7/h1-4,7H'

    rcts_cgr = automol.inchi.graph(rcts_ich)
    prds_cgr = automol.inchi.graph(prds_ich)

    tra = trans.addition(rcts_cgr, prds_cgr)
    inv_tra = autorxn.trans.invert(tra)

    trans_rcts_cgr = autorxn.trans.apply(tra, rcts_cgr)
    assert trans_rcts_cgr == prds_cgr

    trans_prds_cgr = autorxn.trans.apply(inv_tra, prds_cgr)
    assert trans_prds_cgr == rcts_cgr

    rcts_sgr = automol.graph.stereomers(rcts_cgr)[0]
    prds_sgrs = autorxn.trans.stereo_compatible_products(rcts_sgr, tra)
    assert len(prds_sgrs) == 4


def test__abstraction():
    """ test trans.abstraction
    """
    rcts_ich = 'InChI=1S/C7H14.HO/c1-3-5-7-6-4-2;/h3,5H,4,6-7H2,1-2H3;1H'
    prds_ich = 'InChI=1S/C7H13.H2O/c1-3-5-7-6-4-2;/h3,5H,1,4,6-7H2,2H3;1H2'

    rcts_cgr = automol.inchi.graph(rcts_ich)
    prds_cgr = automol.inchi.graph(prds_ich)

    tra = trans.abstraction(rcts_cgr, prds_cgr)
    inv_tra = autorxn.trans.invert(tra)

    trans_rcts_cgr = autorxn.trans.apply(tra, rcts_cgr)
    assert trans_rcts_cgr == prds_cgr

    trans_prds_cgr = autorxn.trans.apply(inv_tra, prds_cgr)
    assert trans_prds_cgr == rcts_cgr


if __name__ == '__main__':
    test__migration()
    test__addition()
    test__abstraction()
