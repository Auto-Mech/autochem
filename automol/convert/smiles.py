""" smiles conversions
"""
# from xml.etree import ElementTree
# import requests
from automol.convert import _pybel


def inchi(smi):
    """ SMILES => InChI

    OpenBabel appears to work now!! (But it fails for isotopes...)
    """
    rdm = _pybel.from_smiles(smi)
    ich = _pybel.to_inchi(rdm)

    return ich


# def chemspider_inchi(smi):
#     """ convert SMILES to an InChI string using ChemSpider
#
#     Note that this function requires an internet connection!
#
#     I don't understand the code in detail -- sumbled my way through it using
#     stackoverflow.
#     """
#     url = ('http://www.chemspider.com/InChI.asmx/SMILESToInChI?smiles={}'
#            .format(smi))
#
#     # read from the URL
#     req = requests.get(url)
#
#     # interpret the XML; req.content yields the site content as a string
#     # use a try-except to fail gracefully if there's an issue
#     try:
#         ele = ElementTree.fromstring(req.content)
#     except ElementTree.ParseError:
#         raise ValueError("Failed to convert SMILES string '{}' to InChI"
#                          .format(smi))
#
#     # there is just a single value, so we can get the result directly
#     ich = ele.text
#
#     return ich
