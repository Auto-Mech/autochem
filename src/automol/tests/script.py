""" script
"""
import automol
import warnings
warnings.filterwarnings("error")

chi = automol.smiles.chi('[H][H]')
geo = automol.chi.geometry(chi)
