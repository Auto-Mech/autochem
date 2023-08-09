"""
Library of physical constants
"""
from qcelemental import constants as qcc


# Physical Constants
NAVO = 6.0221409e+23
RC = 1.98720425864083  # gas constant in cal/(mol.K)
RC_KCAL = 1.98720425864083e-3  # gas constant in kcal/(mol.K)
RC2 = 82.0573660809596  # gas constant in cm^3.atm/(mol.K)
RC_CAL = 1.98720425864083  # gas constant in cal/(mol.K)
RC_ATM = 82.0573660809596  # gas constant in cm^3.atm/(mol.K)
SOL = (qcc.get('speed of light in vacuum') *
       qcc.conversion_factor('meter / second', 'bohr hartree / hbar'))
SOLMS = qcc.get('speed of light in vacuum')
KB = qcc.get('kb')  # The Boltzmann constant (JK$^{-1}$)
H = qcc.get('h')  # The Planck constant (Js)
HBAR = qcc.get('hbar')  # The Planck constant (Js) over 2pi

# Energy Conversion factors
KCAL2CAL = qcc.conversion_factor('kcal/mol', 'cal/mol')
J2CAL = qcc.conversion_factor('J/mol', 'cal/mol')
KJ2CAL = qcc.conversion_factor('kJ/mol', 'cal/mol')
KJ2KCAL = qcc.conversion_factor('kJ/mol', 'kcal/mol')
KJ2EH = qcc.conversion_factor('kJ/mol', 'hartree')
EH2KJ = qcc.conversion_factor('hartree', 'kJ/mol')
EH2J = qcc.conversion_factor('hartree', 'J')
J2EH = qcc.conversion_factor('J', 'hartree')
K2EH = qcc.conversion_factor('kelvin', 'hartree')
KEL2CAL = qcc.conversion_factor('kelvin', 'cal/mol')
WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
KCAL2WAVEN = qcc.conversion_factor('kcal/mol', 'wavenumber')
WAVEN2KEL = 1.4388
KEL2WAVEN = 1.0 / WAVEN2KEL
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')
KCAL2EH = qcc.conversion_factor('kcal/mol', 'hartree')
KCAL2KJ = qcc.conversion_factor('kcal/mol', 'kJ/mol')
WAVEN2EH = qcc.conversion_factor('wavenumber', 'hartree')
EH2WAVEN = qcc.conversion_factor('hartree', 'wavenumber')
EH2KJ = qcc.conversion_factor('hartree', 'kJ/mol')
EH2K = qcc.conversion_factor('hartree', 'kelvin')
WAVEN2J = qcc.conversion_factor('wavenumber', 'hartree')
J2WAVEN = qcc.conversion_factor('hartree', 'wavenumber')

# Length/Angle Conversion Factors
ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')
ANG2M = qcc.conversion_factor('angstrom', 'meter')
M2ANG = qcc.conversion_factor('meter', 'angstrom')
BOHR2M = qcc.conversion_factor('bohr', 'meter')
M2BOHR = qcc.conversion_factor('meter', 'bohr')
BOHR2CM = 5.29177e-9
DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')

# Mass Conversion Factors
# AMU2KG = qcc.conversion_factor('atomic mass unit', 'kilogram')
AMU2EMASS = qcc.conversion_factor("atomic_mass_unit", "electron_mass")
AMU2KG = 1.66053892173e-27
AMU2AU = 1.0 / (9.1093837015e-28 * 6.022e23)

# Pressure Conversion Factors
ATM2TORR = qcc.conversion_factor('atmosphere', 'torr')
TORR2ATM = qcc.conversion_factor('torr', 'atmosphere')
ATM2BAR = qcc.conversion_factor('atmosphere', 'bar')
BAR2ATM = qcc.conversion_factor('bar', 'atmosphere')

# Other Factors
JIFFY = 2.441888e-17
