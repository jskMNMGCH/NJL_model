# Physical constants
const c = 2.99792458e10  # [cm/sec]
const hbarc = 197.32698044404107 # [MeV fm]
const rho0 = 0.17  # fm^-3

# unit conversion factor
const dyncm2_to_MeVfm3 = 1.0/(1.602176634e33)
const MeVfm3_to_dyncm2 = 1.602176634e33
const gcm3_to_MeVfm3 = 1.0/(1.7826619216278976e12)
const MeVfm3_to_gcm3 = 1.7826619216278976e12  # 1e(9+39)*1.602e-19/(3e8)**2  正しくは MeV/fm3/c2 to g/cm3

#  Parameters for the NJL-type quark model.
global Nf = 2.0
global Nc = 3.0
global Lambda = 587.9  # MeV
global Gs = 2.44/Lambda^2
global Gv = Gs
global m = 5.6  # MeV

#  Parameters for the original NJL-model.
global Nf_o = 2.0
global Nc_o = 3.0
global Lambda_o = 587.9  # MeV
global G = 2.44/Lambda_o^2
global m_o = 5.6  # MeV