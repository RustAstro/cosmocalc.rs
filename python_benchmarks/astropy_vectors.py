from astropy.cosmology import FLRW, FlatLambdaCDM, LambdaCDM, Planck18

def flat_universe_distances_no_relativistic_contribution():
    H0 = 69.6
    Om0 = 0.286
    Ode0 = 0.714
    Ob0 = 0.05
    cosmo = FlatLambdaCDM(H0, Om0, Ode0, 0, Ob0=Ob0)

    print('comoving transverse distance to z=3')
    print(cosmo.comoving_transverse_distance(3))
    # 6482.549296743232 Mpc
    print('angular diameter distance to z=3')
    print(cosmo.angular_diameter_distance(3))
    # 1620.637324185808 Mpc
    print('luminosity distance to z=3')
    print(cosmo.luminosity_distance(3))
    # 25930.197186972928 Mpc


def flat_universe_distances_with_radiation_but_no_neutrinos():
    H0 = 69.6
    Om0 = 0.299
    Ode0 = 0.7
    Ob0 = 0.05
    Tcmb0 = 2.7255
    cosmo = FlatLambdaCDM(H0, Om0, Ode0, 2.7255, 0, Ob0=Ob0)

    print('comoving transverse distance to z=3')
    print(cosmo.comoving_transverse_distance(3))
    # 6398.504909397802 Mpc
    print('angular diameter distance to z=3')
    print(cosmo.angular_diameter_distance(3))
    # 1599.6262273494506 Mpc
    print('luminosity distance to z=3')
    print(cosmo.luminosity_distance(3))
    # 25594.01963759121 Mpc


def lookback_time():
   Om0 = 0.27
   Ode0 = 0.73
   Ob0 = 0.044
   H0 = 70.0
   cosmo = FlatLambdaCDM(H0, Om0, Ode0, 2.7255, 0, Ob0=Ob0)
   print('lookback time to z=3')
   print(cosmo.lookback_time(3))
   # 11.641584565145552 Gyr


def comoving_volume():
   Om0 = 0.27
   Ode0 = 0.73
   Ob0 = 0.044
   H0 = 70.0
   cosmo = FlatLambdaCDM(H0, Om0, Ode0, 2.7255, 0, Ob0=Ob0)
   print('hubble distance')
   print(cosmo.hubble_distance)
   # 4282.749400000001 Mpc
   print('comoving volume to z=3')
   print(cosmo.comoving_volume(3))
   # 1179374442443.6943 Mpc3


if __name__=="__main__":
    flat_universe_distances_no_relativistic_contribution()
    flat_universe_distances_with_radiation_but_no_neutrinos()
    lookback_time()
    comoving_volume()
