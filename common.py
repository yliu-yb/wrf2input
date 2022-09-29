
def temperature_from_potential_temperature(potential_temperature, pressure):
    # pressure units:hPa
    # potential_temperature units:K
    # temperature units:K
    return potential_temperature / (1000 / pressure)**0.286

def calc_density(temperature, pressure):
    # temperature units:K
    # pressure units:Pa
    # density units:kgm-3
    Rd = 287 # jkg-1k-1
    return pressure / (Rd * temperature)