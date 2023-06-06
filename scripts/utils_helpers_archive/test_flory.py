import numpy as np
import math
# Compares the model we used (6Rg^2/Rf^2 = 6*25/176) vs a more accurate value from renorm theory (6Rg^2/Rf^2 = 0.952)

def flory_self_avoiding(num_residues, num_segments, bead_radius, base_conn_dist=3.6):
    rg = 1.92 * (num_residues ** 0.6)  # based on Kohn PNAS
    a = math.sqrt((rg ** 2) * 176 / 25 / (num_segments ** 1.2))  # from Flory
    # alternatively, assume <(r_i - r_j)^2> = b^2 * |i - j|^(2v) with v = 6/5
    total_length = a * num_segments
    scale = ((total_length - 2 * bead_radius) / num_segments - 2 * bead_radius) / base_conn_dist
    return int(scale)  # truncating to the nearest smaller integer
    
def true_model(num_residues, num_segments, bead_radius, base_conn_dist=3.6):
    rg = 1.92 * (num_residues ** 0.6)  # based on Kohn PNAS
    a = rg / np.sqrt(0.95 / 6) / (num_segments ** 0.6) # from renorm
    # alternatively, assume 6Rg/Rf=sqrt(0.95)
    total_length = a * num_segments
    scale = ((total_length - 2 * bead_radius) / num_segments - 2 * bead_radius) / base_conn_dist
    return round(scale)  # truncating to the nearest smaller integer
    
r = 7.5
num_res = [125, 244, 178, 318]
for nr in num_res:
    print(f'Testing for n_res = {nr}, n_seg = {nr // 20}, radius = {r}:')
    print(f'\t{flory_self_avoiding(nr, nr // 20, r)}')
    print(f'\t{true_model(nr, nr // 20, r)}')