def apply_inlet(u, v, inlet_indices, U_in):
    for (ix, iy) in inlet_indices:
        u[iy, ix] = U_in
        v[iy, ix] = 0.0

def apply_outlet(u, v, outlet_indices):
    pass

def apply_slip_top_bottom(u, v, top_indices, bottom_indices):
    for (ix, iy) in top_indices:
        v[iy, ix] = 0.0
    for (ix, iy) in bottom_indices:
        v[iy, ix] = 0.0

def apply_no_slip(u, v, wall_indices):
    for (ix, iy) in wall_indices:
        u[iy, ix] = 0.0
        v[iy, ix] = 0.0


