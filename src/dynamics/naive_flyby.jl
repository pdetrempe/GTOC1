using LinearAlgebra, SPICE

export e⃗, hyp_turn_angle, hyp_periapsis, h⃗, hyp_exit_v

# all r⃗ and v⃗ in planet-centered inertial coordinates

e⃗(;r⃗,v⃗,μ_CB) = ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
e⃗(;r⃗,v⃗,CB) = e⃗(r⃗=r⃗,v⃗=v⃗,μ_CB=SPICE.bodvrd(CB,"GM")[1])

hyp_turn_angle(;e) = 2*asin(1/e) # Vallado 4e Eq. 2-28 (p53)
hyp_turn_angle(;r⃗,v⃗,μ_CB) = hyp_turn_angle(e=norm(e⃗(r⃗=r⃗,v⃗=v⃗,μ_CB=μ_CB)))
hyp_turn_angle(;r⃗,v⃗,CB) = hyp_turn_angle(r⃗=r⃗,v⃗=v⃗,μ_CB=SPICE.bodvrd(CB,"GM")[1])

hyp_periapsis(;v⃗∞,turn_angle,μ_CB) = μ_CB/norm(v⃗∞)^2 * (1/cos((π-turn_angle)/2)-1) # Vallado 4e Eq. 12-12 (p959)
hyp_periapsis(;r⃗∞,v⃗∞,μ_CB) = hyp_periapsis(v⃗∞=v⃗∞,turn_angle=hyp_turn_angle(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB=μ_CB),μ_CB=μ_CB)
hyp_periapsis(;v⃗∞,e,μ_CB) = hyp_periapsis(v⃗∞=v⃗∞,turn_angle=hyp_turn_angle(e=e),μ_CB=μ_CB)
hyp_periapsis(;v⃗∞,turn_angle,CB) = hyp_periapsis(v⃗∞=v⃗∞,turn_angle=turn_angle,μ_CB=SPICE.bodvrd(CB,"GM")[1])
hyp_periapsis(;r⃗∞,v⃗∞,CB) = hyp_periapsis(r⃗∞=r⃗∞,v⃗∞=v⃗∞,μ_CB=SPICE.bodvrd(CB,"GM")[1])
hyp_periapsis(;v⃗∞,e,CB) = hyp_periapsis(v⃗∞=v⃗∞,e=e,μ_CB=SPICE.bodvrd(CB,"GM")[1])

h⃗(;r⃗,v⃗) = cross(r⃗, v⃗) # Specific angular momentum

hyp_exit_v(;r⃗∞,v⃗∞,μ_CB) = SPICE.axisar(normalize(h⃗(r⃗=r⃗∞,v⃗=v⃗∞)),hyp_turn_angle(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB=μ_CB)) * v⃗∞
hyp_exit_v(;r⃗∞,v⃗∞,CB) = hyp_exit_v(r⃗∞=r⃗∞,v⃗∞=v⃗∞,μ_CB=SPICE.bodvrd(CB,"GM")[1])

# TODO: add wrapper function to check periapsis against planet surface or "keep-out" and sphere of influence