using LinearAlgebra, SPICE

export e⃗, hyp_turn_angle, hyp_periapsis, hyp_exit_v⃗, hyp_exit_r⃗

GM_CB(μ_CB_or_CB_name) = typeof(μ_CB_or_CB_name) != String ? μ_CB_or_CB_name : bodvrd(μ_CB_or_CB_name,"GM")[1] # if GM provided directly, use it, else retrieve from body name

function e⃗(;r⃗,v⃗,μ_CB_or_CB_name)
    μ_CB = GM_CB(μ_CB_or_CB_name)
    return ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
end

function hyp_turn_angle(;e)
    return 2*asin(1/e) # Vallado 4e Eq. 2-28 (p53)
end

function hyp_periapsis(;v⃗∞,turn_angle,μ_CB_or_CB_name)
    μ_CB = GM_CB(μ_CB_or_CB_name)
    return μ_CB/norm(v⃗∞)^2 * (1/cos((π-turn_angle)/2)-1) # Vallado 4e Eq. 12-12 (p959)
end

# h⃗(;r⃗,v⃗) = cross(r⃗, v⃗) # Specific angular momentum
# ĥ(;r⃗,v⃗) = normalize(h⃗(r⃗=r⃗,v⃗=v⃗)) # Specific angular momentum unit vector

function hyp_exit_r⃗(;r⃗∞,v⃗∞,μ_CB_or_CB_name)
    μ_CB = GM_CB(μ_CB_or_CB_name)
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name) # undefined for circular but we don't care here
    ν̃ = acos((ecc⋅r⃗∞)/(norm(ecc)*norm(r⃗∞))) # Vallado 4e Eq. 2-86, except ignore correction for values over π since all we want is the angle between r⃗ and periapsis to perform the rotation

    return axisar(
        normalize(cross(r⃗∞,v⃗∞)), # construct axis from specific angular momentum vector
        2*ν̃ # Rotate by "true anomaly" to essentially mirror about the ĥ-ê plane
    ) * r⃗∞
end

function hyp_exit_v⃗(;r⃗∞,v⃗∞,μ_CB_or_CB_name)
    return axisar(
    normalize(cross(r⃗∞,v⃗∞)), # construct axis from specific angular momentum vector
    hyp_turn_angle(e=norm(e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name)))
    ) * v⃗∞
end

function hyp_exit_x⃗(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = x⃗∞[1:3]
    v⃗∞ = x⃗∞[4:6]
    exit_x⃗ = Float64[] # we know SPICE.axisar will return double-precision (Float64)
    push!(exit_x⃗, hyp_exit_r⃗(;r⃗∞,v⃗∞,μ_CB_or_CB_name))
    push!(exit_x⃗, hyp_exit_v⃗(;r⃗∞,v⃗∞,μ_CB_or_CB_name))
    return exit_x⃗
end