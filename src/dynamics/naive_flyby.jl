using LinearAlgebra, SPICE

export e⃗, hyp_turn_angle, hyp_periapsis, hyp_exit_v⃗, hyp_exit_r⃗

GM_CB(μ_CB_or_CB_name) = typeof(μ_CB_or_CB_name) != String ? μ_CB_or_CB_name : bodvrd(μ_CB_or_CB_name,"GM")[1] # if GM provided directly (is a number), use it, else retrieve from body name (String)

function e⃗(;r⃗,v⃗,μ_CB_or_CB_name)
    μ_CB = GM_CB(μ_CB_or_CB_name)
    return ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
end

function ν(;r⃗,v⃗,μ_CB_or_CB_name)
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name) # undefined for circular but we don't care here
    ν̃ = acos((ecc⋅r⃗∞)/(norm(ecc)*norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    return r⃗⋅v⃗ > 0 ? ν̃ : 360 - ν̃ # Correct for halfspace
end

function semi_major_axis(;r⃗,v⃗,μ_CB_or_CB_name)
    return 1/(2/norm(r⃗) - norm(v⃗)^2/GM_CB(μ_CB_or_CB_name)) # Vallado 4e Eq. 2-74 (p96)
end

function hyp_anom(;r⃗,v⃗,μ_CB_or_CB_name)
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name) # Copy-paste true anomaly function here to save one internal variable (ecc) and hopefully improve performance there
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗∞)/(e*norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    trueanom = r⃗⋅v⃗ > 0 ? ν̃ : 360 - ν̃ # Correct for halfspace
    return 2*atanh(sqrt((e-1)/(e+1)) * tan(trueanom/2))
end

function hyp_turn_angle(;e)
    return 2*asin(1/e) # Vallado 4e Eq. 2-28 (p53)
end

function hyp_periapsis(;v⃗∞,turn_angle,μ_CB_or_CB_name)
    return GM_CB(μ_CB_or_CB_name)/norm(v⃗∞)^2 * (1/cos((π-turn_angle)/2)-1) # Vallado 4e Eq. 12-12 (p959)
end

# h⃗(;r⃗,v⃗) = cross(r⃗, v⃗) # Specific angular momentum
# ĥ(;r⃗,v⃗) = normalize(h⃗(r⃗=r⃗,v⃗=v⃗)) # Specific angular momentum unit vector

function hyp_exit_r⃗(;r⃗∞,v⃗∞,μ_CB_or_CB_name)
    return vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
        r⃗∞, # To be rotated
        e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name), # periapsis vector (not required to be unit vector)
        π # 180 degrees
    )
end

function hyp_exit_v⃗(;r⃗∞,v⃗∞,μ_CB_or_CB_name)
    return vrotv(
        v⃗∞,
        cross(r⃗∞,v⃗∞), # axis of rotation as specific angular momentum vector
        hyp_turn_angle(e=norm(e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name))) # turn by the hyperbolic turn angle
    )
end

function hyp_exit_x⃗(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = x⃗∞[1:3]
    v⃗∞ = x⃗∞[4:6]
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name)
    exit_x⃗ = Vector{Float64}(undef,6)
    exit_x⃗[1:3] = vrotv( # Mirror the radius vector "at infinity" about the periapsis vector
    r⃗∞, # To be rotated
    ecc, # periapsis vector (not required to be unit vector)
    π # 180 degrees
    )
    exit_x⃗[4:6] = vrotv(
        v⃗∞,
        cross(r⃗∞,v⃗∞), # axis of rotation as specific angular momentum vector
        hyp_turn_angle(e=norm(ecc)) # turn by the hyperbolic turn angle
    )
    return exit_x⃗
end

function hyp_TOF(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = x⃗∞[1:3]
    v⃗∞ = x⃗∞[4:6]
    μ_CB = GM_CB(μ_CB_or_CB_name)
    ecc = ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗∞)/(e*norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    trueanom_in = r⃗⋅v⃗ > 0 ? ν̃ : 360 - ν̃ # Correct for halfspace; true anomaly at SOI entry
    H_in = 2*atanh(sqrt((e-1)/(e+1)) * tan(trueanom_in/2)) # Hyperbolic anomaly at SOI entry
    H_out = -H_in # We know that the hyperbolic anomaly has the same magnitude at entry and exit
    sma = 1/(2/norm(r⃗∞) - norm(v⃗∞)^2/μ_CB) # Vallado 4e Eq. 2-74 (p96)
    return sqrt(-sma^3/μ_CB)*(e*sinh(H_out) - H_out - (e*sinh(H_in) - H_in))
end