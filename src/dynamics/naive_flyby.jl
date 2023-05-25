using LinearAlgebra, SPICE

export e⃗, ν, a, sphere_of_influence, hyp_anom, hyp_turn_angle, hyp_periapsis, hyp_exit_r⃗, hyp_exit_v⃗, hyp_exit_x⃗, hyp_TOF

get_GM(μ_CB_or_CB_name) = typeof(μ_CB_or_CB_name) != String ? μ_CB_or_CB_name : bodvrd(μ_CB_or_CB_name,"GM")[1] # if GM provided directly (is a number), use it, else retrieve from body name (String)

function e⃗(;r⃗,v⃗,μ_CB_or_CB_name)
    μ_CB = get_GM(μ_CB_or_CB_name)
    return ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
end

function ν(;r⃗,v⃗,μ_CB_or_CB_name)
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name)
    ν̃ = acos((ecc⋅r⃗∞)/(norm(ecc)*norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    return r⃗⋅v⃗ > 0 ? ν̃ : 360 - ν̃ # Correct for halfspace
end

function a(;r⃗,v⃗,μ_CB_or_CB_name)
    return 1/(2/norm(r⃗) - norm(v⃗)^2/get_GM(μ_CB_or_CB_name)) # Vallado 4e Eq. 2-74 (p96)
end

function i(;r⃗,v⃗)
    h⃗ = r⃗ × v⃗
    return acos(([0,0,1]⋅normalize(h⃗))) # Vallado 4e Eq. 2-82 (p99)
end

function Ω(;r⃗,v⃗)
    h⃗ = r⃗ × v⃗
    n⃗ = [0,0,1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    RAAN = acos([1,0,0]⋅normalize(n⃗)) # Vallado 4e Eq. 2-84 (p99)
    return n⃗[2] > 0 ? RAAN : 360 - RAAN
end

function ω(;r⃗,v⃗,μ_CB_or_CB_name)
    h⃗ = r⃗ × v⃗
    n⃗ = [0,0,1] × h⃗ # Vallado 4e Eq. 2-83 (p99)
    ecc = e⃗(r⃗=r⃗∞,v⃗=v⃗∞,μ_CB_or_CB_name=μ_CB_or_CB_name)
    AOP = acos((n⃗⋅ecc)/(norm(n⃗)*norm(ecc))) # Vallado 4e Eq. 2-85 (p100)
    return ecc[3] > 0 ? AOP : 360 - AOP
end

function sphere_of_influence(;CB::String,orbiting_body::String)
    orbiting_body_state = spkgeo(bodn2c(orbiting_body),0,base_ref_frame,bodn2c(CB))[1]
    orbiting_body_GM = bodvrd(orbiting_body,"GM")[1]
    CB_GM = bodvrd(CB,"GM")[1]
    orbiting_body_a = 1/(2/norm(orbiting_body_state[1:3]) - norm(orbiting_body_state[4:6])^2/CB_GM) # Vallado 4e Eq. 2-74 (p96)
    return orbiting_body_a*(orbiting_body_GM/CB_GM)^0.4 # Vallado 4e Eq. 12-2 (p948)
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
    return get_GM(μ_CB_or_CB_name)/norm(v⃗∞)^2 * (1/cos((π-turn_angle)/2)-1) # Vallado 4e Eq. 12-12 (p959)
end

# h⃗(;r⃗,v⃗) = r⃗ × v⃗ # Specific angular momentum
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
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
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
        r⃗∞ × v⃗∞, # axis of rotation as specific angular momentum vector
        hyp_turn_angle(e=norm(ecc)) # turn by the hyperbolic turn angle
    )
    return exit_x⃗
end

function hyp_TOF(;x⃗∞,μ_CB_or_CB_name)
    r⃗∞ = x⃗∞[1:3]
    v⃗∞ = x⃗∞[4:6]
    μ_CB = get_GM(μ_CB_or_CB_name)
    ecc = ((norm(v⃗)^2 - μ_CB/norm(r⃗))*r⃗ - (r⃗⋅v⃗)*v⃗)/μ_CB # Vallado 4e Eq. 2-78 (p98)
    e = norm(ecc)
    ν̃ = acos((ecc⋅r⃗∞)/(e*norm(r⃗∞))) # Vallado 4e Eq. 2-86 (p100)
    trueanom_in = r⃗⋅v⃗ > 0 ? ν̃ : 360 - ν̃ # Correct for halfspace; true anomaly at SOI entry
    H_in = 2*atanh(sqrt((e-1)/(e+1)) * tan(trueanom_in/2)) # Hyperbolic anomaly at SOI entry
    H_out = -H_in # We know that the hyperbolic anomaly has the same magnitude at entry and exit
    sma = 1/(2/norm(r⃗∞) - norm(v⃗∞)^2/μ_CB) # Vallado 4e Eq. 2-74 (p96)
    return sqrt(-sma^3/μ_CB)*(e*sinh(H_out) - H_out - (e*sinh(H_in) - H_in))
end