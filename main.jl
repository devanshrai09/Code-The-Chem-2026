using Plots
include("physics.jl")

hi = 1000.0
ho = 500.0

ri = 0.005      # 5 mm inner radius
ro = 0.006      # 6 mm outer radius

L = 20.0        # 2.71 m heat exchanger
kw = 15.0        # W/(m·K)

Ti = 90.0        # °C, hot water
To = 20.0        # °C, cold water

m_dot_inner = 0.01       # kg/s
Cp_inner = 4180.0       # J/(kg·K)

m_dot_outer = 0.015       # kg/s
Cp_outer = 4180.0       # J/(kg·K)

N = 200                 # 200 segments
Ri, Rwall, Ro = cylindrical_resistances(
    hi,
    ho,
    ri,
    ro,
    L,
    kw
)

println("Inner convection resistance = ", Ri)
println("Wall conduction resistance = ", Rwall)
println("Outer convection resistance = ", Ro)

Rtotal = total_thermal_resistance(Ri, Rwall, Ro)

println("Total thermal resistance = ", Rtotal, " K/W")

Ui, Uo = overall_heat_transfer_coefficient(
    Rtotal,
    ri,
    ro,
    L
)

println("Overall U based on inner area = ", Ui, " W/(m²·K)")
println("Overall U based on outer area = ", Uo, " W/(m²·K)")

Ai = 2 * π * ri * L
Ao = 2 * π * ro * L

println("Ui * Ai = ", Ui * Ai) #verification that UiAi=UoAo
println("Uo * Ao = ", Uo * Ao)  




Q, Twi, Two = heat_transfer_and_wall_temperatures(
    Ti,
    To,
    Ri,
    Rwall,
    Ro
)

println()
println("Heat transfer rate = ", Q, " W")
println("Inner wall temperature = ", Twi, " °C")
println("Outer wall temperature = ", Two, " °C")





Ci = heat_capacity_rate(
    m_dot_inner,
    Cp_inner
)

Co = heat_capacity_rate(
    m_dot_outer,
    Cp_outer
)

println()
println("Inner fluid heat capacity rate = ", Ci, " W/K")
println("Outer fluid heat capacity rate = ", Co, " W/K")


Ti_out, To_out = outlet_temperatures(
    Ti,
    To,
    Q,
    Ci,
    Co
)  #whenever this function is called then for a segment ti_out is different for ti_in cuz usne heat lose ya gain ki hai due to radial heat transfer

println()
println("Inner fluid outlet temperature = ", Ti_out, " °C")
println("Outer fluid outlet temperature = ", To_out, " °C")

Q_inner, Q_outer, error_inner, error_outer =
    energy_balance_check(
        Ti,
        Ti_out,
        To,
        To_out,
        Ci,
        Co,
        Q
    )

println()
println("----- Energy Balance Check -----")

println("Calculated Q       = ", Q, " W")
println("Heat lost by inner = ", Q_inner, " W")
println("Heat gained by outer = ", Q_outer, " W")

println()
println("Inner-side error = ", error_inner, " W")
println("Outer-side error = ", error_outer, " W")




ΔL = segment_length(L, N)

println()
println("Number of segments = ", N)
println("Length of each segment = ", ΔL, " m")





Ti_next, To_next, Twi, Two, Qsegment =
    simulate_segment(
        Ti,
        To,
        hi,
        ho,
        ri,
        ro,
        ΔL,
        kw,
        Ci,
        Co
    )

println()
println("----- One Segment Results -----")

println("Segment heat transfer = ", Qsegment, " W")

println("Inner fluid before = ", Ti, " °C")
println("Inner fluid after  = ", Ti_next, " °C")

println()

println("Outer fluid before = ", To, " °C")
println("Outer fluid after  = ", To_next, " °C")

println()

println("Inner wall temperature = ", Twi, " °C")
println("Outer wall temperature = ", Two, " °C")


Ti_profile,
To_profile,
Twi_profile,
Two_profile,
Q_profile = simulate_parallel_flow(
    Ti,
    To,
    hi,
    ho,
    ri,
    ro,
    L,
    kw,
    Ci,
    Co,
    N
)  
println()
println("----- Full Heat Exchanger Results -----")

println("Inner fluid inlet  = ", Ti_profile[1], " °C")
println("Inner fluid outlet = ", Ti_profile[end], " °C")

println()

println("Outer fluid inlet  = ", To_profile[1], " °C")
println("Outer fluid outlet = ", To_profile[end], " °C")

println()

println("Total heat transferred = ", sum(Q_profile), " W") #note that this value is smaller than the one we calculated 
 #assuming L as whole length because we thought about only the inlet and out temp is the drivng force but in reality ◬T decreases as we go across x direction
 #The segmented result is lower because the temperature difference is not actually constant throughout the exchanger.

 x = segment_positions(L, N)

println()
println("First position = ", x[1], " m")
println("Last position = ", x[end], " m")
println("Number of points = ", length(x))


Q_inner,
Q_outer,
Q_segments,
error_inner,
error_outer = validate_simulation(
    Ti_profile,
    To_profile,
    Ci,
    Co,
    Q_profile
)

println()
println("----- Full Simulation Validation -----")

println("Heat lost by inner fluid = ", Q_inner, " W")
println("Heat gained by outer fluid = ", Q_outer, " W")
println("Sum of segment heat transfer = ", Q_segments, " W")

println()
println("Inner error = ", error_inner, " W")
println("Outer error = ", error_outer, " W")


segment = 50

Twi = Twi_profile[segment]
Qsegment = Q_profile[segment]

r, T_radial = radial_temperature_profile(
    Twi,
    Qsegment,
    ri,
    ro,
    kw,
    ΔL,
    100
)
println()
println("----- Radial Temperature Validation -----")

println("Calculated inner wall temperature = ", T_radial[1])
println("Stored inner wall temperature     = ", Twi)

println("Calculated outer wall temperature = ", T_radial[end])
println("Stored outer wall temperature     = ", Two_profile[segment])


# Convert radius from meters to millimeters
r_mm = r .* 1000

# Create radial temperature plot
p_radial = plot(
    r_mm,
    T_radial,
    label = "Temperature Through Tube Wall",
    xlabel = "Radius (mm)",
    ylabel = "Temperature (°C)",
    title = "Radial Temperature Profile - Segment $(segment)",
    linewidth = 3,
    marker = :none
)

# Mark inner and outer wall surfaces
scatter!(
    p_radial,
    [ri * 1000, ro * 1000],
    [Twi, Two_profile[segment]],
    label = "Inner and Outer Wall Surfaces",
    markersize = 6
)

display(p_radial)

x_mid = segment_midpoints(L, N)

# --------------------------------
# AXIAL TEMPERATURE PROFILE
# --------------------------------

p_axial = plot(
    x,
    Ti_profile,
    label = "Inner Fluid",
    xlabel = "Heat Exchanger Length (m)",
    ylabel = "Temperature (°C)",
    title = "Axial Temperature Profile",
    linewidth = 3
)

# Outer fluid
plot!(
    p_axial,
    x,
    To_profile,
    label = "Outer Fluid",
    linewidth = 3
)

# Inner wall
plot!(
    p_axial,
    x_mid,
    Twi_profile,
    label = "Inner Wall",
    linewidth = 2,
    linestyle = :dash
)

# Outer wall
plot!(
    p_axial,
    x_mid,
    Two_profile,
    label = "Outer Wall",
    linewidth = 2,
    linestyle = :dash
)

display(p_axial)


# combined graph for radial as well as axial graph
p_combined = plot(
    p_radial,
    p_axial,
    layout = (1, 2),
    size = (1200, 500)
)

display(p_combined)