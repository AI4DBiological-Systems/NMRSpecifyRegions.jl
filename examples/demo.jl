
import NMRHamiltonian
import NMRSignalSimulator

import NMRSpecifyRegions

using DataDeps
import Tar

using LinearAlgebra
import PyPlot
#import JSON3

include("./helpers/data.jl")
include("./helpers/SH.jl")
include("./helpers/utils.jl")

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

#molecule_entries = ["L-Methionine"; "L-Phenylalanine"; "DSS"; "Ethanol"; "L-Isoleucine"]
#molecule_entries = ["alpha-D-Glucose"; "beta-D-Glucose"; "DSS"; "D2O"]
molecule_entries = ["alpha-D-Glucose"; "DSS"]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462


u_offset = 0.2 #in units ppm.
Δcs_padding = 0.02 #in units ppm.
min_window_cs = 0.06 #in units ppm.

λ0 = 4.0

### end inputs.

### set up.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

As, Rs = runSH(molecule_entries)

#dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)
dummy_SSFID = NMRSignalSimulator.SpinSysParamsType2(0.0)
# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(4.0)

Bs = NMRSignalSimulator.fitclproxies(As, dummy_SSFID, λ0;
    names = molecule_entries)

### end set up.


## frequency locations. For plotting.
ΩS_ppm = getPsnospininfo(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(combinevectors(ΩS_ppm))


u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

# This is the frequency range that we shall work with.
P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

## get intervals.
ΩS0 = getΩS(As)
ΩS0_ppm = getPs(ΩS0, hz2ppmfunc)

Δcs_padding = 0.1 # units are in ppm.
Δsys_cs = initializeΔsyscs(As, Δcs_padding)
exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_entries, ΩS0_ppm, Δsys_cs; min_dist = Δcs_padding)



q = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs)

# evaluate at the plotting positions.
y = q.(U_rad)


band_inds, band_inds_set = NMRSpecifyRegions.getcostinds(exp_info, P)

U_cost = U[band_inds]
P_cost = P[band_inds]

y_cost = y[band_inds]


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(y), label = "data spectrum")
PyPlot.plot(P_cost, real.(y_cost), "x", label = "positions in band_inds")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("positions against data spectrum, real part")
