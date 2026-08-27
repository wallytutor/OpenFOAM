import math

f    = 2/360    # Fraction of total pipe
D    = 50e-6    # sizeDistribution.value
rho  = 2500     # constantProperties.rho0
mdot = 0.01     # mass flow rate
N    = 100      # nParticle

n = mdot / (rho * math.pi * (D**3) / 6)

print(f"parcelsPerSecond {f * n / N:.0f};")
print(f"nParticle        {N:.0f};")