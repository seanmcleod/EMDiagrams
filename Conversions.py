import math
from Atmosphere import ISAtmosphere

def PitotTotalPressure(mach, p):
    if mach < 0:
        return p
    if mach < 1:
        return p * ((1 + 0.2*mach*mach) ** 3.5)
    else:
        #return p * 166.92158 * (mach ** 7) / ((7*mach*mach - 1) ** 2.5)
        return p * 166.92158009316827 * (mach ** 7) / ((7*mach*mach - 1) ** 2.5)
        #factor = (1.2**3.5)*(2.5**2.5)*(2.4**2.5)
        #return p * factor * (mach ** 7) / ((7*mach*mach - 1) ** 2.5)
        
def ImpactPressure(mach, altitude):
    pressure = ISAtmosphere.Pressure(altitude)
    return PitotTotalPressure(mach, pressure) - pressure


def MachtoCAS(mach, altitude):
    pressure = ISAtmosphere.Pressure(altitude)
    sealevel_pressure = ISAtmosphere.sealevel_pressure
    sealevel_density = ISAtmosphere.sealevel_density

    pt = PitotTotalPressure(mach, pressure)
    A = ((pt - pressure)/sealevel_pressure + 1) ** (1/3.5)

    return math.sqrt(7*sealevel_pressure/sealevel_density*(A-1))

def TAStoCAS(tas, altitude):
    mach = ISAtmosphere.Mach(tas, altitude)
    cas = MachtoCAS(mach, altitude)
    return cas

# JSBSim version
def CAStoMach2(cas, altitude):
    pressure = ISAtmosphere.Pressure(altitude)
    pt = pressure + ISAtmosphere.sealevel_pressure*( ((1+cas*cas*ISAtmosphere.sealevel_density/(7*ISAtmosphere.sealevel_pressure)) ** 3.5) - 1)

    if pt/pressure < 1.89223:
        # Mach < 1
        return math.sqrt(5*((pt/pressure)**(1/3.5) - 1))
    else:
        # Mach >= 1
        mach = math.sqrt(0.77666*pt/pressure)
        delta = 1
        target = pt/(166.92158*pressure)

        iter = 0
        while delta > 1e-5 and iter < 10:
            m2 = mach * mach
            m6 = m2 * m2 * m2
            delta = mach*m6/((7*m2-1)**2.5) - target
            diff = 7*m6*(2*m2 - 1)/((7*m2 - 1)**3.5)
            mach -= delta/diff
            iter += 1

        return mach

# JSBSim version
def CAStoMach3(cas, altitude):
    pressure = ISAtmosphere.Pressure(altitude)
    pt = pressure + ISAtmosphere.sealevel_pressure*( ((1+cas*cas*ISAtmosphere.sealevel_density/(7*ISAtmosphere.sealevel_pressure)) ** 3.5) - 1)

    if pt/pressure < 1.89223:
        # Mach < 1
        return math.sqrt(5*((pt/pressure)**(1/3.5) - 1))
    else:
        # Mach >= 1
        asl = ISAtmosphere.sealevel_speed_of_sound
        pt = pressure + ISAtmosphere.sealevel_pressure*(((166.921580* math.pow(cas/asl, 7)) / math.pow(7 * pow(cas/asl, 2) - 1, 2.5)) - 1)
        mach = math.sqrt(0.77666*pt/pressure)
        delta = 1
        target = pt/(166.92158*pressure)

        iter = 0
        while delta > 1e-5 and iter < 10:
            m2 = mach * mach
            m6 = m2 * m2 * m2
            delta = mach*m6/((7*m2-1)**2.5) - target
            diff = 7*m6*(2*m2 - 1)/((7*m2 - 1)**3.5)
            mach -= delta/diff
            iter += 1

        return mach

'''
Based on the formulas in the US Air Force Aircraft Performance Flight Testing Manual (AFFTC-TIH-99-01).
In particular sections 4.6 to 4.8.
The subsonic and supersonic Mach number equations are used with the simple substitutions of (Vc/asl) for M 
and Psl for P. However, the condition for which the equations are used is no longer subsonic (M < 1) or 
supersonic (M > 1) but rather calibrated airspeed being less or greater than the speed of sound ( asl ), 
standard day, sea level (661.48 knots). 
'''
def CAStoMach(cas, altitude):
    pressure = ISAtmosphere.Pressure(altitude)

    if cas < ISAtmosphere.sealevel_speed_of_sound:
        # Bernoulli’s compressible equation (4.11)
        qc = ISAtmosphere.sealevel_pressure * (math.pow(1 + 0.2 * math.pow(cas/ISAtmosphere.sealevel_speed_of_sound, 2), 3.5) - 1)
    else:
        # Rayleigh's supersonic pitot equation (4.16) 
        qc = ISAtmosphere.sealevel_pressure * (((166.9215801 * math.pow(cas/ISAtmosphere.sealevel_speed_of_sound, 7)) / math.pow(7 * math.pow(cas/ISAtmosphere.sealevel_speed_of_sound, 2) - 1, 2.5)) - 1)

    # Solving for M in equation (4.11), also used as initial condition for supersonic case
    mach = math.sqrt(5 * (math.pow(qc/pressure + 1, 2/7) - 1))

    if mach > 1:
        # Iterate equation (4.22) since M appears on both sides of the equation
        for i in range(7):
            mach = 0.88128485 * math.sqrt((qc/pressure + 1) * math.pow(1 - 1/(7*mach*mach), 2.5))

    return mach

def CAStoTAS(cas, altitude):
    mach = CAStoMach(cas, altitude)
    return mach * ISAtmosphere.SpeedOfSound(altitude)

def ktTofps(knots):
    return knots * 1.687664

def fpsTokt(fps):
    return fps / 1.687664


# Bertrand's

def MachFromImpactPressure(qc, p):
  A = qc/p+1
  # Bernoulli’s compressible equation solved for M (4.12)
  # Also used as initial condition for supersonic case
  M = math.sqrt(5.0*(math.pow(A, 1./3.5)-1))

  if M > 1.0:
    for i in range(0, 25):
      # (4.17)
      #M = 0.881285*math.sqrt(A*math.pow(1-1.0/(7.0*M*M),2.5))
      #M = 0.88128485*math.sqrt(A*math.pow(1-1.0/(7.0*M*M),2.5))
      M = 0.8812848543473311 * math.sqrt(A*math.pow(1-1.0/(7.0*M*M),2.5))

  return M

def VcalibratedFromMach(mach, p):
  asl = ISAtmosphere.sealevel_speed_of_sound
  psl = ISAtmosphere.sealevel_pressure
  qc = PitotTotalPressure(mach, p) - p

  return asl * MachFromImpactPressure(qc, psl)

def MachFromVcalibrated(vcas, p):
  asl = ISAtmosphere.sealevel_speed_of_sound
  psl = ISAtmosphere.sealevel_pressure
  qc = PitotTotalPressure(vcas/asl, psl) - psl

  return MachFromImpactPressure(qc, p)


