import math

__author__ = 'Todd Bergman'

"""
THIS SCRIPT DOES NOT ANALYZE GLOBAL STABILITY OF WALL SETTLEMENT

This script calculates the demand capcity ration for all external stability limiting states: Eccentricity,
Overturning, Sliding, and Bearing Restistance. 

Assumptions:
- The wall has an embedment depth of 2 ft
- The water depth below the leveling pad is at least 30ft [from RKCI]

All Outputs = ['Strength 1a', 'Strength 1b-1', 'Strength IV-1'] - Each being a separate load combination
 """

'Definitions'
W1 = []     # Reinforced Fill Weight including pavement DL (lb/ft)
W2 = []     # Sloped backfill weight over reinforced area (lbs/ft)
W3 = []     # Flat backfill weight over reinforced area (lbs/ft)
Mr = []     # Sum of Resisting Moments without Live Load (lb-ft/ft)
Mr2 = []    # Sum of Resisting Moments including live load (lbs-ft/ft)
M_0 = []    # Sum of Overturning Moments (lb-ft/ft)
e = []      # Eccentricity for Overturning (ft) (FLDOT)
Fr = []     # Resisting Force (Sum of factored resisting forces)(lb/ft)
Fd = []     # Driving Force (Sum of factored horizontal forces) (lb/ft)
q1v = []    # Surcharge Vertical Force Over Reinforced Area (lbs/ft)
Fv = []     # Retained Backfill Vertical Resultant (lbs/ft)
Fh = []     # Retained Backfill Horizontal Resultant (lbs/ft)
qv = []     # Total resultant vertical surcharge force (lbs/ft)
qh = []     # Total resultant horizontal surcharge force (lbs/ft)
Rv = []     # Sum of factored vertical forces acting within reinforced soil mass without live load (q1L) used in sliding
# CDR calculation (lbs/ft)
Rv2 = []    # Sum of factored vertical forces acting within (lbs/ft)
e2 = []     # Eccentricity for Bearing (ft)
qvb = []    # Factored Uniform Bearing Pressure (lbs/sq. ft) [AASHTO 11.6.3.2-1] [FHWA NHI-10-024 Eq 4-18]
qn = []     # Nominal Bearing Resistance of a soil layer (lbs/sq. ft) [AASHTO EQ. 10.6.3.1.2-a-1]
L_eff = []  # Effective Foundation Width (ft)

Overturning_CDR = []    # Overturning Capacity-Demand Ratio
Sliding_CDR = []        # Sliding Capacity-Demand Ratio
bearing_CDR = []        # Bearing Capacity-Demand Ratio

'INPUT PARAMETERS'

'Wall Geometry'
H_0 = 10        # Height of Wall above grade (ft) - Field Observation
D = 2           # Wall Embedment Depth (ft) - Assumption based on design standards
H = H_0 + D     # Total Height of the wall including embeddment
L = 8.5         # Reinforcement Strap Length - Field Observation
beta = 26.2     # beta is zero with a Horizontal Backslope (degrees)
lam = 16        # Distance from the back of the wall to the to of the slope (ft)
d = 30          # water depth below leveling pad (ft) - From RKCI Geotech Engineer

'Soil Parameters'
UW_rf = 125     # Unit Weight of Reinforced Soil - 125 pcf for Cement Stabilized per TxDOT Design Sheet
UW_bf = 125     # Unit Weight of Backfill Soil - 125 pcf per TxDOT Design Sheet
phi_bf = 45     # Backfill soil angle of internal friction  - 45 degrees for Cement stabilized per TxDot Design Sheet
wc_fs = 0.15    # Water content at Foundation Soil (%)

DUW_fs = 105                    # Dry Unit Weight of Soil at Foundation (pcf) - From RKCI Boring Log for 1960-NE-1
UW_fs = DUW_fs * (1 + wc_fs)    # Unit Weight of Foundation Soil - Calculated from RKCI's reported Dry Unit Weight
print('Unit Weight of Foundation Soil: ', UW_fs, '[Calculation Based on Dry Unit Weight from RKCI Borings]\n')

phi_fs = 30     # Foundation soil angle of friction - 30 degrees per AASHTO standards [AASHTO 11.10.5.3 - Sliding]
c_fs = 250      # Cohesion of Foundation Soil - 200 psf per RKCI Addendum Letter to Geotechnical Report 01
phi_u = min(phi_fs, phi_bf)      # Base angle of internal friction for sliding - lesser of phi_rf or phi_fs [AASHTO 11.10.5.3 - Sliding]

if d > 1.5 * L:
    CW = 1  # Modification based on ground water level [AASHTO Table 10.6.3.1.2a.2]
else:
    CW = 0.5 * d / 3 / L

'Applied Loading'
q1 = 450    # Uniform Surcharge load over reinforced soil mass (psf) - PAVEMENT 150 pcf [AASHTO 3.11.6.1]
q2 = 450    # Surcharge load behind reinforced soil mass (psf)

'Load Factors [Table 3.4.1-1 AASHTO 2007]'
lf_EV = (1, 1.35, 1.35)     # EARTH VERTICAL Strength 1a, Strength 1b, Strength IV [Table 3.4.1-2 AASHTO 2007]
lf_EH = (1.5, 1.5, 1.5)     # EARTH HORIZONTAL Strength 1a, Strength 1b, Strength IV [Table 3.4.1-2 AASHTO 2007]
lf_LL = (1.75, 1.75, 0)     # LIVE LOAD VERTICAL Strength 1a & Strength 1b [Table 3.4.1-1 AASHTO 2007]
lf_LLh = (1.75, 1.75, 0)    # LIVE LOAD HORIZONTAL Strength 1a & Strength 1b [Table 3.4.1-1 AASHTO 2007]
lf_ESV = (.75, .75, .75)    # EARTH SURCHARGE VERTICAL (Min LF - load increases stability) [Table 3.4.1-2 AASHTO 2007]
lf_ESH = (1.5, 1.5, 1.5)    # EARTH SURCHARGE HORIZONTAL (Max LF - load decreases stability) [Table 3.4.1-2 AASHTO 2007]

'Resistance Factors [Table 11.5.6-1 AASHTO 2007'
phi_sliding = 1.0   # External Stability Resistance Factors for Sliding
phi_bearing = 0.65  # External Stability Resistance Factors for Bearing Resistance

e_max = L / 4  # Max eccentricity for MSE walls with soil foundations. [AASHTO 11.6.6.3 (i.e. e_max = L/3)]

'Bearing Resistance Factors'
if phi_fs == 0:
    Nc = 5.14
else:
    Nc = (((math.exp(math.pi * math.tan(math.radians(phi_fs)))) *
           (math.tan(math.pi / 4 + (math.radians(phi_fs)) / 2)) ** 2) - 1) * (1 / (math.tan(math.radians(phi_fs))))

Nq = ((math.exp(math.pi * math.tan(math.radians(phi_fs)))) * (math.tan(math.pi / 4 + (math.radians(phi_fs)) / 2)) ** 2)

Ng = (((math.exp(math.pi * math.tan(math.radians(phi_fs)))) *
       (math.tan(math.pi / 4 + (math.radians(phi_fs)) / 2)) ** 2) + 1) * 2 * (math.tan(math.radians(phi_fs)))

print('Nc =', Nc)           # Cohesion Term (undrained loading) Bearing Capacity Factor [AASHTO Table 10.6.3.1.2a-1]
print('Nq =', Nq)           # Surcharge Term Bearing Capacity Factor [AASHTO Table 10.6.3.1.2a-1]
print('Ng =', Ng, '\n')     # Unit Weight Term Bearing Capacity Factor [AASHTO Table 10.6.3.1.2a-1]

if beta == 0:

    h = '---'
    I = '---'
    W2 = '---'
    W3 = '---'
    Fv = '---'
    Fh = '---'
    qv = '---'
    qh = '---'

    ' Active Coefficient of Earth Pressure for Vertical Wall with Horizontal Backslopes '
    Kab = math.tan(math.radians(45 - phi_bf / 2)) ** 2  # Active Coefficient of Earth Pressure [Eq.4-1  FHWA NHI-10-024]
    Ft = 0.5 * Kab * UW_rf * H ** 2  # Retained Backfill Resultant
    qt = Kab * q2 * H  # Uniform Surcharge Resultant FLDOT

    for x in range(3):
        q1v.append(q1 * L * lf_LL[x])
        W1.append(UW_rf * H * L * lf_EV[x])
        Rv = W1
        Rv2.append(Rv[x] + q1v[x])
        Fd.append(Ft * lf_EH[x] + qt * lf_LLh[x])
        Fr.append(phi_sliding * Rv[x] * math.tan(math.radians(phi_u)))
        Mr.append(W1[x] * (L / 2))
        Mr2.append(Mr[x] + (q1 * L * (L / 2)) * lf_LLh[x])
        M_0.append(lf_EH[x] * Ft * (H / 3) + lf_LLh[x] * qt * (H / 2))

else:

    if lam >= L:
        I = '---'
        W3 = '---'

        h = H + (L * math.tan(math.radians(beta)))  # Wall height for backfill stress calculations (ft)
        ' Active Coefficient of Earth Pressure for Vertical Walls with Surcharge Slope '
        # gamma for active coefficient of earth pressure calculation
        gamma = (1 + math.sqrt((math.sin(math.radians(phi_bf + beta)) * math.sin(math.radians(phi_bf - beta))) / (
            math.sin(math.radians(90 + beta)) * math.sin(math.radians(90 - beta))))) ** 2
        # Active Coefficient of Earth Pressure
        Kab = (math.sin(math.radians(90 + phi_bf)) ** 2) / (
            gamma * (math.sin(math.radians(90)) ** 2) * (math.sin(math.radians(90 - beta))))
        qt = Kab * q2 * h  # Uniform Surcharge Resultant FLDOT
        Ft = 0.5 * Kab * UW_rf * h ** 2  # Retained Backfill Resultant

        for x in range(3):
            q1v.append(q1 * L * lf_LL[x])
            W1.append(UW_rf * H * L * lf_EV[x])
            W2.append((0.5 * UW_bf * (h - H) * L * lf_EV[x]))
            Fv.append(lf_EH[x] * Ft * math.cos(math.radians(90 - beta)))
            Fh.append(lf_EH[x] * Ft * math.cos(math.radians(beta)))
            qv.append(lf_LL[x] * qt * math.cos(math.radians(90 - beta)))
            qh.append(lf_LL[x] * qt * math.cos(math.radians(beta)))
            Rv.append(W1[x] + W2[x] + Fv[x] + qv[x])
            Rv2.append(W1[x] + W2[x] + Fv[x] + q1v[x] + qv[x])
            Mr.append(W1[x] * (L / 2) + W2[x] * (L * 2 / 3) + Fv[x] * L + qv[x] * L)
            Mr2.append(Mr[x] + (q1v[x] * (L / 2)))
            M_0.append(Fh[x] * (h / 3) + qh[x] * h / 2)

    else:
        h = H + (lam * math.tan(math.radians(beta)))  # Wall height for backfill stress calculations
        I = math.degrees(math.atan((h - H) / (2 * H)))  # Resultant Earth Pressure Inclination (deg)
        ' Active Coefficient of Earth Pressure for Vertical Walls with Broken Backslope '
        # gamma for active coefficient of earth pressure calculation
        gamma = (1 + math.sqrt((math.sin(math.radians(phi_bf + I)) * math.sin(math.radians(phi_bf - I))) / (
            math.sin(math.radians(90 + I)) * math.sin(math.radians(90 - I))))) ** 2
        # Active Coefficient of Earth Pressure
        Kab = (math.sin(math.radians(90 + phi_bf)) ** 2) / (
            gamma * (math.sin(math.radians(90)) ** 2) * (math.sin(math.radians(90 - I))))
        qt = Kab * q2 * h  # Uniform Surcharge Resultant FLDOT
        Ft = 0.5 * Kab * UW_rf * h ** 2  # Retained Backfill Resultant

        for x in range(3):
            q1v.append(q2 * (L - lam) * lf_LL[x])
            W1.append(UW_rf * H * L * lf_EV[x])
            W2.append((0.5 * UW_rf * (lam * math.tan(math.radians(beta))) * lam * lf_EV[x]))
            W3.append((UW_rf * (lam * math.tan(math.radians(beta))) * (L - lam) * lf_EV[x]))
            Fv.append(lf_EH[x] * Ft * math.cos(math.radians(90 - I)))
            Fh.append(lf_EH[x] * Ft * math.cos(math.radians(I)))
            qv.append(lf_LL[x] * qt * math.cos(math.radians(90 - I)))
            qh.append(lf_LL[x] * qt * math.cos(math.radians(I)))
            Rv.append(W1[x] + W2[x] + W3[x] + Fv[x] + q1v[x] + qv[x])
            Rv2 = Rv
            Mr.append((W1[x] * (L / 2)) + (W2[x] * (lam * 2 / 3)) + (W3[x] * (lam + (0.5 * (L - lam)))) + (Fv[x] * L) +
                      (qv[x] * L))
            Mr2.append(Mr[x] + (q1v[x] * (lam + (0.5 * (L - lam)))))
            M_0.append(Fh[x] * (h / 3) + qh[x] * h / 2)

    for x in range(3):
        Fd.append(Fh[x] + qh[x])

for x in range(3):
    Fr.append(Rv[x] * math.tan(math.radians(phi_u)))
    e.append((L / 2) - ((Mr[x] - M_0[x]) / Rv[x]))
    e2.append((L / 2) - ((Mr2[x] - M_0[x]) / Rv2[x]))
    L_eff.append(L - 2 * e2[x])
    # Factored Uniform Bearing Pressure [AASHTO 11.6.3.2-1] [FHWA NHI-10-024 Eq 4-18]
    qvb.append(Rv2[x] / (L - 2 * e2[x]))
    # Nominal Bearing Resistance of a soil layer [AASHTO EQ. 10.6.3.1.2-a-1]
    qn.append(phi_bearing * ((c_fs * Nc) + UW_fs * D * Nq + (0.5 * (L_eff[x]) * UW_fs * Ng * CW)))
    Overturning_CDR.append(Mr[x] / M_0[x])
    Sliding_CDR.append(Fr[x] / Fd[x])
    bearing_CDR.append(qn[x] / qvb[x])

Eccentricity_CDR = (max(e) / e_max)
Overturning_CDR = min(Overturning_CDR)
Sliding_CDR = min(Sliding_CDR)
bearing_CDR = min(bearing_CDR)

print('qvb =', qvb, 'Factored Uniform Bearing Pressure (lbs/sq. ft)')
print('qn =', qn, 'Nominal Bearing Resistance of a soil layer (lbs/sq. ft)')
print('\nh = ', h, 'Effective Wall Height (for Surcharge Slope/Broken Backslope) (ft)')
print('\nW1 = ', W1, 'Reinforced Fill Weight including pavement DL (lb/ft)')
print('W2 = ', W2, 'Sloped backfill weight over reinforced area (lbs/ft)')
print('W3 = ', W3, 'Flat backfill weight over reinforced area (lbs/ft)')
print('\nq1v =', q1v, 'Surcharge Vertical Force Over Reinforced Area (lbs/ft)')
print('\nalpha =', I, 'Resultant Earth Pressure Inclination (deg)')
print('Ft = ', Ft, 'Retained Backfill Resultant (lbs/ft)')
print('Fv = ', Fv, 'Retained Backfill Vertical Resultant (lbs/ft)')
print('Fh = ', Fh, 'Retained Backfill Horizontal Resultant (lbs/ft)')
print('\nqt = ', qt, 'Uniform Surcharge Resultant (lbs/ft)')
print('qv = ', qv, 'Total resultant vertical surcharge force (lbs/ft)')
print('qh = ', qh, 'Total resultant horizontal surcharge force (lbs/ft)')
print('\nFd =', Fd, 'Driving Force (Sum of factored horizontal forces) (lb/ft)')
print('Fr =', Fr, 'Resisting Force (Sum of factored resisting forces)(lb/ft)')
print('\nRv =', Rv, 'Sum of factored vertical forces acting within reinforced soil mass without live load (q1L) used '
                    'in sliding CDR calculation (lbs/ft)')
print('Rv2 =', Rv2, 'Sum of factored vertical forces acting within (lbs/ft)')
print('\nMr =', Mr, 'Sum of Resisting Moments without Live Load (lb-ft/ft)')
print('Mr2 =', Mr2, 'Sum of Resisting Moments including live load (lbs-ft/ft)')
print('M0 =', M_0, 'Sum of Overturning Moments (lb-ft/ft)')
print('\ne =', e, 'Eccentricity for Overturning (ft)')
print('e2 =', e2, 'Eccentricity for Bearing (ft)')
print('L_eff =', L_eff, 'Effective Foundation Width (ft)')
print('e_max =', e_max, 'Max eccentricity for MSE walls with soil foundations (ft)')
print('\nKab = ', Kab, 'Active Coefficient of Earth Pressure')
print('\nEccentricity_CDR =', Eccentricity_CDR)
print('Overturning_CDR =', Overturning_CDR)
print('Sliding_CDR =', Sliding_CDR)
print('Bearing_CDR =', bearing_CDR)

'ECCENTRICITY'
print('\nECCENTRICITY')
if Eccentricity_CDR <= 1:
    print('Wall Passes Limiting Eccentricity with a CDR of', Eccentricity_CDR, ', which is less than 1.\n')
else:
    print('Wall Fails in Limiting Eccentricity with a CDR of,', Eccentricity_CDR, '\n')

'OVERTURNING'
print('OVERTURNING')
if Overturning_CDR >= 1:
    print('Wall Passes Overturning with a CDR of', Overturning_CDR, ', which is greater than 1.\n')
else:
    print('Wall Fails in Overturning\n')

'SLIDING'
print('SLIDING')
if Sliding_CDR >= 1:
    print('Wall Passes in Sliding with a CDR of', Sliding_CDR, ', which is greater than 1.\n')
else:
    print('Wall Fails in Sliding with a CDR of,', Sliding_CDR, '\n')

'BEARING RESISTANCE'
print('BEARING RESISTANCE')
if bearing_CDR >= 1:
    print('Wall Passes in Bearing with a CDR of', bearing_CDR, ', which is greater than 1.\n')
else:
    print('Wall Fails in Bearing with a CDR of,', bearing_CDR, '\n')
