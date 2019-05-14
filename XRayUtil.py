import xraylib
import numpy as np

class Spectrum(object):
  """Object to describe spectrum"""

  def __init__(self, energy, intensity, label):
    """
    :param energy: Energy points [keV]
    :param intensity: Intensity [a.u.]
    :param label: Spectrum label
    """
    super(Spectrum, self).__init__()
    if not len(intensity) == len(energy):
      raise ValueError(
        'Parameters energy and intensity must have equal length')
    self.energy = energy
    self.intensity = intensity
    self.label = label

  def AddMatter(self, chemical_name, density, l, label=None):
    """
    Add absorption matter

    :param chemical_name: Chemical formula of compund to be paresd in xraylib
    :param density: density [g/cm3]
    :param l: thickness [cm]
    :param label: matter label
    """

    cross_section = []
    for e in self.energy:
      cross_section.append(xraylib.CS_Total_CP(chemical_name, e))
    mu = np.array(cross_section) * density
    absorp = np.exp(-mu * (l))
    self.intensity = self.intensity * absorp

    if label is None:
      self.label = '{} + {} matter, {} cm'.format(self.label, chemical_name, l)
    else:
      self.label = '{} + {}'.format(self.label, label)

class Detector(object):
  """Object to describe CCD detector"""

  def __init__(self, model):
    """
    :param model: model of the detector
    """
    super(Detector, self).__init__()
    if model == 'HAMTSU_101U':
      self.efficiency = 0.126433511211
      self.scintillator_matter = 'Gd2O2S'
      self.scintillator_density = 7.3
      self.l = 10e-4
      self.label = 'model: HAMTSU_101U'
    elif  model == 'HAMTSU_102U':
      self.efficiency = 0.126433511211
      self.scintillator_matter = 'Gd2O2S'
      self.scintillator_density = 7.3
      self.scintillator_thickness = 20e-4
      self.label = 'model: HAMTSU_102U'
    elif model == 'XIMEA_11':
      self.efficiency = 6.55605
      self.scintillator_matter = 'Gd2O2S'
      self.scintillator_density = 7.3
      self.scintillator_thickness = 22e-4
      self.label = 'madel: XIMEA_11'
    elif model == 'VARIAN_2520DX':
      self.efficiency = 0.559224396812
      self.scintillator_matter = 'CsI'
      self.scintillator_density = 4.51
      self.scintillator_thickness = 22e-4
      self.label = 'model: VARIAN_2520DX'
    else:
      raise ValueError(
        'Not implemented yet')

def GenerateSpectrum(i, u, z, num_points=100, spec_type='XRAYTUBE'):
  """
  Intensity of Characteristic Lines by formula R=BIA(UA-UK)1,5

  :param i: (float): Current in [a]
  :param u: (float): Voltage [keV]
  :param z: (int): The sequence number of the item
  :param num_points: number of points
  :param spec_type: type of XRay spevtrum can be one of 'XRAYTUBE', 'BRELUNG' or 'CHARLINES' 
  Returns:
    objec Spectrum
  """

  widht = {24:(1.97, 2.39, 70, 30), #Cr
           29:(2.4, 2.98, 100, 50), # Cu
           42:(6.42, 6.66, 500, 150), # Mo
           47:(8.6, 8.9, 150, 50), # Ag
           74:(8.6, 8.9, 500, 200)} # W

  energies = np.linspace(u * 1e-3, u, num_points)
  spectrum = np.zeros_like(energies)

  if spec_type == 'BRELUNG' or spec_type == 'XRAYTUBE':
    c = 3.0
    lambda_0 = 12.398 / u  # Critical wavelength [A]
    lambdas = 12.398 / energies
    spectrum += c * c / lambda_0 * i * 0.001 * z * (lambdas - lambda_0) / lambdas ** 3
  
  if spec_type == 'CHARLINES' or spec_type == 'XRAYTUBE':
    w1, w2, k1, k2 = widht[z]
    w1 *= 2e-1
    w2 *= 2e-1

    spec = np.zeros_like(energies)
  
    energy = xraylib.LineEnergy(z, xraylib.KA1_LINE)
    intensity = k1 * i * 0.001 * np.power(u - energy, 1.5)
    indexes = np.where((energies < energy + w1) & (energies > energy - w1))
    int1 = np.linspace(0.0, intensity, len(indexes[0]) / 2)
    int2 = np.linspace(intensity, 0.0, len(indexes[0]) - len(indexes[0]) / 2)
    spec[indexes] = np.concatenate((int1, int2), axis=0)

    energy = xraylib.LineEnergy(z, xraylib.KA2_LINE)
    intensity = k1 * i * 0.001 * np.power(u - energy, 1.5)
    indexes = np.where((energies < energy + w2) & (energies > energy - w2))
    int1 = np.linspace(0.0, intensity, len(indexes[0]) / 2)
    int2 = np.linspace(intensity, 0.0, len(indexes[0]) - len(indexes[0]) / 2)
    spec[indexes] = np.concatenate((int1, int2), axis=0)

    energy = xraylib.LineEnergy(z, xraylib.KB1_LINE)
    intensity = k2 * i * 0.001 * np.power(u - energy, 1.5)
    indexes = np.where((energies < energy + w1) & (energies > energy - w1))
    int1 = np.linspace(0.0, intensity, len(indexes[0]) / 2)
    int2 = np.linspace(intensity, 0.0, len(indexes[0]) - len(indexes[0]) / 2)
    spec[indexes] = np.concatenate((int1, int2), axis=0)

    #energy = xraylib.LineEnergy(z, xraylib.KB2_LINE)
    #intensity = k2 * i * 0.001 * np.power(u - energy, 1.5)
    #indexes = np.where((energies < energy + w2) & (energies > energy - w2))
    #int1 = np.linspace(0.0, intensity, len(indexes[0]) / 2)
    #int2 = np.linspace(intensity, 0.0, len(indexes[0]) - len(indexes[0]) / 2)
    #spec[indexes] = np.concatenate((int1, int2), axis=0)

    spectrum += spec

  return Spectrum(energy=energies, intensity=spectrum, label='I={} mA, U={} kV, Z={}'.format(i, u, z))

def CreateDetector(model):
  """
  Return created object of detector

  :param model: name of detector model 
  :return: detector
  """
  return Detector(model)

def GetRegisteredValue(spectrum, detector):
  """
  Return registred signal

  :param spectrum: X-Ray spectra
  :param detector: X-Ray detector
  :return: registred signal intensity
  """
  ab_eff = []
  for e in spectrum.energy:
    cross_section = xraylib.CS_Total_CP(detector.scintillator_matter, e)
    ab_eff.append(1.0 - np.exp(- cross_section * 
                  detector.scintillator_density * 
                  detector.scintillator_thickness))

  reg_val = np.sum(ab_eff * spectrum.intensity) / len(ab_eff)
  reg_val *= detector.efficiency
  return reg_val
