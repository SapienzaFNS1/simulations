import numpy as np
import matplotlib.pyplot as plt

#Definiamo intanto una funzione che ci permette di conoscere l'energia finale di un fotone che diffonde su una particella massiva (effetto Compton), noti la sua energia iniziale e l'angolo a cui viene diffuso:

def get_Ef(Ei, theta, M):
  return Ei / (1 + Ei/M * (1 - np.cos(theta)))


#Generiamo un certo numero di fotoni con una delta di Dirac come  distribuzione energetica, centrata sul valore $E_0$:

def generate_energy_fixed(E0=500,Nphotons=100000): # energy in keV
  Ei = np.ones(Nphotons) * E0 # sarà un array con 10'000 entries, tutte uguali a E0
  return Ei

def generate_energy_gaus(E0=500,sigma=50,Nphotons=100000): # energy in keV
  Ei = np.random.normal(E0,sigma,Nphotons)
  return Ei

def generate_energy_uniform(Emin=400,Emax=600,Nphotons=100000): # energy in keV
  Ei = np.random.uniform(Emin,Emax,Nphotons)
  return Ei

# distribuizione angolare degli angoli secondo Klein-Nishima non relativistico
# data e' l'array degli angoli theta
# Dal nostro punto di vista, si tratta di definire una funzione che inverta la cumulativa di una distribuzione. Nella fattispecie, per noi la distribuzione sarà l'istogramma
# di una certa quantità (theta, che è un array - l'argomento data della definizione della funzione qua sotto), interpolato opportunamente (vedi la chiamata a interpolate.interp1d).
import scipy.interpolate as interpolate # ci serve per interp1d
def inverse_transform_sampling(data, weights, n_bins=40, n_samples=1000):
    hist, bin_edges = np.histogram(data, weights=weights, bins=n_bins, density=True)
    cum_values = np.zeros(bin_edges.shape)
    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
    r = np.random.rand(n_samples)
    return inv_cdf(r)

def getTheta_KleinNishima(Nphotons):
  thetas = np.linspace(0, 2*np.pi, 10000) # x della funzione di distribuzione: da 0 a pi a passi di 1/10000
  thetadist = 1 + np.power(np.cos(thetas), 2) # y, cioè il valore della funzione di distribuzione
  theta_1pluscos2 = inverse_transform_sampling(thetas, thetadist, n_bins=len(thetas), n_samples=Nphotons)
  theta_unif = np.random.uniform(0, 2*np.pi, Nphotons)
  
  # Plottiamo l'istogramma dell'array  𝜃 , per assicurarci che tutto sia andato a buon fine:
  plt.hist(theta_1pluscos2, bins=100, label='$\\frac{dN}{d\\theta}\\propto 1+\\cos^2{\\theta}$')
  plt.hist(theta_unif, bins=100, label='$\\frac{dN}{d\\theta}\\propto \\mathrm{cost}$', alpha=0.3)
  plt.xlabel('$\\theta$')
  plt.ylabel('# of photons')
  plt.legend()
  plt.savefig("comp_angles.pdf")
  return theta_1pluscos2 


def plotOne(Ef,theta):
    fig, ax = plt.subplots(1, 2, figsize=(15,5))

    ax[0].hist(theta, bins=100)
    ax[0].set_xlabel('$\\theta$')
    ax[0].set_ylabel('# of photons')

    ax[1].hist(Ef, bins=100)
    ax[1].set_xlabel('$E_{f}$ [keV]')
    ax[1].set_ylabel('# of photons')

    fig.savefig("compton.pdf")

def plotTwo(Ef1,Ef2):
  fig, ax = plt.subplots(1, 1, figsize=(10,10))

  ax.hist(Ef_1pluscos2, bins=100, label='$\\frac{dN}{d\\theta}\\propto 1+\\cos^2{\\theta}$')
  ax.hist(Ef_unif, bins=100, label='$\\frac{dN}{d\\theta}\\propto \\mathrm{cost}$', alpha=0.3)
  ax.set_xlabel('$E_{i}$ [keV]')
  ax.set_ylabel('# of photons')
  ax.legend()

  fig.savefig("compton_comp.pdf")
  
  
if __name__ == '__main__':

    Ei = 500 #keV
    Nphotons = 100000
    
    theta_unif = np.random.uniform(0, 2*np.pi, Nphotons)
    theta_kn   = getTheta_KleinNishima(Nphotons)
    theta = theta_unif
    
    M = 511 # keV
    #Ei = generate_energy_fixed(Ei,Nphotons)   # Delta di Dirac at Ei
    Ei = generate_energy_gaus(Ei,0.05*Ei,Nphotons) # Gaussian centered at Ei=500 keV and resolution=5%
    #Ei = generate_energy_uniform(Ei-0.2*Ei,Ei+0.2*Ei,Nphotons) # Gaussian centered at Ei=500 keV and resolution=5%
    #Ef = get_Ef(Ei, theta, M=M)
    
    #plotOne(Ef,theta)

    ## confrontiamo le distribuzioni di energia nei due casi
    Ef_unif = get_Ef(Ei, theta_unif, M=M)
    Ef_1pluscos2 = get_Ef(Ei, theta_kn, M=M)
    plotTwo(Ef_unif,Ef_1pluscos2)
    
    
