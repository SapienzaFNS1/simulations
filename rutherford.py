# Proviamo a visualizzare alcune quantità di base calcolabili nella trattazione classica dell'esperimento di Geiger e Marsden. Importiamo innanzitutto i moduli necessari:

import numpy as np
import matplotlib.pyplot as plt

# Lavoriamo per semplicità nel sistema internazionale di unità di misura: definiamo quindi le costanti fondamentali,

e = 1.6e-19 # C
epsilon_0 = 8.85418782e-12 # SI units

# Fissiamo i parametri del problema:
# 1. la particella incidente è una particella alpha (nucleo di elio)
# 2. il nucleo incidente è un nucleo d'oro (numero atomico Z=79)
# 3. la particella incidente ha una energia cinetica T di 5 MeV

z = 2
Z = 79
T = 5.5e6*e # moltiplicato per "e" in modo da convertire eV in J

# Possiamo quindi calcolare la costante A del calcolo classico,

A = z * Z * np.power(e, 2) / (4 * np.pi * epsilon_0 * T)

# Consideriamo i diversi valori possibili dell'angolo di scattering theta, escludendo strategicamente il valore theta=0 (per cui sappiamo che la sezione d'urto diverge):

theta = np.linspace(1e-4, np.pi, 1000)

# Calcoliamo ora tre grandezze:
# 1. il parametro d'impatto b corrispondente a ogni dato valore di theta
# 2. la distanza minima fra la particella alpha e il nucleo (cioè il valore di r, rmin, nel punto di massimo avvicinamento fra alpha e nucleo)
# 3. la sezione d'urto differenziale dsigma/dOmega, che è funzione del solo angolo di scattering theta

def b(A, theta):
  return A/2 * 1/np.tan(theta/2)

def r(A, theta):
  _b = b(A=A, theta=theta)
  return A/2 * (1 + np.sqrt(1 + 4*np.power(_b/A, 2)))

def dSigma_dOmega(A, theta):
  return np.power(A, 2) / (16 * np.power(np.sin(theta/2), 4))

# Visualizziamo queste quantità. Partiamo dalla relazione fra parametro d'impatto e angolo di scattering:

def plot_b():
    plt.plot(theta*180/np.pi, b(A=A, theta=theta))
    plt.yscale('log')
    plt.ylabel('impact parameter [m]')
    plt.xlabel('scattering angle $\\theta$ (deg)')
    plt.savefig("impact_parameter.pdf")
    
    # Per valori elevati del parametro d'impatto, l'effetto del campo coulombiano del nucleo è pressoché nullo - per cui l'angolo di scattering è quasi zero.
    # Per valori sempre più piccoli del parametro d'impatto, invece, la particella alpha viene "rimbalzata all'indietro" (back-scattering).
    # La misura sperimentale di questi back-scattering è stata una delle prove eclatanti del fatto che la carica positiva dell'atomo è concentrata nel suo nucleo.
    
    
# Vediamo ora la distanza r:
def plot_r():
    plt.plot(theta*180/np.pi, r(A=A, theta=theta))
    plt.yscale('log')
    plt.ylabel('minimum distance from nucleus [m]')
    plt.xlabel('scattering angle $\\theta$ (deg)')
    plt.savefig('closest_approach.pdf')
    
    # Questo ci dice, ancora intuitivamente, che le particelle alpha che arrivano più vicine al nucleo atomico sono quelle che vengono rimbalzate all'indietro.

# La sezione d'urto differenziale è infine:
def plot_sigma():
    plt.plot(theta*180/np.pi, dSigma_dOmega(A=A, theta=theta))
    plt.yscale('log')
    plt.ylabel('${d\\sigma}/{d\\Omega}$ [m$^{2}$/sr]')
    plt.xlabel('scattering angle $\\theta$ (deg)')
    plt.savefig('sigma.pdf')
    
# Notiamo l'andamento divergente per angoli circa uguali a zero, dovuto al fatto che la forza di Coulomb agisce anche a distanza.


# Infine, proviamo a calcolare quante particelle alpha osserverebbe un rivelatore rettangolare di 2x1 cm posto a due metri di distanza, e a un angolo di 30 gradi, da un foglio d'oro di 1 mm, su cui incidono 1 milione di particelle alpha al secondo.

def calc_rutherford_detector():
    # I dati sono la superficie del rivelatore, la sua distanza (e quindi l'angolo solido deltaOmega coperto dal rivelatore),
    # il rate di particelle incidenti, e lo spessore e la densità del foglio d'oro:

    # posizione del rivelatore (determina la sezione d'urto di Rutherford)
    theta_detector = 30 * np.pi / 180
    print(f'La sezione d\'urto differenziale vale {dSigma_dOmega(A=A, theta=theta_detector):.3g} m2/sr')

    # angolo solido "visto" dal rivelatore
    S = 2e-2 * 1e-2 # m^2
    dist = 2 # m
    deltaOmega = S / np.power(dist, 2)
    print(f'Il rivelatore copre un angolo solido di {deltaOmega:.3g} sr')

    # numero di bersagli (= nuclei d'oro)
    rho = 19.320 # g/cm^3, da Wikipedia
    d = 1e-3 # m
    Navogadro = 6e23
    mass_no = 197 # numero di massa dell'oro, in g/mol, da Wikipedia
    N_nuclei = (rho*1e6) * Navogadro / mass_no * d # 1e6 per avere g/m^3 
    print(f'Le particelle alpha incontrano {N_nuclei:.3g} nuclei d\'oro per cm2 di sezione del fascio')
    
    
    # numero di particelle alpha incidenti nell'unità di tempo
    dNalpha_dt = 1e6

    dNrivelate_dt = deltaOmega * dSigma_dOmega(A=A, theta=theta_detector) * dNalpha_dt * N_nuclei
    print(f'Il rivelatore osserverà {dNrivelate_dt:.2f} particelle alpha al secondo')
    return (deltaOmega,dNalpha_dt,N_nuclei)
    
    ## Se spostiamo il rivelatore, che succede? Plottiamo il numero di eventi osservati nell'unità di tempo, in funzione dell'angolo a cui piazziamo il rivelatore:

def plot_rutherford_detector():
    theta = np.linspace(1e-4, 2*np.pi, 1000)
    #theta = 30*np.pi/180
    deltaOmega,dNalpha_dt,N_nuclei = calc_rutherford_detector()
    rate = deltaOmega * dSigma_dOmega(A=A, theta=theta) * dNalpha_dt * N_nuclei # un vettore di 1000 elementi, con i 1000 rate per ciascun valore dell'angolo theta

    plt.plot(theta*180/np.pi, rate)
    plt.ylim(1e-2, 1e9) # tagliamo l'asse y altrimenti non si vede niente
    plt.yscale('log') # come sopra
    plt.ylabel('rate of detected alpha particles [Hz]')
    plt.xlabel('scattering angle [deg]')
    plt.savefig("rutherford_scattering_angle.pdf")

if __name__ == '__main__':

    # Grafico 1. il parametro d'impatto b corrispondente a ogni dato valore di theta
    # plot_b()

    # Grafico 2. la distanza minima fra la particella alpha e il nucleo (cioè il valore di r, rmin, nel punto di massimo avvicinamento fra alpha e nucleo)
    # plot_r()
    
    # Grafico 3. la sezione d'urto differenziale dsigma/dOmega, che è funzione del solo angolo di scattering theta
    # plot_sigma()

    # Stampiamo i dati intermedi e il numero di particelle viste dal rivelatore dell'esperimento di Rutherford-Geiger e Marsden
    # con numeri realistici
    # calc_rutherford_detector()

    # grafichiamo il n. di eventi osservati nell'unita' di tempo:
    plot_rutherford_detector()
