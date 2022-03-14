import numpy as np # ci permette di lavorare su array
import matplotlib.pyplot as plt # ci permette di fare i plot


def higgspt(NHiggs=100000):
    ## approximate the Higgs boson pT distribution with a Landau with mu=10 GeV and width=15 GeV
    mu = 15 # GeV
    width = 20 # GeV
    from scipy.stats import moyal
    ## Landau distribution
    pTH = moyal.rvs(mu,width,NHiggs)

    # Plottiamo l'istogramma dell'array ptH, per assicurarci che tutto sia andato a buon fine:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(pTH, bins=100, label='$\\frac{dN}{dp_{T}^{Higgs}}$')
    plt.xlabel('$p_{T}^{Higgs}$')
    plt.ylabel('# of Higgs bosons')
    plt.legend()
    plt.savefig("higgspt.pdf")
    return pTH

# distribuzione tipica dei pi0 nei jet a LHC
def pi0pt(N=100000):
    ## approximate the Higgs boson pT distribution with a Landau with mu=10 GeV and width=15 GeV
    mu = 2 # GeV
    width = 5 # GeV
    from scipy.stats import moyal
    pTpi0 = moyal.rvs(mu,width,N)
    return pTpi0
    
# relazione dell'angolo minimo di apertura per una particella di massa M e energia E
# nota bene: qui l'energia e' sostituita dall'energia trasversa perche' a LHC l'impulso longitudinale e' ignoto, perche' non si conosce il momento longitudinale dei quark all'interno del protone
def deltaTheta_min(pts,M,particle):

    Ets = Et = np.sqrt(M*M + np.power(pts,2)) 
    delta_thetas = 180./np.pi * 2 * np.arcsin(M/Ets) # in gradi

    ylabel = "Higgs bosons" if particle == "higgs" else "$\\pi^{{0}} mesons$"
    xlabel = "Higgs" if particle == "higgs" else "\\pi^{{0}}"

    ## faccio il grafico della distribuzione dei pT della particella che decade
    fig, ax = plt.subplots(1, 2, figsize=(15,5))
    ax[0].hist(pts, bins=100)
    ax[0].set_xlabel('$p_{{T}}^{{{part}}}$'.format(part=xlabel))
    ax[0].set_ylabel('# of {part}'.format(part=ylabel))

    ## faccio il grafico dell'angolo MINIMO di apertura 
    ax[1].hist(delta_thetas, bins=100)
    ax[1].set_xlabel('$\\Delta\\theta_{\\gamma\\gamma} [^{\\circ}]$')
    ax[1].set_ylabel('# of {part}'.format(part=ylabel))

    ## linea che rappresenta un esempio realistico di dimensione angolare di una cella calorimetrica
    plt.axvline(5, label='typical $\\Delta\\theta$ calorimeter cell', linestyle='dotted', color='green')
    plt.legend()
    fig.savefig("openangle_{p}.pdf".format(p=particle))

    
if __name__ == '__main__':

  ### Higgs boson -> 2 photons
  #pts = higgspt()
  #MH = 125 # GeV
  #deltaTheta_min(pts,MH,"higgs")
  
  ### pi0 -> 2 photons
  pts = pi0pt()
  Mpi0 = 0.135 # GeV
  deltaTheta_min(pts,Mpi0,"pi0")


    
    
