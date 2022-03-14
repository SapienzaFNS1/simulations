import numpy as np # ci permette di lavorare su array
import matplotlib.pyplot as plt # ci permette di fare i plot


# Generiamo un vettore contenente gli angoli possibili di una data particella nel riferimento del centro di massa:
theta_star = np.linspace(0, np.pi, 10000) # theta_star sarà un array di 10'000 elementi equispaziati fra 0 e pi

# theta_star, cioè l'angolo nel sistema del centro di massa fra la direzione di moto di una determinata particella e l'asse coincidente con il boost di Lorentz (cioè l'asse parallelo a quello lungo il quale, nel riferimento del laboratorio, si muove il centro di massa);

# betaCOM = pTot(LAB)/ETot(LAB), cioè il boost di Lorentz del centro di massa rispetto al laboratorio;

# beta_star=p*/E*, cioè il boost di Lorentz della particella considerata calcolato nel riferimento del centro di massa Sottolineiamo che queste ultime due quantità sono in generale diverse: nel caso generale, il sistema del centro di massa è composto da 𝑁 particelle, ma noi ci chiediamo come cambiano le coordinate di una di esse nel passaggio da un sistema di riferimento all'altro

def tantheta(theta_star, betaCOM, betaParticleInCOM):
  gammaCOM = 1 / np.sqrt(1 - pow(betaCOM, 2))
  return np.sin(theta_star) / (gammaCOM * (betaCOM/betaParticleInCOM + np.cos(theta_star)))

# Definiamo poi una funzione, più comoda, che ci ritorni proprio l'angolo (cioè l'arcotangente della tangente). Dovremo prendere il risultato dell'operazione "arcotangente", che è espresso fra  −𝜋/2  e  𝜋/2 , e riportarlo fra  0  e  𝜋 .
def get_theta(theta_star, betaCOM, betaParticleInCOM):
  # atan(y/x) è definita tra -pi/2 e pi/2, atan2(y, x) invece è definita su tutto il piano
  gammaCOM = 1 / np.sqrt(1 - pow(betaCOM, 2))
  return np.arctan2(np.sin(theta_star), (gammaCOM * (betaCOM/betaParticleInCOM + np.cos(theta_star))))

def plotThetas(betaCOM_values,betaParticleInCOM_values, theta_star):
    fig = plt.figure(figsize=(9, 9), dpi=80, facecolor='w', edgecolor='k') # definiamo le dimensioni e l'aspetto della figura
    for betaCOM in betaCOM_values: # loop sui valori di betaCOM considerati
        for betaParticleInCOM in betaParticleInCOM_values: # loop sui valori di betaParticleInCOM considerati
            # aggiungiamo la serie di punti (x, y) dove x sono i 10'000 valori di theta_star considerati, e y i valori theta corrispondenti
            items = plt.plot(theta_star, get_theta(theta_star, betaCOM, betaParticleInCOM), label=f'$\\beta_{{COM}}={betaCOM:.2f}$, $\\beta^*={betaParticleInCOM:.2f}$')
            if betaCOM <= betaParticleInCOM:
                # dove ce lo aspettiamo, aggiungiamo anche la linea verticale corrispondente al valore per cui la particella è a pi/2 nel laboratorio
                plt.axvline(np.arccos(-betaCOM/betaParticleInCOM), label='$\\theta^*=\\arccos\\left(-\\beta_{COM}/\\beta^*\\right)$', linestyle='dotted', color=items[-1].get_color())
    plt.xlabel('$\\theta^*$')
    plt.ylabel('$\\theta$')
    # semplifichiamo la lettura del plot includendo una riga orizzontale a pi/2
    plt.axhline(np.pi/2, label='$\\theta=\\pi/2$', linestyle='dashed', color='black')
    plt.legend()
    #plt.show()
    plt.savefig("thetas_lab.pdf")



##  Angolo fra due particelle emesse in un decadimento a due corpi ##

# Immaginiamo di avere una particella A che decade in una particella B e in una particella C, attraverso il processo
# 𝐴→𝐵+𝐶. 
# Ragioniamo innanzitutto nel sistema di riferimento del centro di massa. Per la conservazione dell'impulso spaziale, le due particelle saranno emesse ad angoli  𝜃∗𝐵  e  𝜃∗𝐶=𝜋−𝜃∗𝐵 .
# Ipotizziamo che tutti gli angoli siano equiprobabili (vedi l'esempio sull'effetto Compton per un caso generale), e scriviamo quindi due array che rappresentino gli angoli:

def openingAngle(betaCOM_values,betaParticleInCOM_values):
    thetaB_star = np.linspace(0, np.pi, 10000) # theta_star sarà un array di 10'000 elementi equispaziati fra 0 e pi
    thetaC_star = np.pi - thetaB_star # back-to-back

    fig = plt.figure(figsize=(9, 9), dpi=80, facecolor='w', edgecolor='k') # definiamo le dimensioni e l'aspetto della figura
    for betaCOM in betaCOM_values: # loop sui valori di betaCOM considerati
        for betaParticleInCOM in betaParticleInCOM_values: # loop sui valori di betaParticleInCOM considerati
            # calcoliamo thetaB e thetaC
            thetaB = get_theta(thetaB_star, betaCOM, betaParticleInCOM)
            thetaC = get_theta(thetaC_star, betaCOM, betaParticleInCOM)
            # calcoliamo deltatheta
            delta_theta = np.arccos(np.cos(thetaB - thetaC))
            # grafichiamolo in funzione di thetaB_star
            items = plt.plot(thetaB_star, delta_theta, label=f'$\\beta_{{COM}}={betaCOM:.2f}$, $\\beta^*={betaParticleInCOM:.2f}$')
    
    plt.xlabel('$\\theta_B^*$')
    plt.ylabel('$\\Delta\\theta$')
    plt.axvline(np.pi/2., label='$\\theta^*=\\pi/2$', linestyle='dashed', color='green')
    #plt.ylim(-10,10)
    plt.legend()
    plt.savefig("deltaTheta.pdf")



    
if __name__ == '__main__':

    # Definiamo alcuni casi specifici di  betaCOM e beta_star :
    betaCOM_values = [0.6]
    betaParticleInCOM_values = [0.3, 0.59, 0.6, 0.61, 0.7, 0.9]

    # Plottiamo i valori di  𝜃  che otterremmo scansionando tutti i valori possibili di  𝜃∗  - ovvero, plottiamo la funzione  𝜃=𝑓(𝜃∗) . Aggiungiamo al plot, per semplificare l'interpretazione dei risultati:
    # - una riga verticale, nei casi in cui  𝛽𝐶𝑂𝑀≤𝛽∗ , per rappresentare il punto in cui  tan𝜃  va a infinito;
    # - una riga orizzontale a  𝜃=𝜋/2 , per ricordarci da dov'è che si inizia a parlare di particella emessa all'indietro.
    plotThetas(betaCOM_values,betaParticleInCOM_values,theta_star)
    
    # Calcoliamo ora gli angoli nel sistema di riferimento del laboratorio e grafichiamo la loro differenza. Facciamo ciò in tutti gli "scenari" (diversi valori di  𝛽∗=𝛽∗𝐵=𝛽∗𝐶  e  𝛽𝐶𝑂𝑀) :
    # betaCOM_values = [0.60]
    # betaParticleInCOM_values = [0.10,0.25,0.50,0.60,0.90,0.99]
    # openingAngle(betaCOM_values,betaParticleInCOM_values)
    




    
    
