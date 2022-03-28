import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
norm.cdf(1.96)

# troviamo l'area sottesa dalle due PDF Gaussiane. Siccome sono probabilita' (normalizzate a 1), l'area e' la frazione di paticelle che si incontrano in ciascuna direzione (x e y)
# nella collisione dei fasci
def solve(m1,m2,std1,std2):
  a = 1/(2*std1**2) - 1/(2*std2**2)
  b = m2/(std2**2) - m1/(std1**2)
  c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
  return np.roots([a,b,c])

def calc_overlap(m1,std1,m2,std2,plot=False):
    
    #m1 = 2.5
    #std1 = 1.0
    #m2 = 5.0
    #std2 = 1.0
    
    #Get point of intersect
    result = solve(m1,m2,std1,std2)

    if plot:
        #Get point on surface
        #x = np.linspace(-5,9,10000)
        x = np.linspace(m1-5*std1,m2+5*std2,10000)
        plot1=plt.plot(x,norm.pdf(x,m1,std1),label="electron beam")
        plot2=plt.plot(x,norm.pdf(x,m2,std2),label="positron beam")
        plot3=plt.plot(result,norm.pdf(result,m1,std1),'o')
        plt.xlabel('beam profile along one direction [m]')
        plt.ylabel('probability density function')
        plt.legend()
        
    r = result[0]
    if plot:
        #Plots integrated area
        olap = plt.fill_between(x[x>r], 0, norm.pdf(x[x>r],m1,std1),alpha=0.3)
        olap = plt.fill_between(x[x<r], 0, norm.pdf(x[x<r],m2,std2),alpha=0.3)
    
    # integrate
    area = norm.cdf(r,m2,std2) + (1.-norm.cdf(r,m1,std1))

    if plot:
        print("Area under curves ", area)
        #plt.show()
        plt.savefig("simple_overlap.pdf")

    return area


# luminosita' del collider circolare e+ e-, dato il numero di pacchetti nbunch, il numero di elettroni Ne, il numero di positroni Np, il raggio R del collider.
# A = superficie di base del cilindro con cui approssimiamo il pacchetto
def luminosity(nbunch, Ne, Np, R, A):

    c = 3e8 # m/s
    freq = nbunch * c / (2*np.pi*R) # fequenza di collisione = nbunch / periodo
    inst_lumi = Ne * Np / A * freq
    print(f'la freq di collisione di {freq:.2E} s-1')
    print(f'La luminosita` istantanea di {inst_lumi:.2E} m-2s-1')
    return inst_lumi

# grafichiamo la luminosita' in funzione della distanza tra i fasci
def plot_lumi(nbunch, Ne, Np, R, sigma):

    # un vettore di 1000 elementi, con i 1000 rate per ciascun valore della distanza tra i fasci
    dist = np.linspace(1e-10, 2e-5, 1000)

    lumi = []
    for d in dist:
        overlap = np.power(calc_overlap(0,sigma,d,sigma),2)
        lumi.append(luminosity(nbunch, Ne*overlap, Np*overlap, R, A))
    lumis = np.array(lumi)
    
    plt.plot(dist, lumis)
    plt.xlabel("distance between beam axes [m]")
    plt.ylabel("instantaneous luminosity [$m^{-2}s^{-1}$]")
    plt.savefig("lumi_vs_dist.pdf")
    
if __name__ == '__main__':

    # caso collisione tra cilindretti
    A = 0.1e-6 ## superficie di base del cilindretto 0.1 mm^2
    Ne = 1e12 # numero di elettroni per pacchetto
    Np = 1e9 # numero di positroni per pacchetto
    R = 2e3 # 2 km
    nbunch = 4 # numero di pacchetti di ciascun tipo
    
    # controlliamo il valore della liminosita' istantanea
    # luminosity(nbunch, Ne, Np, R, A)

    sigma=3.5e-6 # 3.5 microns per ogni pacchetto (uguali dimensioni)
    # mettiamo la posizione di un fascio a 0, poiche' conta solo la distanza: 1 micron tra i due
    m1=0
    m2=7e-6
    
    # semplice controllo che il calcolo dell'overlap sia giusto
    #calc_overlap(2.5,1,5,1)
    
    # approssimazione rozza: meglio sarebbe calcolare l'overlap di 2 gaussiane 2D
    overlap = np.power(calc_overlap(m1,sigma,m2,sigma,plot=True),2)
    print(f'Frazione di overlap e` di {overlap:.3f}')

    #plot_lumi(nbunch,Ne,Np,R,sigma)
    
