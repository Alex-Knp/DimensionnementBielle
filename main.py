from math import *
import numpy as np

tau = @valeur taux compression@ #[-]
D = @valeur alesage@ #[m]
C = @valeur course@ #[m]
L = @valeur longueur bielle@ #[m]
mpiston = @valeur masse piston@ #[kg]
mbielle = @valeur masse bielle@ #[kg]
Q = @valeur chaleur emise par fuel par kg de melange admis@ #[J/kg_inlet gas]


def myfunc(rpm, s, theta, thetaC, deltaThetaC):


    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


def t_compute(peak_force):
    """
    :param peak_force: force maximale exercée sur la bielle
    :type peak_force : int
    :rtype: int
    :return: dimension t de la section en I
    """
    #il faut résoudre un 4e degré de la forme a*t^4 + b*t^2 + c = 0
    t = None

    a = -1/peak_force
    b = 1/(495*10**7)
    c = ((L**2)*6)/((pi**2)*131*(10**11))

    roots = np.roots([a,0,b,0,c])

    for val in roots:
        if(val > 0 and np.isreal([val]) == [True]):
            t = val
            break

    return t

def f_pied(theta, p_theta, rpm):
    """
    :param theta: angle moteur
    :param p_theta: pression dans le cylindre en theta
    :param rpm: vitesse moteur
    :return: force totale appliquée sur le pied de la bielle
    """
    return (pi*(D**2)/4)*p_theta - mpiston*(C/2)*((6*rpm)**2)*cos(theta)



def f_tete(theta, p_theta, rpm):
    """
    :param theta: angle moteur
    :param p_theta: pression dans le cylindre en theta
    :param rpm: vitesse moteur
    :return: force totale appliquée sur la tête de la bielle
    """
    return -(pi*(D**2)/4)*p_theta + (mpiston+mbielle)*(C/2)*((6*rpm)**2)*cos(theta)



def volume(theta):
    """
    :param theta: angle moteur
    :return: volume du cylindre
    """
    vc = (pi*(D**2)/4)*C
    beta = 2*L/C

    return (vc/2)*(1-cos(theta)+beta-sqrt(beta**2-sin(theta)**2))+vc/(tau-1)

def q_compute(theta , thetaC, deltaThetaC):
    """/!\ vérifier que Q est bien la variable qu'il faut"""

    return Q*0.5*(1-cos(pi*((theta-thetaC)/deltaThetaC)))