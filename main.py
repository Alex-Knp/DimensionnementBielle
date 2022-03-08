from math import *
import numpy as np

"""
tau = @valeur taux compression@ #[-]
D = @valeur alesage@ #[m]
C = @valeur course@ #[m]
L = @valeur longueur bielle@ #[m]
mpiston = @valeur masse piston@ #[kg]
mbielle = @valeur masse bielle@ #[kg]
Q = @valeur chaleur emise par fuel par kg de melange admis@ #[J/kg_inlet gas]
"""


def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    """
    fonction main
    :param rpm: vitesse moteur
    :type rpm: int
    :param s: pression turbo
    :type s: float
    :param theta: angles moteur
    :type theta: list[int]
    :param thetaC: angle d'allumage
    :type thetaC: int
    :param deltaThetaC: durée de combustion
    :type deltaThetaC: int
    :return:
    """
    V_output = []
    Q_output = []
    dVdt = []
    dQdt = []
    F_pied_output = []
    F_tete_output = []
    p_output = []

    for t in range(len(theta)):
        V_output[t] = volume(t)
        Q_output[t] = q_compute(t, thetaC, deltaThetaC)
        dVdt[t] = dvdt_compute(t)
        dQdt[t] = dqdt_compute(t, thetaC, deltaThetaC)

    p_output = p_theta(s, theta, dVdt, V_output, dQdt,thetaC)

    for t in range(len(theta)):
        F_pied_output[t] = f_pied(t, p_output[t], rpm)
        F_tete_output[t] = f_tete(t, p_output[t], rpm)

    max_pied = max(max(F_pied_output), abs(min(F_pied_output)))
    max_tete = max(max(F_tete_output), abs(min(F_tete_output)))
    peak_force = max(max_tete, max_pied)
    t = t_compute(peak_force)

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


def p_theta(s, theta, dVdt, V_output, dQdt, thetaC):
    # partie la plus dure du devoir -> il faut intégrer numériquement
    """
    calcule la pression dans le cylindre à l'angle moteur theta
    :param theta: angle moteur
    :param s: taux de suralimentation
    :return: pression dans le clindre
    """

    gamma = 1.3
    p = [0]*361
    for angle in theta:
        if (angle <= thetaC):
            p[angle + 180]= (s * 100000 * ((volume(-180) / volume(angle)) ** gamma))

    rungekutta(360-thetaC, p[thetaC])

    return p


def t_compute(peak_force):
    """
    calcule la dimension de la bielle en fonction de la force maximale à
    laquelle elle va être soumise
    :param peak_force: force maximale exercée sur la bielle
    :type peak_force : int
    :rtype: int
    :return: dimension t de la section en I
    """
    # il faut résoudre un 4e degré de la forme a*t^4 + b*t^2 + c = 0
    t = None

    a = -1 / peak_force
    b = 1 / (495 * 10 ** 7)
    c = ((L ** 2) * 6) / ((pi ** 2) * 131 * (10 ** 11))

    roots = np.roots([a, 0, b, 0, c])

    for val in roots:
        if val > 0 and np.isreal([val]) == [True]:
            t = val
            break

    return t


def f_pied(theta, ptheta, rpm):
    """
    calcule la focre exercée sur le pied de la bielle à l'angle moteur theta
    :param theta: angle moteur
    :param ptheta: pression dans le cylindre en theta
    :param rpm: vitesse moteur
    :return: force totale appliquée sur le pied de la bielle
    """
    return (pi * (D ** 2) / 4) * ptheta - mpiston * (C / 2) * ((6 * rpm) ** 2) * cos(rad(theta))


def f_tete(theta, ptheta, rpm):
    """
    calcule la force exercée sur la tête de la bienne à l'angle moteur theta
    :param theta: angle moteur
    :param ptheta: pression dans le cylindre en theta
    :param rpm: vitesse moteur
    :return: force totale appliquée sur la tête de la bielle
    """
    return -(pi * (D ** 2) / 4) * ptheta + (mpiston + mbielle) * (C / 2) * ((6 * rpm) ** 2) * cos(rad(theta))


def volume(theta):
    """
    calcule le volume du cylindre à l'angle moteur theta
    :param theta: angle moteur
    :return: volume du cylindre
    """
    vc = (pi * (D ** 2) / 4) * C
    beta = 2 * L / C

    return (vc / 2) * (1 - cos(rad(theta)) + beta - sqrt(beta ** 2 - sin(rad(theta)) ** 2)) + vc / (tau - 1)


def q_compute(theta, thetaC, deltaThetaC):
    """/!\ vérifier que Q est bien la variable qu'il faut"""

    return Q * 0.5 * (1 - cos(rad(pi * ((theta - thetaC) / deltaThetaC))))


def dvdt_compute(t):
    vc = (pi * (D ** 2) / 4) * C
    beta = 2 * L / C

    return (vc/2)*(sin(rad(t))+(sin(rad(t))*cos(t))/(sqrt(beta**2-(sin(rad(t)))**2)))


def dqdt_compute(t, thetaC, deltaThetaC):


    return (pi*Q)/(2*deltaThetaC)*sin(rad((pi/deltaThetaC)*(t-thetaC)))



def rad(t):
    return 360 * t / (2 * pi)
