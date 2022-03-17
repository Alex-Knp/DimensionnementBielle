from math import *
import numpy as np
from matplotlib import pyplot as plt

# moteur Opel Corsa C 1.3 CDTI (70 Hp)
tau = 18  # [-]
D = 0.0696  # [m]
C = 0.0854  # [m]
L = 0.0969  # [m]
mpiston = 0.600  # [kg] APPROXIMER
mbielle = 0.450  # [kg]
Q = 1650000  # [J/kg_inlet gas]


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
    V_output = [0] * 361
    Q_output = [0] * 361
    F_pied_output = [0] * 361
    F_tete_output = [0] * 361
    p_output = [0] * 361

    for t in range(len(theta)):
        V_output[t] = volume(t - 180)
        Q_output[t] = q_compute(t - 180, thetaC, deltaThetaC)

    p_output = p_theta(s, theta, deltaThetaC, thetaC)

    for t in range(len(theta)):
        F_pied_output[t] = f_pied(t - 180, p_output[t], rpm)
        F_tete_output[t] = f_tete(t - 180, p_output[t], rpm)

    max_pied = max(max(F_pied_output), abs(min(F_pied_output)))
    max_tete = max(max(F_tete_output), abs(min(F_tete_output)))
    peak_force = max(max_tete, max_pied)
    t = t_compute(peak_force)

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


def p_theta(s, theta, deltaThetaC, thetaC):
    # partie la plus dure du devoir -> il faut intégrer numériquement
    """
    calcule la pression dans le cylindre à l'angle moteur theta
    :param theta: angle moteur
    :param s: taux de suralimentation
    :return: pression dans le clindre
    """

    p = rungekutta(s*100000, theta, thetaC, deltaThetaC)

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
    c = ((L * 2) * 6) / ((pi * 2) * 131 * (10 ** 11))

    roots = np.roots([a, 0, b, 0, c])
    print(roots)

    for val in roots:
        if val > 0 and np.isreal([val]) == [True]:
            t = np.real(val)
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
    return (pi * (D * 2) / 4) * ptheta - mpiston * (C / 2) * ((6 * rpm) * 2) * cos(rad(theta))


def f_tete(theta, ptheta, rpm):
    """
    calcule la force exercée sur la tête de la bienne à l'angle moteur theta
    :param theta: angle moteur
    :param ptheta: pression dans le cylindre en theta
    :param rpm: vitesse moteur
    :return: force totale appliquée sur la tête de la bielle
    """
    return -(pi * (D * 2) / 4) * ptheta + (mpiston + mbielle) * (C / 2) * ((6 * rpm) * 2) * cos(rad(theta))


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
    if (theta < -thetaC): return 0

    if (theta > -thetaC + deltaThetaC): return Q * 10 ** -5 * 0.5 * (
            1 - cos(pi * (rad((deltaThetaC)) / rad(deltaThetaC))))

    return Q * 10 ** -5 * 0.5 * (1 - cos(pi * (rad((theta + thetaC)) / rad(deltaThetaC))))


def dvdt_compute(t):
    vc = (pi * (D ** 2) / 4) * C
    beta = 2 * L / C

    return (vc / 2) * (sin(rad(t)) + (sin(rad(t)) * cos(rad(t))) / (sqrt(beta ** 2 - (sin(rad(t))) ** 2)))


def dqdt_compute(t, thetaC, deltaThetaC):
    if (t < -thetaC or t > -thetaC + deltaThetaC): return 0

    return (pi * Q * 10 ** -2) / (2 * rad(deltaThetaC)) * sin((pi / rad(deltaThetaC)) * rad((t + thetaC)))


def rad(t):
    return 2 * pi * t / 360


def fun(p, theta, thetaC, deltaThetaC):
    return (-1.3 * p / volume(theta) * dvdt_compute(theta) + 0.3 * dqdt_compute(theta, thetaC, deltaThetaC) / volume(
        theta))


def rungekutta(p0, t, thetaC, deltaThetaC):
    P = [0] * len(t)
    P[0] = p0

    for j in range(len(t) - 1):
        K1 = fun(P[j], t[j], thetaC, deltaThetaC)
        K2 = fun(P[j] + K1 / 2, t[j] + 1 / 2, thetaC, deltaThetaC)
        K3 = fun(P[j] + K2 / 2, t[j] + 1 / 2, thetaC, deltaThetaC)
        K4 = fun(P[j] + K3, t[j] + 1, thetaC, deltaThetaC)
        P[j + 1] = P[j] + (K1 + 2 * K2 + 2 * K3 + K4) / 6

    return P


def euler_expl(r, theta, thetaC, deltaThetaC):
    P = [0] * len(theta)

    P[0] = r

    for k in range(len(theta) - 1):
        P[k + 1] = P[k] + fun(P[k], theta[k], thetaC, deltaThetaC)

    return P

#Tests
theta = [0] * 361
for i in range(361):
    theta[i] = i - 180

V, q, Fp, Ft, p, t = myfunc(3000, 1, theta, 5, 60)

dqdt = [0] * 361
dvdt = [0] * 361
dr = [0] * 361
ds = [0] * 361

for i in theta:
    dqdt[i - 180] = dqdt_compute(i, 5, 60)
    dvdt[i - 180] = dvdt_compute(i)
    dr[i - 180] = 0.3 * dqdt_compute(i, 5, 60) / volume(i) + (100000 * ((volume(-180) / volume(i) ** 1.3)))
    ds[i - 180] = (100000 * ((volume(-180) / volume(i) ** 1.3)))

plt.figure(figsize=(6, 3))

plt.plot(theta, p)
plt.axvline(x=-5, linestyle='-.', color='red', label="Début explosion", linewidth=0.5)
plt.axvline(x=80, linestyle='-.', color='green', label="Fin explosion", linewidth=0.5)

plt.show()

print("Q = ", q)
print("p = ", p)
print("t = ", t)
