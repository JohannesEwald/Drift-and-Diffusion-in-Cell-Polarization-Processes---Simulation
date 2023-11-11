import math
import numpy as np

# a, b, ce = 1, 1, 1
# f, g, h = 20, -20, 2 * math.log(10)          # +-+ rechtsdrift


def calculate(a, b, ce, f, g, h):

    amp = a * math.exp(+f / 2)
    apm = a * math.exp(-f / 2)
    bmp = b * math.exp(+g / 2)
    bpm = b * math.exp(-g / 2)
    cba = ce * math.exp(+h / 2)
    cab = ce * math.exp(-h / 2)


    p_A_plus = amp / (apm + amp)
    p_B_plus = bmp / (bpm + bmp)

    Z = [apm, amp, bpm, bmp, cab, cba]
    # print(Z)

    M = [[0, 0, 0, 0, 0],
         [0, -amp - bmp, apm + (1 - p_B_plus) * cba, bpm + (1 - p_A_plus) * cab, 0],
         [0, amp, -apm - cba - bmp, p_A_plus * cab, bpm],
         [0, bmp, p_B_plus * cba, -bpm - cab - amp, apm],
         [0, 0, bmp, amp, -bpm - apm]]
    M = np.array(M)

    A = 2 * math.sqrt(p_B_plus * cba * bpm * amp * p_A_plus * cab * apm * bmp)
    alpha = 1 / 2 * math.log((p_B_plus * cba * bpm * amp) / (p_A_plus * cab * apm * bmp))

    B = 2 * math.sqrt((1 - p_B_plus) * cba * amp * (1 - p_A_plus) * cab * bmp)
    beta = 1 / 2 * math.log(((1 - p_B_plus) * cba * amp) / ((1 - p_A_plus) * cab * bmp))

    E = math.sqrt(M[2, 2] * M[3, 3])
    epsilon = 1 / 2 * math.log(M[2, 2] / M[3, 3])

    F = math.sqrt(apm * amp * bpm * bmp)
    phi = 1 / 2 * math.log((apm * amp) / (bpm * bmp))

    # print(A, B, E, F, alpha, beta, epsilon, phi)

    J_d = - 1 / a_1(0, 0, A, B, E, F, alpha, beta, epsilon, phi, M) * ((M[1, 1] + M[4, 4]) * A * math.sinh(alpha) + M[4, 4] * B * E * math.sinh(
        beta - epsilon) + B * F * math.sinh(phi) * math.cosh(beta))
    J_c = - 1 / a_1(0, 0, A, B, E, F, alpha, beta, epsilon, phi, M) * ((M[1, 1] + M[4, 4]) * A * math.cosh(alpha) + M[4, 4] * B * E * math.cosh(
        beta - epsilon) + B * F * math.sinh(phi) * math.sinh(beta) + 2 * M[4, 4] * A * B / F * math.cosh(
        alpha - beta + phi) - 2 * M[1, 1] * M[4, 4] * ((A / F) ** 2))

    DD = - 1 / a_1(0, 0, A, B, E, F, alpha, beta, epsilon, phi, M) * ((M[1, 1] + M[4, 4]) * A * math.cosh(alpha) + M[4, 4] * B * E * math.cosh(
        beta - epsilon) + B * F * math.sinh(phi) * math.sinh(beta))

    return [J_d, J_c, DD]


# def M_tilde(c=0, d=0):
#     return [[-amp - bmp, apm + (1 - p_B_plus) * cba * math.exp(c + d), bpm + (1 - p_A_plus) * cab * math.exp(c - d), 0],
#             [amp, -apm - cba - bmp, p_A_plus * cab * math.exp(c - d), bpm],
#             [bmp, p_B_plus * cba * math.exp(c + d), -bpm - cab - amp, apm],
#             [0, bmp, amp, -bpm - apm]]


def a_2(c, d, A, B, E, F, alpha, beta, epsilon, phi, M):
    a_2_hat = E ** 2 + M[1, 1] * M[4, 4] - 2 * (M[1, 1] + M[4, 4]) * E * math.cosh(epsilon) - 2 * F * math.cosh(phi)
    x = a_2_hat - B * math.exp(c) * math.cosh(beta + d) - math.exp(2 * c) * (A ** 2) / (F ** 2)
    return x


def a_1(c, d, A, B, E, F, alpha, beta, epsilon, phi, M):
    a_1_hat = 2 * M[1, 1] * M[4, 4] * E * math.cosh(epsilon) - (M[1, 1] + M[4, 4]) * (
            E ** 2 - F * math.cosh(phi)) - 2 * E * F * math.cosh(epsilon) * math.cosh(phi)
    x = a_1_hat - 2 * A * math.exp(c) * math.cosh(alpha + d) + M[4, 4] * B * math.exp(c) * math.cosh(
        beta + d) - B * E * math.exp(c) * math.cosh(beta - epsilon + d) - math.exp(2 * c) * (
                    A * B / F * math.cosh(alpha - beta + phi) - (M[1, 1] + M[4, 4]) * (A ** 2) / (F ** 2))
    return x


def a_0(c, d, A, B, E, F, alpha, beta, epsilon, phi, M):
    a_0_hat = M[1, 1] * M[4, 4] * (E ** 2) + M[1, 1] * E * F * math.cosh(epsilon + phi) + M[4, 4] * E * F * math.cosh(epsilon - phi) + (F ** 2) * (math.sinh(phi) ** 2)
    x = a_0_hat + (M[1, 1] + M[4, 4]) * A * math.exp(c) * math.cosh(alpha + d) + M[4, 4] * B * E * math.exp(
        c) * math.cosh(beta - epsilon + d) + B * F * math.exp(c) * math.sinh(phi) * math.sinh(beta + d) + M[
            4, 4] * math.exp(2 * c) * (A * B / F * math.cosh(alpha - beta + phi) - M[1, 1] * (A ** 2) / (F ** 2))
    return x


# y = calculate(a, b, ce, f, g, h)
# print(y)

# CURRENTS
# print("CURRENTS:")
#
# J_d = - 1 / a_1(0, 0) * ((M[1, 1] + M[4, 4]) * A * math.sinh(alpha) + M[4, 4] * B * E * math.sinh(
#     beta - epsilon) + B * F * math.sinh(phi) * math.cosh(beta))
# J_c = - 1 / a_1(0, 0) * ((M[1, 1] + M[4, 4]) * A * math.cosh(alpha) + M[4, 4] * B * E * math.cosh(
#     beta - epsilon) + B * F * math.sinh(phi) * math.sinh(beta) + 2 * M[4, 4] * A * B / F * math.cosh(
#     alpha - beta + phi) - 2 * M[1, 1] * M[4, 4] * ((A / F) ** 2))
#
# print(J_d)
# print(J_c)
