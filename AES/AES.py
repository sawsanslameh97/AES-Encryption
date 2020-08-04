import numpy as np
from collections import deque
from typing import Deque

base_coeff = [1, 0, 0, 0, 1, 1, 0, 1, 1]  # x**8 + x**4 + x**3 +x+ 1
base = np.poly1d(base_coeff)


def Text_toArray(text):
    arrFromtext = text.split(",")
    arrFromtext = np.array(arrFromtext)
    state = np.reshape(arrFromtext, (-1, 4))
    twoDimArray = state.transpose()

    return twoDimArray


def strTopoly(list):
    base_array_int = [[int(num, 16) for num in row] for row in list]  # int

    base_array_bin = [[int(bin(num)[2:]) for num in row] for row in base_array_int]  # bin

    base_array_bin_arr = [[[int(x) for x in str(num)] for num in row] for row in base_array_bin]  # arr

    base_array_poly = [[np.poly1d(num) for num in row] for row in base_array_bin_arr]

    return base_array_poly


def polyToHex(poly):
    int_coeff = poly.coeffs
    list_int_coeff = [int(coef) for coef in int_coeff]

    # convert to binary num
    int_num = int("".join(str(x) for x in list_int_coeff), 2)

    hexnum = hex(int_num)[2:]
    return hexnum


def strTobin(str):
    int_num = int(str, 16)
    bin_num = bin(int_num)[2:]
    return bin_num


def strNumTopoly(str):
    int_num = int(str)
    bin_num = bin(int_num)[2:]

    eq_coeff = [int(i) for i in bin_num]
    equ = np.poly1d(eq_coeff)

    remainder = modArray(equ)
    return remainder


def modArray(equ):
    _, remainder = np.polydiv(equ, base)
    quo_coeff_after_mod = divmod(remainder.coeffs, 2)[1]
    remainder = np.poly1d(quo_coeff_after_mod)
    return remainder


def xor(state, key):
    states = ""
    for row in range(len(state)):
        for col in range(len(state)):
            states += hex(int(state[row][col], 16) ^ int(key[row][col], 16))[2:] + ","

    states = Text_toArray(states[:-1]).transpose()
    return states


def inverse():
    equ_coeff = [1, 0, 1, 1, 0, 1, 0, 1]  # X**7 + x**5+x**4+x**2+1
    equ = np.poly1d(equ_coeff)

    print("the equation is ")
    print(equ)
    if (equ.order < base.order):
        inv = findInverse(equ, base)
        print("---------------------")
        print("The inv :")
        print(inv)
    else:
        print("you have error")


def findInverse(equ, base=base):
    mono_coeff = [1]
    mono_equ = np.poly1d(mono_coeff)  # mono equation

    rounds = 0
    a0 = a00 = 1
    a1 = a11 = 0
    a2 = a22 = base

    b0 = 0
    b1 = 1
    b2 = equ

    if equ == mono_equ:
        return equ

    while True:
        rounds = rounds + 1
        quotient, remainder = np.polydiv(a2, b2)

        # convert the coeff to 0 , 1

        quo_coeff_after_mod = divmod(quotient.coeffs, 2)[1]
        rem_coeff_after_mod = divmod(remainder.coeffs, 2)[1]

        quotient = np.poly1d(quo_coeff_after_mod)
        remainder = np.poly1d(rem_coeff_after_mod)

        # convert the array of coeff to array of int
        coeff_int_list = [int(coef) for coef in b2.coef]

        if coeff_int_list == [0]:
            return "no inverse"

        elif (coeff_int_list == [1]):
            return b1

        else:
            a00 = b0
            a11 = b1
            a22 = b2

            b0 = a0 - (quotient * b0)
            b1 = a1 - (quotient * b1)
            b2 = remainder

            a0 = a00
            a1 = a11
            a2 = a22

            b1_coeff_after_mod = divmod(b1.coeffs, 2)[1]

            b1 = np.poly1d(b1_coeff_after_mod)

            continue


def sBox(num):
    num63 = '01100011'  # 63 in binary

    if num == "00":
        num = format(int(num), '#010b')[2:]  # make it 8bit num
        return hex(int(num, 2) ^ int(num63, 2))[2:]

    # convert str to hex to bin
    bin_num = bin(int(num, 16))[2:]

    # convert hex num to list of coeff
    eq_coeff = [int(i) for i in str(bin_num)]
    equ = np.poly1d(eq_coeff)

    # find inverse
    inv = findInverse(equ, base)
    inv_poly_coeff = inv.coeffs
    list_int_coeff = [int(coef) for coef in inv_poly_coeff]

    # convert to binary num
    int_num = int("".join(str(x) for x in list_int_coeff), 2)
    binary_num = format(int_num, '#010b')[2:]  # make it 8bit num

    # start Xor the num with 3
    binary_after_xor = ""
    for i in range(0, 8):
        firstbit = i
        secondbit = (i - 4) % 8
        thirdbit = (i - 5) % 8
        fourthbit = (i - 6) % 8
        fifthbit = (i - 7) % 8

        result_xor = int(binary_num[firstbit], 2) ^ int(binary_num[secondbit], 2) ^ int(binary_num[thirdbit], 2) ^ int(
            binary_num[fourthbit], 2) ^ int(binary_num[fifthbit], 2) ^ int(num63[i], 2)

        binary_after_xor = binary_after_xor + str(result_xor)

    hexnum = hex(int(binary_after_xor, 2))[2:]
    return hexnum


def mixcol(state):
    base_array = [["02", "03", "01", "01"], ["01", "02", "03", "01"], ["01", "01", "02", "03"],
                  ["03", "01", "01", "02"]]

    base_array = strTopoly(base_array)
    state = strTopoly(state)

    # mul arr
    mul_array = np.tensordot(base_array, state, axes=(1, 0))

    # mod for element
    for row in range(0, len(mul_array)):
        for item in range(0, len(mul_array)):
            mul_array[row][item] = modArray(mul_array[row][item])  # mod
            mul_array[row][item] = polyToHex(mul_array[row][item])  # hex

    return mul_array


def generateKey(key, rounds):
    rc = ["1", "2", "4", "8", "10", "20", "40", "80", "1b"]

    key = key.transpose()
    w0 = key[0]
    w1 = key[1]
    w2 = key[2]
    w3 = key[3]

    w3 = deque(w3)
    w3.rotate(-1)
    w3 = list(w3)

    g = []
    for num in w3:
        g.append(sBox(num))

    g[0] = int(g[0], 16) ^ int(rc[rounds - 1], 16)

    g[0] = hex(g[0])[2:]

    w4 = [int(w0[num], 16) ^ int(g[num], 16) for num in range(len(g))]
    w5 = [int(w4[num]) ^ int(w1[num], 16) for num in range(len(g))]
    w6 = [int(w5[num]) ^ int(w2[num], 16) for num in range(len(g))]
    w7 = [int(w6[num]) ^ int(w3[num], 16) for num in range(len(g))]

    key_final = []

    for num in range(len(g)):
        w4[num] = hex(w4[num])[2:]
        w5[num] = hex(w5[num])[2:]
        w6[num] = hex(w6[num])[2:]
        w7[num] = hex(w7[num])[2:]

    key_final.append(w4)
    key_final.append(w5)
    key_final.append(w6)
    key_final.append(w7)
    key_final = np.reshape(key_final, (-1, 4))
    key_final = key_final.transpose()

    return key_final


def main():
    base_coeff = [1, 0, 0, 0, 1, 1, 0, 1, 1]  # x**8 + x**4 + x**3 +x+ 1
    base = np.poly1d(base_coeff)

    plaintext = "a5,b6,ca,d4,17,8f,6a,b2,23,3c,d5,61,3a,b4,43,4c"
    key = "23,D5,3A,8F,6A,17,D6,CA,61,B4,4C,B7,3F,B1,1C,D8"
    print("please write the round")
    rounds = int(input())

    state = Text_toArray(plaintext)
    key = Text_toArray(key)
    print("state", state)
    print("-------------")

    key = generateKey(key, rounds)

    state = xor(state, key)
    print("xor", state)
    print("-------------")
    # s box
    statess = []
    for num in state:
        for col in num:
            statess.append(sBox(col))
    state = np.reshape(statess, (-1, 4))
    print("S-box", state)
    print("-------------")

    # shift rows
    i = 0
    for i in range(len(state)):
        d = deque(state[i])
        d.rotate(-i)
        state[i] = list(d)
    print("After Shift", state)
    print("-------------")

    # mix Col
    state = mixcol(state)
    print("The matrix After Mix columns", state)
    print("-------------")


if __name__ == "__main__": main()