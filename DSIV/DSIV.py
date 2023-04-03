from math import log2, floor
import numpy as np
from scipy.spatial.distance import hamming
from prettytable import PrettyTable
import heapq


class Node:
    def __init__(self, freq, symbol, left=None, right=None):
        self.freq = freq
        self.symbol = symbol
        self.left = left
        self.right = right
        self.huff = ''

    def __lt__(self, nxt):
        return self.freq < nxt.freq


def print_nodes(node, val=''):
    # huffman code for current node
    new_val = val + str(node.huff)

    if node.left:
        print_nodes(node.left, new_val)
    if node.right:
        print_nodes(node.right, new_val)
    if not node.left and not node.right:
        print(f"{node.symbol} -> {new_val}, lc = {node.freq * len(new_val)}\tli(kraft)= 2^{len(new_val)} = {2**-len(new_val)}")


def calculate_entropy(output_enable, xlist=[]):
    x_list = []
    n_t = 0
    if len(xlist) > 0:
        x_list = xlist
    else:
        t_i_list = list(map(float, input("Enter τi (Symbol duration)\t").split()))
        x_list = list(map(float, input("Enter space separated Bit Source:\t").split()))
        n_t = len(t_i_list)

    n_x = len(x_list)
    h_x = 0
    t = 0
    if n_x <= 0:
        return
    for i in range(n_x):
        if x_list[i] != 0:
            h_x += x_list[i] * log2(x_list[i])
        if n_t > 0:
            t += t_i_list[i] * x_list[i]
    else:
        h_x *= -1
    h_0 = log2(n_x)
    R = h_0 - h_x
    r = R / h_0

    if output_enable:
        print(f"Entropie der Quelle H(X) = -Σ_(i=1)^n P(xi) * ld(P(x_i))bit = {round(h_x, 3)}bit\n"
              f"Entscheidungsgehalt H_0 = ld(n) =  {round(log2(n_x),3)} bit\n"
              f"Redundanz R = H_0 - H(X) =  {round(R,3)} bit\n"
              f"relative Redundanz r = R/H_0 =  {round(r, 3)} = {round(r*100, 3)}%")
        print("Informationsgehalt der Zeichen(-ld(x_i)):")
        for i, item in enumerate([-log2(i) for i in x_list]):
            print(f'l(x{i+1}) = {round(item,3)} bit\t')


        if n_t > 0:
            print(f"Mittlere Symboldauer τ =Σ_(i=1)^n τi * P(xi)= {t} s\n"
                  f"Informationsfluss ΦI =H(X)/t =  {h_x/t} bit/s")
    return h_x


def calculate_entropy_dual(n=0, p_x_y=0) -> None:
    H_X_com_Y = 0.0
    H_X_sem_Y = 0.0
    p_x = []
    p_y = []

    if n == 0 and p_x_y == 0:
        n = int(input("Enter n\t"))
        p_x_y = []
        for i in range(n):
            p_x_y.append(list(map(float, input(f"Enter space separated Bit Source for x{i+1} to y1 to y{n} ").split())))

    for i in range(n):
        p_x.append(round(sum(p_x_y[i]), 4))
    for i in range(n):
        tmp = 0.0
        for j in range(n):
            tmp += p_x_y[j][i]
        p_y.append(round(tmp, 4))

    H_X = calculate_entropy(False, p_x)
    H_Y = calculate_entropy(False, p_y)

    for i in range(n):
        for j in range(n):
            if p_x_y[i][j] != 0:
                H_X_com_Y += p_x_y[i][j] * log2(p_x_y[i][j])
                H_X_sem_Y += p_x_y[i][j] * log2(p_x_y[i][j] / (p_x[i] * p_y[j]))
    H_X_com_Y *= -1
    H_X_vert_Y = H_X_com_Y - H_Y
    H_Y_vert_X = H_X_com_Y - H_X

    table = PrettyTable(["P(Xi,Yj)", *['y' + str(i+1) for i in range(n)], "P(Xi)"])
    table.align = "l"

    for i in range(n):
        table.add_row([f"x{i+1}", *p_x_y[i], p_x[i]])
    table.add_row(["P(yj)", *p_y, "1"])
    print(table)
    print(f'Entropie H(X) = -Σ_(i=1)^m P(xi) * ld(P(x_i))bit= {round(H_X,4)} bit\n'
          f'Entropie H(Y) = -Σ_(i=1)^n P(yi) * ld(P(y_i))bit =  {round(H_Y,4)} bit\n'
          f'Verbundentropie H(X,Y) =-Σ_(i=1)^m Σ_(j=1)^n P(xi,yj) * ld(P(xi,yj)) bit= {round(H_X_com_Y,4)} bit\n'
          f'Äquivokation H(X|Y) = -Σ_(i=1)^m Σ_(j=1)^n P(xi,yj) * ld(P(xi|yj)) bit = {round(H_X_vert_Y,4)} bit\n'
          f'Irrelevanz H(Y|X) = -Σ_(i=1)^m Σ_(j=1)^n P(xi,yj) * ld(P(yj|xi)) bit =  {round(H_Y_vert_X,4)} bit\n'
          f'Transinformation (BSC: Kanalkapazität) H(X;Y) = Σ_(i=1)^m Σ_(j=1)^n P(xi,yj) * ld(P(xi,yj)/(P(xi)* P(yj))) bit = {round(H_X_sem_Y, 4)} bit\n')
    return


def bit_symmetric_channel() -> None:
    p = float(input("Enter error propability:\t"))
    q = float(input("Enter Symbol probability:\t"))

    print(f'\nx1y1 = (1 - p) * q = (1 - {p}) * {q} = {(1-p)*q}')
    print(f'x1y2 = p * q = {p} * {q} = {p*q}')
    print(f'x2y1 = p * ( 1 - q) = {p} * (1 - {q}) = {p*(1-q)}')
    print(f'x2y2 = (1-p)*(1-q) = (1 - {p}) * (1 - {q}) = {(1-p)*(1-q)}\n')

    calculate_entropy_dual(2, [[(1-p)*q, p*q], [p*(1-q), (1-p)*(1-q)]])
    return


def binary_markoff_source() -> None:
    left_iterator = [1, 2, 1, 2]
    right_iterator = [1, 1, 2, 2]

    i = 0
    while i < 4:
        try:
            temp = float(input(f"Enter P(x{left_iterator[i]}|x{right_iterator[i]}) or enter to skip variable:\t"))

        except ValueError:
            print("skipped")
            i += 1
        else:
            match i:
                case 0:
                    px1_x1 = temp
                    i += 2
                case 1:
                    px1_x1 = 1 - temp
                    i += 1
                case 2:
                    px1_x2 = temp
                    i += 2
                case 3:
                    px1_x2 = 1-temp
                    i += 1


    p = np.array([[1-px1_x1, -px1_x2],[1,1]])
    k = np.array([[0], [1]])
    px = np.linalg.solve(p, k)
    hx_x1 = round(-px1_x1 * log2(px1_x1) - (1 - px1_x1) * log2(1-px1_x1),4)
    hx_x2 = round(-px1_x2 * log2(px1_x2) - (1 - px1_x2) * log2(1 - px1_x2),4)
    hx_inf = round(float(px[0] * hx_x1 + px[1] * hx_x2),4)
    print(f'P(x1) = {round(float(px[0]), 4)}\nP(x2) = {round(float(px[1]), 4)}')
    print(f'H(X|x_1) = -P(x1|x1) * ld(P(x1|x1) - (1 - P(x1|x1)) * ld(1-P(x1|x1) = {hx_x1}')
    print(f'H(X|x_2)= -P(x1_x2) * ld(P(x1|x2)) - (1 - P(x1|x2)) * ld(1 - P(x1|x2) = {hx_x2}')
    print(f'H∞(X) = P(x1) * H(X|x1) + P(x2) * H(X|x2) = {hx_inf}')
    print(f'H0 = ld(n) = 1')
    print(f'Die Quelle hat {f"ein Gedächtnis, weil H∞(X) < H0" if hx_inf < 1 else "kein Gedächtnis, weil H∞(X) >= H0" }')
    return


def huffman_code() -> None:
    freq = list(map(float, input("Enter space separated Bit Source from x_1 to x_n:\t").split()))
    chars = [f'x{i + 1}' for i in range(len(freq))]
    nodes = []
    for x in range(len(chars)):
        heapq.heappush(nodes, Node(freq[x], chars[x]))

    while len(nodes) > 1:
        left = heapq.heappop(nodes)
        right = heapq.heappop(nodes)

        left.huff = 0
        right.huff = 1

        newNode = Node(left.freq + right.freq, left.symbol + right.symbol, left, right)

        heapq.heappush(nodes, newNode)

    print_nodes(nodes[0])


def codegen_matrix() -> None:
    # gen_matrix = [[1, 0, 0, 0, 1, 0, 1], [0, 1, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1, 0], [0, 0, 0, 1, 0, 1, 1]]
    # info = [0, 1, 0, 1]

    gen_matrix = []
    code = []
    n = int(input("Enter height of generator matrix:\t"))
    for i in range(n):
        gen_matrix.append(list(map(int, input(f"Enter space separated Bits for generator matrix line {i} ").split())))
    print(f'Generator Matrix:\n{gen_matrix}')
    while True:
        info = list(map(int, [*input("Enter the message, press enter to return\t")]))
        if len(info) != len(gen_matrix):
            print("Dimensional error")
            pass
        if len(info) == 0:
            return
        for i in range(len(gen_matrix[0])):
            temp = 0
            for j in range(len(info)):
                temp += (gen_matrix[j][i] == 1 and info[j] == 1)
            code.append(temp % 2)
        print(code)
        code = []


def crc() -> None:
    n_code = int(input("Enter code length:\t"))
    n_info = int(input("Enter information length:\t"))
    gen_poly = list(map(int, [*input("Enter generator polynom \t")]))
    info_vector = []
    unsystematic_code_vector = []
    systematic_code_vector = [[0 for j in range(n_code)]for i in range(2**n_info)]
    unsystematic_generator_matrix = []
    systematic_generator_matrix = []
    systematic_table = PrettyTable(["Nr.", "Informationsvektor x", "Codewort c"])
    unsystematic_table = PrettyTable(["Nr.", "Informationsvektor x", "Codewort c"])

    hamming_table = PrettyTable([' ', *[i for i in range(n_code)]])
    unsys_gen_poly = list(gen_poly)
    while len(unsys_gen_poly) < n_code:
        unsys_gen_poly.append(0)
    unsys_gen_poly = np.array(unsys_gen_poly)

    for i in range(n_info):
        unsystematic_generator_matrix.append([unsys_gen_poly])
        unsys_gen_poly = np.roll(unsys_gen_poly, 1)

    for i in range(2**n_info):
        info_vector.append([int(x) for x in list('{0:0b}'.format(i))])

        unsystematic_code_vector.append((np.polymul(info_vector[i], gen_poly) % 2).tolist())  # Polynom multiplication in respect to GF(2)
        while len(info_vector[i]) < n_info:  # fill with zeroes
            info_vector[i].insert(0, 0)
        while len(unsystematic_code_vector[i]) < n_code:  # fill with zeroes
            unsystematic_code_vector[i].insert(0, 0)
        unsystematic_table.add_row([i, str(info_vector[i]), str(unsystematic_code_vector[i])])
        a = int("".join(str(x) for x in systematic_code_vector[i][n_info:]), 2)
        systematic_code_vector[int("".join(str(x) for x in unsystematic_code_vector[i][:n_info]), 2)] = unsystematic_code_vector[i]

    for i in range(2**n_info):
        systematic_table.add_row([i, str(info_vector[i]), str(systematic_code_vector[i])])

    print("\nNicht-systematische Generatormatrix:")
    for row in unsystematic_generator_matrix:
        print(row)

    print("\nNicht-systematische Codetabelle:")
    print(unsystematic_table)

    print("\nSystematische Generatormatrix: Unsystemattische Generator Matrix mit Gauss in die Form Einheitsmatrix | P bringen")



    print("\nSystematische Codetabelle:")
    print(systematic_table)

    calculate_hamming_distance(unsystematic_code_vector)
    return


def calculate_hamming_distance(code_vector=[]) -> None:
    if len(code_vector) == 0:
        n_code = int(input("Enter code length:\t"))
        for i in range(n_code):
            code_vector.append(list(map(int, input(f"Enter space separated Code for line{i}:\t").split())))
    else:
        n_code = len(code_vector[0])
    hamming_table = PrettyTable([' ', *[f'C{i}' for i in range(n_code)]])
    hamming_distance = 10000
    hamming_distance_table = [[0 for j in range(n_code)] for i in range(n_code)]
    for i in range(n_code):
        for j in range(n_code):
            temp = hamming(code_vector[i], code_vector[j]) * n_code
            hamming_distance_table[i][j] = temp
            if hamming_distance > temp > 0:
                hamming_distance = temp

    for i, row in enumerate(hamming_distance_table):
        row.insert(0, f'C{i}')
        hamming_table.add_row([*row])
    print(f'Hamming Table:\n{hamming_table}')
    print(f'\nHamming-Distance d_min= {hamming_distance}')
    print(f't_erk = ⌊d_min-1⌋ = {floor(hamming_distance-1)}')
    print(f't_korr = ⌊(d_min-1)/2⌋ = {floor((hamming_distance-1)/2)}')


def show_menu() -> None:
    print("\n"
          "0 Entropy\n"
          "1 Dual Entropy\n"
          "2 Bit Symmetric Channel BSC\n"
          "3 Binary Markoff Source\n"
          "4 Huffman Code\n"
          "5 Code Generation with Code Matrix\n"
          "6 CRC Calculation\n"
          "7 Hamming Distance Calculation\n"

          "Press 'c' to close")


def main() -> None:
    choice = 'u'
    while choice != 'c':
        show_menu()
        choice = input()

        match choice:
            case '0':
                calculate_entropy(True)
            case '1':
                calculate_entropy_dual()
            case '2':
                bit_symmetric_channel()
            case '3':
                binary_markoff_source()
            case '4':
                huffman_code()
            case '5':
                codegen_matrix()
            case '6':
                crc()
            case '7':
                calculate_hamming_distance()
            case 'c':
                pass
            case _:
                print("Input Error")


if __name__ == "__main__":
    main()

