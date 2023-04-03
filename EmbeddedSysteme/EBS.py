from math import log2, sqrt, ceil, floor
from prettytable import PrettyTable


def manchester_code() -> None:
    message = list(map(int, [*input("Nachricht eingeben:\t")]))
    manchester_message = []
    for element in message:
        manchester_message.append(False ^ element)
        manchester_message.append(True ^ element)
    print(f'Manchester Code: {manchester_message}')


def hamming_code_generator() -> None:
    message = list(map(int, [*input("Nachricht eingeben:\t")]))
    hamming_message = list(message)
    n = len(message)
    k = 0
    while 2 ** k < k + n + 1:
        k += 1
    N = k + n
    table = PrettyTable([i + 1 for i in range(N)])
    for i in range(N):
        if i == 0 or ceil(log2(1 + i)) == floor(log2(1 + i)):
            hamming_message.insert(i, f'p{i + 1}')
    table.add_row(hamming_message)
    print(table)

    print("Nummern der Bits mit 1")
    first_match = True
    for i, item in enumerate(hamming_message):
        if item == 1:
            if first_match:
                first_match = False
                parity = i + 1
            else:
                parity = parity ^ i + 1
            print(f'1 bei Stelle {i + 1:02} --> {format(i + 1, "08b")}')
    print(f'Pariätsbits: {format(parity, "08b")}')
    parity_bits = [item for item in bin(parity)[::-1]]
    j = 0
    for i in range(N):
        if isinstance(hamming_message[i], str):
            hamming_message[i] = parity_bits[j]
            j += 1
    table.add_row(hamming_message)

    print(table)

def hamming_code_check() -> None:
    encoded_message = list(map(int, [*input("Nachricht eingeben:\t")]))
    N = len(encoded_message)
    table = PrettyTable([i + 1 for i in range(N)])
    table.add_row(encoded_message)
    parity_bits = []
    for i in range(ceil(sqrt(N))):
        parity_bits.append(encoded_message[(2 ** i) - 1])
    parity_bits.reverse()
    for i in range(N):
        if i == 0 or ceil(log2(1 + i)) == floor(log2(1 + i)):
            encoded_message[i] = f'p{i + 1}'
    table.add_row(encoded_message)

    first_match = True
    for i, item in enumerate(encoded_message):
        if item == 1:
            if first_match:
                first_match = False
                message_parity = i + 1
            else:
                message_parity = message_parity ^ i + 1
            print(f'1 bei Stelle {i + 1:02} --> {format(i + 1, "08b")}')
            message_parity_list = [item for item in bin(message_parity)[::-1]]
            del message_parity_list[-2:]
            message_parity_list = list(map(int, message_parity_list))
    message_parity_list.reverse()
    print(table)
    print("Parität der übertragenen Nachricht:\t", *parity_bits)
    print("Berechnete Parität der Nachricht:\t", *message_parity_list)

    error = int("".join(str(x) for x in message_parity_list), 2) ^ int("".join(str(x) for x in parity_bits), 2)
    if error == 0:
        print("Nachricht korrekt")
    else:
        print(f'Fehler in Bit {error}')


def crc_validation() -> None:
    gen_poly = list(map(int, [*input("Enter generator polynom:\t")]))
    message = list(map(int, [*input("Nachricht eingeben:\t")]))
    crc_code = [(np.polymul(message, gen_poly) % 2).tolist()]
    print(f"CRC code: {crc_code}")


def four_b_five_b(mode: bool) -> None:
    dict4b5b= {"0000": "11110", "0100": "01010", "1000":"10010", "1100":"11010", "0001":"01001", "0101":"01011",
               "1001":"10011", "1101":"11011", "0010":"10100", "0110":"01110", "1010":"10110", "1110":"11100",
               "0011":"10101", "0111":"01111", "1011":"10111", "1111":"11101"}

    print(f'{"Encode 4B/5B" if mode else "Decode 5B/4B"}')
    message = list(map(str, input(f"Enter Message in {4 if mode else 5} Bit Packets, each space separated)\t").split()))
    output = []
    for item in message:
        if mode:
            output.append(f'{dict4b5b[item]}, ')
        else:
            output.append(list(dict4b5b.keys())[list(dict4b5b.values()).index(item)])
    print(output)

def show_menu() -> None:
    print("\n"
          "0 Manchester Code\n"
          "1 Hamming-Code Generator\n"
          "2 Hamming-Code Check\n"
          "3 CRC Validation\n"
          "4 4B/5B Encoder\n"
          "5 5B/4B Decoder\n"

          "Press 'c' to close")


def main() -> None:
    choice = 'u'
    while choice != 'c':
        show_menu()
        choice = input()

        match choice:
            case '0':
                manchester_code()
            case '1':
                hamming_code_generator()
            case '2':
                hamming_code_check()
            case '3':
                crc_validation()

            case '4':
                four_b_five_b(True)
            case '5':
                four_b_five_b(False)

            case 'c':
                pass
            case _:
                print("Input Error")


if __name__ == "__main__":
    main()
