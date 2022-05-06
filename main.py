import numpy
import pyldpc
import bitarray


# Чтение любой строки из файла и перевод в битовую строку
def read_from_file(path):
    with open(path, mode='rb') as file:
        bits = bitarray.bitarray()
        bits.fromfile(file)
        np_bits = numpy.zeros((len(bits))).astype(numpy.int8)
        for i in range(len(bits)):
            np_bits[i] = bits[i]
        return np_bits


# Перевод битовой строки в символы и запись в файл
def write_to_file(path, binary):
    with open(path, mode='wb') as file:
        bits = bitarray.bitarray()
        for i in range(len(binary)):
            bits.append(binary[i])
        bits.tofile(file)


# Добавление нулей до кратности block_size и разделение на блоки длиной block_size
def split(array, block_size):
    if len(array) % block_size != 0:
        array = numpy.append(array, numpy.zeros(block_size - (len(array) % block_size), dtype=numpy.uint8))

    return array.reshape((int(len(array) / block_size), block_size))


class McEliece:
    def __init__(self):
        self.n = 500  # размер выходного блока
        self.k = 0  # размер блока открытого текста
        self.t = 5  # количество ошибок
        self.public_key, self.private_key = self.__keys_generation()

    # Алгоритм декодирования
    # H - parity check matrix (private_key[1][0])
    # c - codeword
    # maxiter - максимальное число итераций для исправления ошибки
    def bit_flipping(self, codeword, maxiter=200):
        import numpy as np
        import warnings

        from random import choice

        n, k = self.private_key[1][0].shape

        Cn = {i: set() for i in range(n)}  # check nodes
        Vn = {i: set() for i in range(k)}  # variable nodes

        for i in range(n):
            for j in range(k):
                if self.private_key[1][0][i][j] == 1:
                    Cn[i].add(j)
                    Vn[j].add(i)

        def fix_errors(c, iteration):
            En = np.zeros(k, dtype=int)  # total error counter
            Sn = np.zeros(n, dtype=int)  # syndrome

            # max depth exceeded -> exit with warning
            if iteration == 0:
                message = "Bit flipping failed with {} iterations.".format(maxiter)
                message += " You may need to increase this number."
                warnings.warn(message)
                return

            # calculate syndrome
            for node, relations in Cn.items():
                Sn[node] = sum(c[list(relations)]) % 2

            # count errors
            for node, check in enumerate(Sn):
                # if check is not satisfied
                if check == 1:
                    for relation in Cn[node]:
                        En[relation] += 1

            # if syndrome is correct -> exit
            if sum(En) == 0:
                return

            # find indexes with max error
            m = 0
            max_indexes = []
            for i, elem in enumerate(En):
                if elem == m and m > 0:
                    max_indexes.append(i)
                elif elem > m:
                    max_indexes = []
                    m = elem
                    max_indexes.append(i)

            # flip one random element of them
            flip = choice(max_indexes)
            c[flip] = 1 if c[flip] == 0 else 0

            # go down recursively
            fix_errors(c, iteration - 1)

        fix_errors(codeword, maxiter)

        return codeword

    # Генерация случайной невырожденной матрицы S(k, k)
    def __s_generation(self):
        while True:
            S = numpy.random.randint(0, 2, (self.k, self.k))

            A = pyldpc.utils.gaussjordan(S)  # преобразование в ступенчатую матрицу
            A = numpy.array(A, dtype=numpy.uint8)

            if (A == numpy.eye(self.k, dtype=numpy.uint8)).all():  # если на главной диагонали 1, то det(S) не ноль
                return S

    # Генерация матриц H и G для LDPC кода, k изменяется
    def __g_generation(self):
        H, G = pyldpc.make_ldpc(self.n, 5, 10, systematic=True)  # G возвращается в транспонированной форме
        self.k = G.shape[1]
        return H, G.T

    # Генерация матрицы перестановок P(n, n)
    def __p_generation(self):
        P = numpy.eye(self.n, dtype=numpy.uint8)  # единичная матрица
        numpy.random.shuffle(P)  # перемешивание с сохранением условия: в каждой строке и столбце один бит 1
        return P

    # Генерация открытого ключа (SGP, t) и закрытого ключа (S, G, P)
    def __keys_generation(self):
        G = self.__g_generation()
        S = self.__s_generation()
        P = self.__p_generation()
        return ((S @ G[1] @ P) % 2, self.t), (S, G, P)  # кортеж (public_key(SGP, t), private_key(S ,G, P))

    # Шифрование одного блока открытым ключом
    def __block_cipher(self, m_block):
        # Случайное число единиц в интервале [0, t) для вектора ошибок
        one_count = numpy.random.randint(0, self.public_key[1])

        # Формирование вектора ошибок z
        z = numpy.ones(one_count, dtype=int)
        z = numpy.append(z, numpy.zeros(self.n - one_count, dtype=numpy.uint8))
        numpy.random.shuffle(z)

        return (m_block @ self.public_key[0] + z) % 2  # m * (S*G*P) + z

    # Дешифрование одного блока закрытым ключом
    def __block_decipher(self, c_block):
        A, inv_p = pyldpc.utils.gaussjordan(self.private_key[2], change=True)  # обратная матрица для P

        c_res = (c_block @ inv_p) % 2  # c_res = c * P^-1
        c_res = numpy.array(c_res).astype(numpy.uint8)

        # Декодирования вектора c_res
        decoded = self.bit_flipping(c_res)
        m_res = decoded[:self.private_key[1][1].shape[0]]

        A, inv_s = pyldpc.utils.gaussjordan(self.private_key[0], change=True)  # обратная матрица для S

        return (m_res @ inv_s) % 2

    def encrypt(self, plain_text):
        encrypted_text = numpy.array([], dtype=numpy.uint8).astype(numpy.uint8)
        blocks_for_cipher = split(plain_text, self.k)  # разбиение на блоки

        for i in range(len(blocks_for_cipher)):  # шифрование блока и конкатенация результата с предыдущей итерацией
            encrypted_text = numpy.append(encrypted_text, self.__block_cipher(blocks_for_cipher[i]))
        return encrypted_text

    def decrypt(self, encrypted_text):
        decrypted_text = numpy.array([]).astype(numpy.uint8)
        blocks_for_decipher = split(encrypted_text, self.n)  # разбиение на блоки

        for i in range(len(blocks_for_decipher)):  # дешифрование блоков, конкатенация результата с предыдущей итерацией
            decrypted_text = numpy.append(decrypted_text, self.__block_decipher(blocks_for_decipher[i]))
        return decrypted_text


def main():
    system = McEliece()
    plain_text = read_from_file("plain_text.txt")
    print("Cipher starting...")
    cipher_text = system.encrypt(plain_text)
    print("Cipher end")
    print("Decipher starting...")
    decrypted_text = system.decrypt(cipher_text)
    write_to_file("decrypted_text.txt", decrypted_text)
    print("Decipher end")


main()
