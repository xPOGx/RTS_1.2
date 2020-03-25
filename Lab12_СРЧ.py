import math
import time
import matplotlib.pyplot as plt
import random
import numpy as np
#Dzhevaha
n_harm = 10; w_freq = 1500; N = 256
#Gorobez
#n_harm = 10; w_freq = 2700; N = 256

x_t = list(range(0, N))


def generate_rand_signal(n):
    xt = list(range(0, n))
    start = time.time()
    time.sleep(1)
    harmonics = []
    for i in range(n_harm):
        A = random.randint(1, 100)
        fi = random.randint(1, 100)
        harmonics.append(list(map(
                        lambda x: A*math.sin(w_freq/(i+1)*xt.index(x) + fi),
            xt)))
    smth = list(zip(*harmonics))
    signal = []
    for el in smth:
        signal.append(sum(el))
    end = time.time() - start - 1
    return signal, harmonics, end


def get_Mx(rand_signal: list):
    start = time.time()
    time.sleep(1)
    Mx = sum(rand_signal)/N
    end = time.time() - start - 1
    return Mx, end


def get_Dx(rand_signal: list, mx: float):
    start = time.time()
    time.sleep(1)
    Dx = 0
    for i in range(N):
        Dx += (rand_signal[i] - mx)**2
    Dx = Dx/(N-1)
    end = time.time() - start - 1
    return Dx, end


def get_Rxx_func(rand_signal: list, mx: float, dx: float):
    start = time.time()
    time.sleep(1)
    func = []
    for tau in range((N - 1) // 2):
        cov_xx = 0
        for t in range((N-1)//2):
            cov_xx += (rand_signal[t] - mx)*(rand_signal[t+tau] - mx)
        cov_xx /= (N - 1)
        Rxx = cov_xx / (math.sqrt(dx) ** 2)
        func.append(Rxx)
    end = time.time() - start - 1
    return func, end


def get_Rxy_func(rand_signal1: list, mx1: float, dx1: float,
            rand_signal2: list, mx2: float, dx2: float):
    start = time.time()
    time.sleep(1)
    func = []
    for tau in range((N - 1) // 2):
        cov_xy = 0
        for t in range((N-1)//2):
            cov_xy += (rand_signal1[t] - mx1) * (rand_signal2[t + tau] - mx2)
        cov_xy /= (N - 1)
        Rxy = cov_xy /(math.sqrt(dx1)*math.sqrt(dx2))
        func.append(Rxy)
    end = time.time() - start - 1
    return func, end


def labs1_2():
    # Calculating
    (res1, harmonics1, t1) = generate_rand_signal(N)
    (res2, harmonics2, t2) = generate_rand_signal(N)
    (res3, harmonics3, t3) = generate_rand_signal(N)

    (Mx1, Mx1_time) = get_Mx(res1)
    (Mx2, Mx2_time) = get_Mx(res2)
    (Mx3, Mx3_time) = get_Mx(res3)
    (Dx1, Dx1_time) = get_Dx(res1, Mx1)
    (Dx2, Dx2_time) = get_Dx(res2, Mx2)
    (Dx3, Dx3_time) = get_Dx(res3, Mx3)

    Rxy = get_Rxy_func(res1, Mx1, Dx1, res2, Mx2, Dx2)[0]
    Rxx = get_Rxx_func(res3, Mx3, Dx3)[0]

    # Plotting
    # Plot no.1: Harmonics no.1 & random signal no.1 -> lab2_harm1.png
    plt.figure(1)
    plt.subplot(311)
    plt.title("Перша гармоніка")
    for i in range(len(harmonics1)):
        plt.plot(x_t, harmonics1[i])

    plt.subplot(313)
    plt.title("Перший рандомний сигнал")
    plt.plot(x_t, res1)
    plt.savefig("./lab12_Агарм.png")
    plt.close(1)

    # Plot no.2: Harmonics no.2 & random signal no.2 -> lab2_harm2.png
    plt.figure(2)
    plt.subplot(311)
    plt.title("Друга гармоніка")
    for i in range(len(harmonics2)):
        plt.plot(x_t, harmonics2[i])

    plt.subplot(313)
    plt.title("Другий рандомний сигнал")
    plt.plot(x_t, res2)
    plt.savefig("./lab12_Бгарм.png")
    plt.close(2)

    # Plot no.3: Random signal no.1 & random signal no.2 -> lab2_compare.png
    plt.figure(3)
    plt.title("Порівняння двох рандомних сигналів")
    plt.subplot(311)
    plt.title("Перший сигнал")
    plt.plot(x_t, res1)
    plt.subplot(312)
    plt.title("Другий сигнал")
    plt.plot(x_t, res2)
    plt.subplot(313)
    plt.title("Кореляція зі Т")
    plt.plot(list(range(len(Rxy))), Rxy)
    plt.savefig("./lab12_порів.png")
    plt.close(3)

    # Plot no.4: Random signal no.1 & random signal no.2 -> lab2_compare.png
    plt.figure(4)
    plt.title("Автокореляція")
    plt.subplot(311)
    plt.title("Третій сигнал")
    plt.plot(x_t, res1)
    plt.subplot(313)
    plt.title("Автокореляція зі Т")
    plt.plot(list(range(len(Rxx))), Rxx)
    plt.savefig("./lab12_автокор.png")
    plt.close(4)

    # Text to write into a file result.txt
    text = [
        "# Перший сигнал ({0:.5} s) Параметри:\n\tMx1 = {0:<10.4}\t Tm1: {1:.5} сек.\n"
        "\tDx1 = {2:<10.6}\t Td1: {3:.5} сек.\n".format(Mx1, Mx1_time,
                                                           Dx1, Dx1_time, t1),
        "# Другий сигнал ({0:.5} s) Параметри:\n\tMx2 = {0:<10.4}\t Tm2: {1:.5} сек.\n"
        "\tDx2 = {2:<10.6}\t Td2: {3:.5} сек.\n".format(Mx2, Mx2_time,
                                                           Dx2, Dx2_time, t2),
        "# Третій сигнал ({0:.5} s) Параметри:\n\tMx3 = {0:<10.4}\t Tm3: {1:.5} сек.\n"
        "\tDx3 = {2:<10.6}\t Td3: {3:.5} сек.\n".format(Mx3, Mx3_time,
                                                           Dx3, Dx3_time, t3)
    ]

    # Writing into a file
    f = open("./lab12_висновки.txt", "w")
    for i in text:
        f.write(i)
    f.close()


if __name__ == '__main__':
    labs1_2()
