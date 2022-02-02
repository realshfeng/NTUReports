import matplotlib.pyplot as plt


def derivative_e(concentration_e, concentration_s, concentration_es):
    return 150 * concentration_es + 600 * concentration_es - 100 * concentration_e * concentration_s


def derivative_s(concentration_e, concentration_s, concentration_es):
    return 600 * concentration_es - 100 * concentration_e * concentration_s


def derivative_es(concentration_e, concentration_s, concentration_es):
    return 100 * concentration_e * concentration_s - 150 * concentration_es - 600 * concentration_es


def derivative_p(concentration_es):
    return 150 * concentration_es


def runge_kutta(start, end, step, concentration_e, concentration_s, concentration_es, concentration_p):
    n = int((end - start) / step)
    time = start

    time_point = []
    e_list = []
    s_list = []
    es_list = []
    p_list = []
    for i in range(1, n + 1):
        # slope for all four substrates
        k1_e = derivative_e(concentration_e, concentration_s, concentration_es)
        k1_s = derivative_s(concentration_e, concentration_s, concentration_es)
        k1_es = derivative_es(concentration_e, concentration_s, concentration_es)
        k1_p = derivative_p(concentration_es)

        k2_e = derivative_e(concentration_e + step * k1_e / 2, concentration_s + step * k1_s / 2,
                            concentration_es + step * k1_es / 2)
        k2_s = derivative_s(concentration_e + step * k1_e / 2, concentration_s + step * k1_s / 2,
                            concentration_es + step * k1_es / 2)
        k2_es = derivative_es(concentration_e + step * k1_e / 2, concentration_s + step * k1_s / 2,
                              concentration_es + step * k1_es / 2)
        k2_p = derivative_p(concentration_es + step * k1_es / 2)

        k3_e = derivative_e(concentration_e + step * k2_e / 2, concentration_s + step * k2_s / 2,
                            concentration_es + step * k2_es / 2)
        k3_s = derivative_s(concentration_e + step * k2_e / 2, concentration_s + step * k2_s / 2,
                            concentration_es + step * k2_es / 2)
        k3_es = derivative_es(concentration_e + step * k2_e / 2, concentration_s + step * k2_s / 2,
                              concentration_es + step * k2_es / 2)
        k3_p = derivative_p(concentration_es + step * k2_es / 2)

        k4_e = derivative_e(concentration_e + step * k3_e, concentration_s + step * k3_s,
                            concentration_es + step * k3_es)
        k4_s = derivative_s(concentration_e + step * k3_e, concentration_s + step * k3_s,
                            concentration_es + step * k3_es)
        k4_es = derivative_es(concentration_e + step * k3_e, concentration_s + step * k3_s,
                              concentration_es + step * k3_es)
        k4_p = derivative_p(concentration_es + step * k3_es)

        concentration_e = concentration_e + (step / 6.0) * (k1_e + 2 * k2_e + 2 * k3_e + k4_e)
        concentration_s = concentration_s + (step / 6.0) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s)
        concentration_es = concentration_es + (step / 6.0) * (k1_es + 2 * k2_es + 2 * k3_es + k4_es)
        concentration_p = concentration_p + (step / 6.0) * (k1_p + 2 * k2_p + 2 * k3_p + k4_p)
        print((step / 6.0) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s))
        if (step / 6.0) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s) >= -0.01:
            print(i)
            print("concentration_e = {}, concentration_s = {}, concentration_es = {}, concentration_p = {}".
                  format(concentration_e, concentration_s, concentration_es, concentration_p))
            print("concentration_p = {}, reaction velocity={}".
                  format(concentration_p, derivative_p(concentration_es)))
            return
        time = time + step
        p_list.append(derivative_p(concentration_es))
        time_point.append(time)
        plt.plot(time_point, p_list)
        plt.xlabel("time/min")
        plt.ylabel("velocity/ug/min")
        plt.show()


start = 0.0
end = 10000.0
step = 0.001

initial_e = 1.0
initial_s = 10.0
initial_es = 0.0
initial_p = 0.0

runge_kutta(start, end, step, initial_e, initial_s, initial_es, initial_p)
