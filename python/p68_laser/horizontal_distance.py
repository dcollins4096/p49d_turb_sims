from dtools.starter1 import *

def pairwise(iterable):
    from itertools import tee

    a, b = tee(iterable)
    _ = next(b, None)
    yield from zip(a, b)


def ho(a,b):
    intersections = []

    for x1, y1 in enumerate(a):
        x2 = len(b)
        y2 = y1
        for x3, (y3, y4) in enumerate(pairwise(b)):
            x4 = x3 + 1

            try:
                t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
                u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
            except ZeroDivisionError:
                continue

            if 0 <= t <= 1.0 and 0 <= u <= 1.0:
                px, py = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
                intersections.append((x1, px, py))
                break
    return intersections

def try2(a,b):
    intersections = ho(a,b)
    plt.clf()
    plt.plot(range(len(a)), a, color="blue")
    plt.plot(range(len(b)), b, color="orange")
    for i, x, y in intersections:
        xs = [i, x]
        ys = [a[i], y]
        plt.plot(xs, ys, "r--")
        #plt.plot(x, y, "r+")
    plt.savefig('plots_to_sort/hor_test_2')
    return intersections
def try1():
    a = [4, 1, 2, 7, 8, 8, 6, 11, 7, 10, 11, 15, 14, 14, 13, 17, 17, 21, 22, 20]
    b = [3, 0, 1, 6, 3, 6, 9, 11, 8, 8, 11, 15, 14, 15, 17, 14, 18, 17, 18, 20]
    intersections = ho(a,b)
    plt.clf()
    plt.plot(range(len(a)), a, color="blue")
    plt.plot(range(len(b)), b, color="orange")
    for i, x, y in intersections:
        xs = [i, x]
        ys = [a[i], y]
        plt.plot(xs, ys, "r--")
        #plt.plot(x, y, "r+")
    plt.savefig('plots_to_sort/hor_test')
    return intersections
