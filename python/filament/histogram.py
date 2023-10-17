import numpy as np

def histogram(scalar, epsilon):
        scalar.sort()
        index = 0
        N = scalar.size
        x = scalar[index]
        hist = [np.array([x, 0, 0])]
        accummulated = 1./N
        index = 1
        while index < N:
                accummulated += 1./N
                x = scalar[index]
                if (accummulated > epsilon):
                        hist[-1][1] = x
                        hist[-1][2] = accummulated
                        hist.append(np.array([x, 0, 0]))
                        accummulated = 0
                index+=1
        if (accummulated > 0):
                hist[-1][1]=x
                hist[-1][2]=accummulated
        else:
                hist.pop()
        return hist

def histogram_weighted(scalar, weight, epsilon):
        ordering = scalar.argsort()
        scalar = scalar[ordering]
        weight = weight[ordering]
        index = 0
        N = scalar.size
        total_weight = weight.sum()
        weight /= total_weight
        x = scalar[index]
        hist = [np.array([x, 0, 0])]
        accummulated = weight[0]
        index = 1
        while index < N:
                accummulated += weight[index]
                x = scalar[index]
                if (accummulated > epsilon):
                        hist[-1][1] = x
                        hist[-1][2] = accummulated
                        hist.append(np.array([x, 0, 0]))
                        accummulated = 0
                index += 1
        if (accummulated > 0):
                hist[-1][1] = x
                hist[-1][2] = accummulated
        else:
                hist.pop()
        return hist

def histogram2D(field1, field2, nx, ny):
        res = np.zeros((nx * ny, 5))

        xmin = field1.min()
        xmax = field1.max()
        deltax = xmax - xmin
        ymin = field2.min()
        ymax = field2.max()
        deltay = ymax - ymin

        weight = 1 / field1.size
        for i in range(nx):
                for j in range(ny):
                        res[ny * i + j, 0] = xmin + deltax * i / nx
                        res[ny * i + j, 1] = ymin + deltay * j / ny
                        res[ny * i + j, 2] = xmin + deltax * (i + 1) / nx
                        res[ny * i + j, 3] = ymin + deltay * (j + 1) / ny
                        #res[i + nx*j, 4] = ((field1 >= binsX[i]) * (field1 < binsX[i + 1]) * (field2 >= binsY[j]) * (field2 < binsY[j + 1])).sum(axis = None) / N
        for k in range(field1.size):
                i = max(0, min(nx - 1, np.floor(nx * (field1[k] - xmin) / deltax).astype(int)))
                j = max(0, min(ny - 1, np.floor(ny * (field2[k] - ymin) / deltay).astype(int)))
                res[ny * i + j, 4] += weight
        return res

def histogram2D_weighted(field1, field2, weight, nx, ny):
        res = np.zeros((nx * ny, 5))

        xmin = field1.min()
        xmax = field1.max()
        deltax = xmax - xmin
        ymin = field2.min()
        ymax = field2.max()
        deltay = ymax - ymin

        total = weight.sum(axis = None)
        weight /= total
        for i in range(nx):
                for j in range(ny):
                        res[ny * i + j, 0] = xmin + deltax * i / nx
                        res[ny * i + j, 1] = ymin + deltay * j / ny
                        res[ny * i + j, 2] = xmin + deltax * (i + 1) / nx
                        res[ny * i + j, 3] = ymin + deltay * (j + 1) / ny
                        #res[i + nx*j, 4] = weight[(field1 >= binsX[i]) * (field1 < binsX[i + 1]) * (field2 >= binsY[j]) * (field2 < binsY[j + 1])].sum(axis = None)
        #indicesX = np.floor(nx * (field1 - xmin) / deltax)
        #indicesY = np.floor(ny * (field2 - ymin) / deltay)
        for k in range(field1.size):
                i = max(0, min(nx - 1, np.floor(nx * (field1[k] - xmin) / deltax).astype(int)))
                j = max(0, min(ny - 1, np.floor(ny * (field2[k] - ymin) / deltay).astype(int)))
                res[ny * i + j, 4] += weight[k]
        return res