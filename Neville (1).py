def Neville(dx, dy, x):
    length = len(dx)
    t= length * [0]
    for k in range(length):
        for i in range(length - k):
            if k == 0:
                t[i] = dy[i]
            else:
                t[i] = ((x - dx[i + k]) * t[i] + \
                        (dx[i] - x) * t[i + 1]) / \
                       (dx[i] - dx[i + k])
    return t[0]
def main():
    dx = [1, 1.2, 1.3, 1.4]#from example
    dy = [0, 0.112463, 0.167996, 0.222709] #from example
    x = 1.28
    print(Neville(dx, dy, x))
main()