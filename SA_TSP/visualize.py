import matplotlib.pyplot as plt
import imageio
import os
if not os.path.isdir(".record"):
    os.mkdir(".record")

f = open("record.txt", "r")
N, T = [int(x) for x in f.readline().split()]
for t in range(T):
    plt.xlabel("x")
    plt.ylabel("y")
    temperture, cost = [float(x) for x in f.readline().split()]
    plt.title(f"T={temperture}, cost={cost}")
    x0, y0 = [int(x) for x in f.readline().split()]
    for i in range(N):
        x1, y1 = [int(x) for x in f.readline().split()]
        plt.plot((x0, x1), (y0, y1), 'ro-')
        x0, y0 = x1, y1
    plt.savefig(f".record/p{t}.jpg");
    plt.cla()

f.close()
frames = []
for t in range(T):
    img = imageio.v2.imread(f'./.record/p{t}.jpg')
    frames.append(img)
imageio.mimsave('./record.gif', frames, fps = 4)
