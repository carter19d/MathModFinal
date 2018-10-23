import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

im = Image.open("toporecolor.png")
darkness = list(im.getdata(0))
darknessA = 1-np.array(darkness, int).reshape(im.size[::-1])
print(im.size)

print(darknessA)
#plt.matshow(darknessA)

density = np.zeros(im.size[::-1])

size = 50

for x in range(size, im.size[0] - size):
    for y in range(size, im.size[1] - size):
        density[y, x] = sum(sum(darknessA[y-size:y+size+1, x-size:x+size+1]))

plt.matshow(density)
plt.show()
