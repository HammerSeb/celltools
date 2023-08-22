import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt


from linalg.basis import vector, basis, standard_basis
from draw import draw as dr

v1 = [1, 0, 0]
v2 = [0, 1, -1]
v3 = [1, 0, 1]

b = basis(v1,v2,v3)

a = vector([1,1,1], b)
c = vector([0.7,-.873, 1.23], b)


f, ax = dr.make_figure()
dr.draw_point(ax, a, 100, "k")
# dr.draw_point(ax, c, 100, "k")
dr.draw_line(ax, a, c, 3, "grey")

dr.draw_basis(ax, b)


plt.show()