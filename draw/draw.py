import numpy as np
from matplotlib import pyplot as plt

from . import vector, basis, standard_basis
from . import Arrow3D

def make_figure(axis="on", xlim=(-2,2), ylim=(-2,2), zlim=(-2,2)):
    """
    create 3D plot instance
    :return: figure, axes
    """
    f = plt.figure()
    ax = f.add_subplot(projection="3d")
    if axis=="off":
        ax.set_axis_off()
    elif axis=="on":
        pass
    else:
        raise ValueError("axis must be \"on\" or \"off\"")
        return
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    f.tight_layout()
    return f, ax

def draw_point(ax, coord, s, c):
    """
    draws point at position given by coords
    :param ax: axes instance
    :param coord: vector instance or list/tuple size 3
    :param s: size for axes.scatter
    :param c: color for axes.scatter
    """
    if isinstance(coord, vector):
        ax.scatter(*coord.global_coord, color=c, s=s)
    else:
        ax.scatter(coord[0], coord[1], coord[2], color=c, s=s)

def draw_line(ax, coord1, coord2, lw, c):
    """
    draws line between positions given by coords1 and coords2
    :param ax: axes instance
    :param coord1: vector instance or list/tuple size 3
    :param coord2: vector instance or list/tuple size 3
    :param lw:  line width for axes.plot
    :param c: color for axes.scatter
    """
    if isinstance(coord1, vector) and isinstance(coord2, vector):
        x = [coord1.global_coord[0], coord2.global_coord[0]]
        y = [coord1.global_coord[1], coord2.global_coord[1]]
        z = [coord1.global_coord[2], coord2.global_coord[2]]
    elif type(coord1) != type(coord2):
        raise ValueError("coordinates must be of same type")
    else:
        x = [coord1[0], coord2[0]]
        y = [coord1[1], coord2[1]]
        z = [coord1[2], coord2[2]]
    ax.plot(x, y, z, c=c, lw=lw)

def draw_vector(ax, vec, offset=None, c="r", lw=3):
    """
    draw vector arrow in direction of vec. Offset specifies origin
    :param ax: axes instance
    :param vec: vector instance
    :param offset: vector instance
    """
    if isinstance(offset, np.ndarray):
        if isinstance(vec,vector):
            ax.arrow3d(offset[0], offset[1], offset[2], vec.global_coord[0], vec.global_coord[1], vec.global_coord[2],
                       mutation_scale=20, arrowstyle="-|>", color=c, linewidth=lw)
        else:
            raise ValueError("input needs to be instance of vector class")
    else:
        if isinstance(vec,vector) and isinstance(offset, vector):
            ax.arrow3d(offset.global_coord[0], offset.global_coord[1], offset.global_coord[2],
                       vec.global_coord[0], vec.global_coord[1], vec.global_coord[2],
                       mutation_scale=20, arrowstyle="-|>", color=c, linewidth=lw)
        else:
            raise ValueError("input needs to be instance of vector class")

def draw_basis(ax, basis, c="grey", lw=3):
    draw_vector(ax,vector(basis[0]), offset=basis.offset, c=c, lw=lw)
    draw_vector(ax,vector(basis[1]), offset=basis.offset, c=c, lw=lw)
    draw_vector(ax,vector(basis[2]), offset=basis.offset, c=c, lw=lw)


def draw_frame(ax, basis, c="black", lw=1.5):
    _origin = vector(basis.offset)
    _a, _b, _c = vector(basis.basis[0]), vector(basis.basis[1]), vector(basis.basis[2])
    draw_line(ax, _origin, _a+_origin, c=c, lw=lw)
    draw_line(ax, _origin, _b+_origin, c=c, lw=lw)
    draw_line(ax, _a+_origin, _a+_b+_origin, c=c, lw=lw)
    draw_line(ax, _b+_origin, _a+_b+_origin, c=c, lw=lw)
    draw_line(ax, _origin, _c+_origin, c=c, lw=lw)
    draw_line(ax, _a+_origin, _a+_c+_origin, c=c, lw=lw)
    draw_line(ax, _b+_origin, _b+_c+_origin, c=c, lw=lw)
    draw_line(ax, _a+_b+_origin, _a+_b+_c+_origin, c=c, lw=lw)
    draw_line(ax, _c+_origin, _a+_c+_origin, c=c, lw=lw)
    draw_line(ax, _c+_a+_origin, _c+_a+_b+_origin, c=c, lw=lw)
    draw_line(ax, _c+_b+_origin, _c+_b+_a+_origin, c=c, lw=lw)
    draw_line(ax, _c+_origin, _b+_c+_origin, c=c, lw=lw)