import math
import numpy as np
from pyx import *

cellsize = 5
origin = [0,0]
ray_enter = 0.5   # x-dimension
normal = [0.7, 0.3]  # both must be positive for this illustration
unit.set(xscale=0.35*cellsize)

# Normalize ray direction vector
normal = np.array(normal)
normal /= math.sqrt((normal**2).sum())

#------------------------------------------------------------------------
# Figure with true covering factor
#------------------------------------------------------------------------

def figure1():
    c = canvas.canvas()

    # Ray properties
    ray_start = np.array([origin[0] + ray_enter*cellsize, origin[1]])
    dx = origin[1] + cellsize - ray_start[0]
    dy = dx * normal[1] / normal[0]
    ray_end = np.array( [origin[0] + cellsize, origin[1]+dy] )
    ray_center = 0.5*(ray_start + ray_end)

    exterior = 2.0
    ray_head = ray_start - exterior*normal
    ray_tail = ray_end + exterior*normal

    perp = np.array([-normal[1], normal[0]])
    beam_width = 0.4 * cellsize

    # Column from r to r+dr
    ray_area = box.polygon([(ray_start[0] - beam_width*perp[0],
                             ray_start[1] - beam_width*perp[1]),
                            (ray_start[0] + beam_width*perp[0],
                             ray_start[1] + beam_width*perp[1]),
                            (ray_end[0] + beam_width*perp[0],
                             ray_end[1] + beam_width*perp[1]),
                            (ray_end[0] - beam_width*perp[0],
                             ray_end[1] - beam_width*perp[1])])
    covering = box.polygon([(ray_start[0], ray_start[1]),
                            (origin[0]+cellsize, origin[1]),
                            (ray_end[0], ray_end[1]),
                            (ray_end[0] + beam_width*perp[0],
                             ray_end[1] + beam_width*perp[1]),
                            (ray_start[0] + beam_width*perp[0],
                             ray_start[1] + beam_width*perp[1])])
    
    c.fill(ray_area.path(), [color.grey(0.7)])
    c.fill(covering.path(), [color.grey(0.4)])

    # Draw the ray
    c.stroke(path.line(ray_head[0], ray_head[1], ray_tail[0], ray_tail[1]),
             [style.linewidth.Thick, deco.earrow(), color.rgb.red])
    c.fill(path.circle(ray_center[0], ray_center[1], 0.1))
    c.text(ray_head[0]-cellsize*0.1, ray_head[1], r"$\gamma$")
    c.text(ray_center[0], ray_center[1], r"$dr$",
           [text.halign.right, text.valign.bottom])

    # Draw the beam boundaries
    c.stroke(path.line(ray_start[0] + 0.1*beam_width*perp[0],
                       ray_start[1] + 0.1*beam_width*perp[1],
                       ray_start[0] - 0.1*beam_width*perp[0],
                       ray_start[1] - 0.1*beam_width*perp[1]),
             [style.linewidth.thick])
    c.stroke(path.line(ray_end[0] + 0.1*beam_width*perp[0],
                       ray_end[1] + 0.1*beam_width*perp[1],
                       ray_end[0] - 0.1*beam_width*perp[0],
                       ray_end[1] - 0.1*beam_width*perp[1]),
             [style.linewidth.thick])
    c.stroke(path.line(ray_head[0] + beam_width*perp[0],
                       ray_head[1] + beam_width*perp[1],
                       ray_tail[0] + beam_width*perp[0],
                       ray_tail[1] + beam_width*perp[1]),
             [style.linestyle.dashed])
    c.stroke(path.line(ray_head[0] - beam_width*perp[0],
                       ray_head[1] - beam_width*perp[1],
                       ray_tail[0] - beam_width*perp[0],
                       ray_tail[1] - beam_width*perp[1]),
             [style.linestyle.dashed])

    # Cell and center
    sqrt2inv = math.sqrt(0.5)
    c.stroke(path.rect(origin[0], origin[1], origin[0]+cellsize,
                       origin[1]+cellsize))
    c.stroke(path.line(origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize))
    c.stroke(path.line(origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize))

    c.text(0,1.05*cellsize,'(a)')
    c.writeEPSfile("covering")
    del c

    return

#------------------------------------------------------------------------
# Figure to illustrate the estimation of the correction factor
#------------------------------------------------------------------------

def figure2():
    c = canvas.canvas()

    # Cell and center
    sqrt2inv = math.sqrt(0.5)
    c.stroke(path.rect(origin[0], origin[1], origin[0]+cellsize,
                       origin[1]+cellsize))
    c.stroke(path.line(origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize))
    c.stroke(path.line(origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize-sqrt2inv*0.05*cellsize,
                       origin[0]+0.5*cellsize+sqrt2inv*0.05*cellsize))
    
    # Ray properties
    ray_start = np.array([origin[0] + ray_enter*cellsize, origin[1]])
    dx = origin[1] + cellsize - ray_start[0]
    dy = dx * normal[1] / normal[0]
    ray_end = np.array( [origin[0] + cellsize, origin[1]+dy] )
    ray_center = 0.5*(ray_start + ray_end)
    
    exterior = 2.0
    ray_head = ray_start - exterior*normal
    ray_tail = ray_end + exterior*normal
    
    # Draw the ray
    c.stroke(path.line(ray_head[0], ray_head[1], ray_tail[0], ray_tail[1]),
             [style.linewidth.Thick, deco.earrow(), color.rgb.red])
    c.fill(path.circle(ray_center[0], ray_center[1], 0.1))
    c.text(ray_head[0]-cellsize*0.1, ray_head[1], r"$\gamma$")

    # Draw the beam boundaries
    perp = np.array([-normal[1], normal[0]])
    beam_width = 0.4 * cellsize
    c.stroke(path.line(ray_head[0] + beam_width*perp[0],
                       ray_head[1] + beam_width*perp[1],
                       ray_tail[0] + beam_width*perp[0],
                       ray_tail[1] + beam_width*perp[1]),
             [style.linestyle.dashed])
    c.stroke(path.line(ray_head[0] - beam_width*perp[0],
                       ray_head[1] - beam_width*perp[1],
                       ray_tail[0] - beam_width*perp[0],
                       ray_tail[1] - beam_width*perp[1]),
             [style.linestyle.dashed])

    # Distances to the cell edges
    c.stroke(path.line(ray_center[0], ray_center[1],
                       ray_center[0], origin[1]),
             [color.rgb.blue])
#             [deco.earrow(), deco.barrow(), color.rgb.blue])
    c.stroke(path.line(ray_center[0], ray_center[1],
                       origin[0] + cellsize, ray_center[1]),
             [color.rgb.blue])
#             [deco.earrow(), deco.barrow(), color.rgb.blue])
    unit.set(xscale=0.2*cellsize)

    label0_pos = [ray_center[0]+0.02*cellsize,
                  0.5*(ray_center[1]+origin[1])-0.2*cellsize]
    label1_pos = [0.5*(ray_center[0] + origin[0] + cellsize) - 0.2*cellsize,
                  ray_center[1] + 0.2 * cellsize]
    
    c.text(label0_pos[0], label0_pos[1], r"$D_{c,0}$")
    c.text(label1_pos[0], label1_pos[1], r"$D_{c,1}$")
    unit.set(xscale=0.35*cellsize)

    line_center0 = [ray_center[0],
                    0.5*(ray_center[1] + origin[1])]
    line_center1 = [0.5*(ray_center[0] + origin[0] + cellsize),
                    ray_center[1]]

    # Account for text width
    label0_pos[0] += 0.07*cellsize
    label0_pos[1] += 0.05*cellsize
    label1_pos[0] += 0.13*cellsize
    label1_pos[1] += 0.02*cellsize
    line_center1[0] += 0.05*cellsize

    dx_path0 = label0_pos[0] - line_center0[0]
    dx_path1 = label1_pos[0] - line_center1[0]
    dy_path0 = -(label0_pos[1] - line_center0[1])
    dy_path1 = -(label1_pos[1] - line_center1[1])
    
    pp = path.curve(label0_pos[0], label0_pos[1],
                    label0_pos[0] + 0.2*dx_path0,
                    label0_pos[1] + 0.5*dy_path0,
                    label0_pos[0] + 0.5*dx_path0,
                    label0_pos[1] + 0.7*dy_path0,
                    line_center0[0], line_center0[1])
    c.stroke(pp, [deco.earrow()])

    pp = path.curve(label1_pos[0], label1_pos[1],
                    label1_pos[0] - 0.5*dx_path1,
                    label1_pos[1] + 0.2*dy_path1,
                    label1_pos[0] - 1*dx_path1,
                    label1_pos[1] + 0.7*dy_path1,
                    line_center1[0], line_center1[1])
    c.stroke(pp, [deco.earrow()])

    # Beam width
    c.stroke(path.line(ray_head[0] + beam_width*perp[0] + normal[0],
                       ray_head[1] + beam_width*perp[1] + normal[1],
                       ray_head[0] - beam_width*perp[0] + normal[0],
                       ray_head[1] - beam_width*perp[1] + normal[1]),
             [style.linewidth.Thick, deco.earrow(), deco.barrow(), color.rgb.blue])

    lpix_pos = [ray_head[0] + beam_width*perp[0] + normal[0],
                ray_head[1] + beam_width*perp[1] + normal[1] + 0.4*cellsize]
    c.text(lpix_pos[0], lpix_pos[1], r"$L_{pix}$")

    lpix_pos[0] += 0.15*cellsize
    lpix_pos[1] -= 0.05*cellsize
    beam_point = [0.5*(ray_head[0] + beam_width*perp[0] + normal[0] +
                       ray_head[0] - 0.2*beam_width*perp[0] + normal[0]),
                  0.5*(ray_head[1] + beam_width*perp[1] + normal[1] +
                       ray_head[1] - 0.2*beam_width*perp[1] + normal[1])]
    dx_path = [lpix_pos[0]-beam_point[0], lpix_pos[1]-beam_point[1]]
    
    pp = path.curve(lpix_pos[0], lpix_pos[1],
                    lpix_pos[0] + 0.5*dx_path[0],
                    lpix_pos[1] - 0.2*dx_path[1],
                    lpix_pos[0] + 0.7*dx_path[0],
                    lpix_pos[1] - 0.7*dx_path[1],
                    beam_point[0], beam_point[1])
    c.stroke(pp, [deco.earrow()])

    c.text(0,1.05*cellsize,'(b)')
    c.writeEPSfile("calc_fc")
    del c

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------

figure1()
figure2()
