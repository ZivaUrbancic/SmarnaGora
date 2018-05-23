from __future__ import division
import math
#import random
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch


class Trikotnik:

    def __init__(self, triterka):
        self.oglisca = triterka
        a, b, c = triterka
        self.sosedi = []
        self.counterclockwise = oriented_counterclockwise(triterka)

    def doloci_sosede(self, T):
        self.sosedi = []
        for trikotnik in T:
            if self != trikotnik:
                if imata_skupni_rob(self, trikotnik):
                    self.sosedi.append(trikotnik)


def imata_skupni_rob(jaz, trikotnik):
    # Function returns weather given triangles have a common edge (are neighbours).
    if len(set(jaz.oglisca).intersection(set(trikotnik.oglisca))) == 2:
        return True
    return False


def informacija_za_menjavo(jaz, trikotnik):
    # Given two triangles this function returns points a, b, c and d where a and b lie on th common edge, and c and d
    # are the remaining vertices.
    presek = set(jaz.oglisca).intersection(set(trikotnik.oglisca))
    if len(presek) == 2:
        a = presek.pop()
        b = presek.pop()
        c = (set(jaz.oglisca).difference(trikotnik.oglisca)).pop()
        d = (set(trikotnik.oglisca).difference(jaz.oglisca)).pop()
        return a, b, c, d
    else:
        return None


def points(triangulacija):
    tocke = set()
    for triangle in triangulacija:
        tocke = tocke.union(set(triangle))
    return tocke


def generify(S):
    # This function adds murmur to our input.
    noviS = []
    for trikotnik in S:
        novi = []
        for i in range(3):
            a, b = trikotnik[i]
            a = a * (1 + math.sin(1 / np.random.randint(1, high=200)))
            b = b * (1 + math.sin(1 / np.random.randint(1, high=200)))
            novi.append((a, b))
        noviS.append(tuple(novi))
    return noviS


def switch(trikotnik, sosed, T):
    # This function takes two triangles and a triangulation as input and it returns a triangulation we get if we flip
    #  the common edge
    if sosed not in T:
        return T
    novi_T = [t for t in T]
    if not pogoj_za_menjavo(trikotnik, sosed):
        return T
    a, b, tocka1, tocka2 = informacija_za_menjavo(trikotnik, sosed)
    novi_T.remove(trikotnik)
    novi_T.remove(sosed)
    dodani1 = Trikotnik((a, tocka1, tocka2))
    dodani2 = Trikotnik((b, tocka1, tocka2))
    novi_T.append(dodani1)
    novi_T.append(dodani2)
    for t in sosed.sosedi:
        t.doloci_sosede(novi_T)
    for t in trikotnik.sosedi:
        t.doloci_sosede(novi_T)
    dodani1.doloci_sosede(novi_T)
    dodani2.doloci_sosede(novi_T)
    return novi_T


def oriented_counterclockwise(trikotnik):
    # Checks weather the given triangle is counter clockwise oriented
    (x1, y1), (x2, y2), (x3, y3) = trikotnik
    return ((x2 - x1) * (y2 + y1) + (x3 - x2) * (y3 + y2) + (x1 - x3) * (y1 + y3)) < 0


def pogoj_za_menjavo(trikotnik, sosed):
    # Checks if we should perform edge-flip on two given triangles
    if sosed not in trikotnik.sosedi:
        return False
    else:
        a, b, c, d = informacija_za_menjavo(trikotnik, sosed)
        if incircle((a, b, c), d):
            return True
        else:
            return False


def incircle(triangle, tocka):
    # Checks weather tocka lies in the triangle
    if oriented_counterclockwise(triangle):
        (a1, a2), (b1, b2), (c1, c2) = triangle
    else:
        (a1, a2), (c1, c2), (b1, b2) = triangle
    a, b = tocka
    matrika = [[a1, a2, a1**2 + a2**2, 1], [b1, b2, b1**2 + b2**2, 1], [c1, c2, c1**2 + c2**2, 1],
                   [a, b, a**2 + b**2, 1]]
    if np.linalg.det(matrika) > 0:
        return True
    return False


def is_optimal(neko_ime):
    # Checks if triangulation is Delaunay and it returns triangles to perform edge-flip on if it is not.
    for trikotnik in neko_ime:
        for sosed in trikotnik.sosedi:
            if sosed in neko_ime and pogoj_za_menjavo(trikotnik, sosed):
                return trikotnik, sosed
    return None


def optimize(T_):
    # Main function. It takes a triangulation and returns it optimized version - Delaunay triangulation.
    triangulacija_s_trikotniki = []
    for trikotnik in T_:
        triangulacija_s_trikotniki.append(Trikotnik(trikotnik))
    for trikotnik in triangulacija_s_trikotniki:
        trikotnik.doloci_sosede(triangulacija_s_trikotniki)
    par = is_optimal(triangulacija_s_trikotniki)
    while par != None:
        trikotnik, sosed = par
        triangulacija_s_trikotniki = [t for t in switch(trikotnik, sosed, triangulacija_s_trikotniki)]
        par = is_optimal(triangulacija_s_trikotniki)
    triangulacija = [t.oglisca for t in triangulacija_s_trikotniki]
    narisi(T_, triangulacija)
    return triangulacija


def narisi(T_, T):
    # Function I use to plot triangulations.
    tocke = points(T)
    xi = [x for (x, _) in tocke]
    yi = [y for (_, y) in tocke]
    Mx = math.ceil(max(xi))
    mx = math.floor(min(xi))
    My = math.ceil(max(yi))
    my = math.floor(min(yi))
    fig = plt.figure()
    bx = fig.add_subplot(121)
    xrange = [mx - 1, Mx + 1]
    yrange = [my - 1, My + 1]
    bx.set_title("Initial triangulation")
    bx.set_xlim(*xrange)
    bx.set_ylim(*yrange)
    bx.set_xticks(range(*xrange + [xrange[-1]]))
    bx.set_yticks(range(*yrange + [yrange[-1]]))
    bx.set_aspect(1)
    for trikotnik in T_:
        t_i = Polygon(trikotnik)
        patch_i = PolygonPatch(t_i, color=np.random.rand(3, ))
        bx.add_patch(patch_i)
    ax = fig.add_subplot(122)
    ax.set_title("Delauney triangulation")
    ax.set_xlim(*xrange)
    ax.set_ylim(*yrange)
    ax.set_xticks(range(*xrange + [xrange[-1]]))
    ax.set_yticks(range(*yrange + [yrange[-1]]))
    ax.set_aspect(1)
    for trikotnik in T:
        t_i = Polygon(trikotnik)
        patch_i = PolygonPatch(t_i, color=np.random.rand(3, ))
        ax.add_patch(patch_i)
    plt.show()


#test1 = [((0, 0), (5, -1), (7, -5)), ((5, -1), (7, -5), (9, 4)), ((0, 0), (5, -1), (9, 4)), ((0, 0), (3, 9), (9, 4))]
#test3 = [[(-122, -225), (209, -219), (7, -205)], [(-122, -225), (7, -205), (-14, -203)],
#          [(209, -219), (7, -205), (-14, -203)], [(-122, -225), (-14, -203), (169, -163)],
#          [(209, -219), (-14, -203), (169, -163)], [(-122, -225), (169, -163), (-107, -35)],
#          [(209, -219), (169, -163), (131, -30)], [(169, -163), (-107, -35), (131, -30)],
#          [(-107, -35), (131, -30), (98, 8)], [(209, -219), (131, -30), (182, 140)], [(-107, -35), (98, 8), (182, 140)],
#          [(131, -30), (98, 8), (182, 140)], [(-122, -225), (-107, -35), (-202, 150)],
#          [(-107, -35), (182, 140), (-202, 150)], [(182, 140), (-202, 150), (5, 156)],
#          [(-202, 150), (5, 156), (-134, 157)], [(182, 140), (5, 156), (169, 158)], [(5, 156), (-134, 157), (169, 158)],
#          [(-202, 150), (-134, 157), (-14, 173)], [(-134, 157), (169, 158), (-14, 173)],
#          [(-202, 150), (-14, 173), (3, 185)], [(169, 158), (-14, 173), (3, 185)]]
# test5 = [[(-1046, 139), (-1019, 337), (-1014, -957)], [(-1046, 139), (-1019, 337), (-1004, 701)],
#          [(-1019, 337), (-1014, -957), (-1004, 701)], [(-1014, -957), (-1004, 701), (-982, 978)],
#          [(-1014, -957), (-982, 978), (-974, 492)], [(-1014, -957), (-974, 492), (-905, -307)],
#          [(-982, 978), (-974, 492), (-905, -307)], [(-982, 978), (-905, -307), (-891, 852)],
#          [(-1014, -957), (-905, -307), (-837, -732)], [(-905, -307), (-891, 852), (-837, -732)],
#          [(-891, 852), (-837, -732), (-812, 431)], [(-1014, -957), (-837, -732), (-771, -750)],
#          [(-837, -732), (-812, 431), (-771, -750)], [(-812, 431), (-771, -750), (-770, -124)],
#          [(-771, -750), (-770, -124), (-766, -241)], [(-812, 431), (-770, -124), (-750, -154)],
#          [(-771, -750), (-766, -241), (-750, -154)], [(-770, -124), (-766, -241), (-750, -154)],
#          [(-982, 978), (-891, 852), (-739, 822)], [(-891, 852), (-812, 431), (-739, 822)],
#          [(-812, 431), (-750, -154), (-739, 822)], [(-771, -750), (-750, -154), (-714, 794)],
#          [(-750, -154), (-739, 822), (-714, 794)], [(-771, -750), (-714, 794), (-695, -428)],
#          [(-714, 794), (-695, -428), (-691, 682)], [(-695, -428), (-691, 682), (-689, -202)],
#          [(-695, -428), (-689, -202), (-673, 88)], [(-691, 682), (-689, -202), (-673, 88)],
#          [(-771, -750), (-695, -428), (-672, -369)], [(-695, -428), (-673, 88), (-672, -369)],
#          [(-1014, -957), (-771, -750), (-671, -1046)], [(-771, -750), (-672, -369), (-671, -1046)],
#          [(-691, 682), (-673, 88), (-660, 510)], [(-673, 88), (-672, -369), (-660, 510)],
#          [(-672, -369), (-671, -1046), (-660, 510)], [(-671, -1046), (-660, 510), (-652, -537)],
#          [(-660, 510), (-652, -537), (-650, 91)], [(-671, -1046), (-652, -537), (-610, 78)],
#          [(-660, 510), (-650, 91), (-610, 78)], [(-652, -537), (-650, 91), (-610, 78)],
#          [(-982, 978), (-739, 822), (-540, 809)], [(-739, 822), (-714, 794), (-540, 809)],
#          [(-714, 794), (-691, 682), (-540, 809)], [(-691, 682), (-660, 510), (-540, 809)],
#          [(-671, -1046), (-610, 78), (-540, 809)], [(-660, 510), (-610, 78), (-540, 809)],
#          [(-671, -1046), (-540, 809), (-537, 12)], [(-540, 809), (-537, 12), (-514, 592)],
#          [(-671, -1046), (-537, 12), (-407, -160)], [(-540, 809), (-514, 592), (-407, -160)],
#          [(-537, 12), (-514, 592), (-407, -160)], [(-540, 809), (-407, -160), (-382, 622)],
#          [(-671, -1046), (-407, -160), (-378, -682)], [(-407, -160), (-382, 622), (-378, -682)],
#          [(-382, 622), (-378, -682), (-348, -310)], [(-671, -1046), (-378, -682), (-307, -791)],
#          [(-382, 622), (-348, -310), (-307, -791)], [(-378, -682), (-348, -310), (-307, -791)],
#          [(-382, 622), (-307, -791), (-277, -76)], [(-540, 809), (-382, 622), (-275, 550)],
#          [(-382, 622), (-277, -76), (-275, 550)], [(-671, -1046), (-307, -791), (-274, -944)],
#          [(-307, -791), (-277, -76), (-274, -944)], [(-277, -76), (-275, 550), (-274, -944)],
#          [(-275, 550), (-274, -944), (-251, 161)], [(-274, -944), (-251, 161), (-247, -599)],
#          [(-251, 161), (-247, -599), (-236, -310)], [(-275, 550), (-251, 161), (-226, 309)],
#          [(-251, 161), (-236, -310), (-226, 309)], [(-247, -599), (-236, -310), (-216, 199)],
#          [(-236, -310), (-226, 309), (-216, 199)], [(-982, 978), (-540, 809), (-204, 759)],
#          [(-540, 809), (-275, 550), (-204, 759)], [(-275, 550), (-226, 309), (-204, 759)],
#          [(-226, 309), (-216, 199), (-204, 759)], [(-982, 978), (-204, 759), (-188, 1036)],
#          [(-216, 199), (-204, 759), (-188, 1036)], [(-247, -599), (-216, 199), (-168, 716)],
#          [(-216, 199), (-188, 1036), (-168, 716)], [(-274, -944), (-247, -599), (-167, -198)],
#          [(-247, -599), (-168, 716), (-167, -198)], [(-274, -944), (-167, -198), (-161, -393)],
#          [(-168, 716), (-167, -198), (-161, -393)], [(-168, 716), (-161, -393), (-141, -80)],
#          [(-188, 1036), (-168, 716), (-114, -17)], [(-168, 716), (-141, -80), (-114, -17)],
#          [(-161, -393), (-141, -80), (-114, -17)], [(-274, -944), (-161, -393), (-114, -798)],
#          [(-161, -393), (-114, -17), (-114, -798)], [(-188, 1036), (-114, -17), (-109, 295)],
#          [(-114, -17), (-114, -798), (-109, 295)], [(-114, -798), (-109, 295), (-105, -173)],
#          [(-188, 1036), (-109, 295), (-81, 874)], [(-114, -798), (-105, -173), (-81, 874)],
#          [(-109, 295), (-105, -173), (-81, 874)], [(-114, -798), (-81, 874), (-31, -384)],
#          [(-81, 874), (-31, -384), (-30, 709)], [(-31, -384), (-30, 709), (-29, -250)],
#          [(-114, -798), (-31, -384), (12, -418)], [(-31, -384), (-29, -250), (12, -418)],
#          [(-30, 709), (-29, -250), (12, -418)], [(-188, 1036), (-81, 874), (26, 808)],
#          [(-81, 874), (-30, 709), (26, 808)], [(-30, 709), (12, -418), (26, 808)], [(-188, 1036), (26, 808), (26, 956)],
#          [(12, -418), (26, 808), (41, 512)], [(26, 808), (26, 956), (41, 512)], [(-274, -944), (-114, -798), (42, -703)],
#          [(-114, -798), (12, -418), (42, -703)], [(12, -418), (41, 512), (42, -703)], [(41, 512), (42, -703), (59, -205)],
#          [(26, 956), (41, 512), (84, -681)], [(-274, -944), (42, -703), (84, -681)], [(41, 512), (59, -205), (84, -681)],
#          [(42, -703), (59, -205), (84, -681)], [(-274, -944), (84, -681), (108, -709)],
#          [(26, 956), (84, -681), (108, -709)], [(26, 956), (108, -709), (148, 244)],
#          [(108, -709), (148, 244), (152, -137)], [(108, -709), (152, -137), (165, -628)],
#          [(148, 244), (152, -137), (165, -628)], [(148, 244), (165, -628), (207, -430)],
#          [(-274, -944), (108, -709), (218, -747)], [(108, -709), (165, -628), (218, -747)],
#          [(165, -628), (207, -430), (218, -747)], [(26, 956), (148, 244), (226, 364)],
#          [(148, 244), (207, -430), (226, 364)], [(207, -430), (218, -747), (226, 364)],
#          [(-671, -1046), (-274, -944), (227, -906)], [(-274, -944), (218, -747), (227, -906)],
#          [(218, -747), (226, 364), (227, -906)], [(-671, -1046), (227, -906), (293, -938)],
#          [(226, 364), (227, -906), (293, -938)], [(26, 956), (226, 364), (294, 543)],
#          [(226, 364), (293, -938), (294, 543)], [(26, 956), (294, 543), (314, 674)],
#          [(293, -938), (294, 543), (314, 674)], [(293, -938), (314, 674), (350, 501)],
#          [(293, -938), (350, 501), (371, -46)], [(293, -938), (371, -46), (376, -728)],
#          [(26, 956), (314, 674), (407, 780)], [(314, 674), (350, 501), (407, 780)],
#          [(350, 501), (371, -46), (407, 780)], [(371, -46), (376, -728), (407, 780)],
#          [(293, -938), (376, -728), (462, -655)], [(376, -728), (407, 780), (462, -655)],
#          [(293, -938), (462, -655), (504, -694)], [(407, 780), (462, -655), (504, -694)],
#          [(407, 780), (504, -694), (515, 210)], [(504, -694), (515, 210), (525, -18)],
#          [(-188, 1036), (26, 956), (531, 889)], [(26, 956), (407, 780), (531, 889)], [(407, 780), (515, 210), (531, 889)],
#          [(515, 210), (525, -18), (531, 889)], [(525, -18), (531, 889), (534, 672)], [(525, -18), (534, 672), (540, 562)],
#          [(531, 889), (534, 672), (540, 562)], [(293, -938), (504, -694), (549, -857)],
#          [(504, -694), (525, -18), (549, -857)], [(525, -18), (540, 562), (549, -857)],
#          [(531, 889), (540, 562), (577, -151)], [(540, 562), (549, -857), (577, -151)],
#          [(-671, -1046), (293, -938), (579, -1034)], [(293, -938), (549, -857), (579, -1034)],
#          [(549, -857), (577, -151), (579, -1034)], [(577, -151), (579, -1034), (580, -740)],
#          [(531, 889), (577, -151), (595, -240)], [(577, -151), (580, -740), (595, -240)],
#          [(579, -1034), (580, -740), (595, -240)], [(-188, 1036), (531, 889), (603, 1019)],
#          [(531, 889), (595, -240), (603, 1019)], [(579, -1034), (595, -240), (623, -723)],
#          [(595, -240), (603, 1019), (623, -723)], [(603, 1019), (623, -723), (747, 470)],
#          [(579, -1034), (623, -723), (751, -997)], [(623, -723), (747, 470), (751, -997)],
#          [(747, 470), (751, -997), (796, 125)], [(751, -997), (796, 125), (796, -938)],
#          [(751, -997), (796, -938), (807, -941)], [(796, 125), (796, -938), (807, -941)],
#          [(603, 1019), (747, 470), (857, 131)], [(747, 470), (796, 125), (857, 131)], [(796, 125), (807, -941), (857, 131)],
#          [(807, -941), (857, 131), (862, -812)], [(603, 1019), (857, 131), (882, 439)],
#          [(857, 131), (862, -812), (882, 439)], [(862, -812), (882, 439), (941, -264)],
#          [(862, -812), (941, -264), (944, -538)], [(807, -941), (862, -812), (1003, -572)],
#          [(882, 439), (941, -264), (1003, -572)], [(862, -812), (944, -538), (1003, -572)],
#          [(941, -264), (944, -538), (1003, -572)], [(882, 439), (1003, -572), (1005, -310)],
#          [(751, -997), (807, -941), (1009, -760)], [(807, -941), (1003, -572), (1009, -760)],
#          [(1003, -572), (1005, -310), (1009, -760)]]
#print(optimize(test1))
