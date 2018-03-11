#!/usr/bin/env python3
import numpy as np
import cv2


H = 60
W = 2000
with open('ranks/index.html', 'w') as index:
    index.write('<HTML><BODY>\n')
    for l in open('interest', 'r'):
        gene, probe = l.strip().split('\t')

        image = np.zeros((H, W), dtype=np.uint8)
        with open("ranks/%s_%s" % (gene, probe), 'r') as f:
            for v in f:
                v = int(round(W * float(v.strip())))
                cv2.line(image, (v, 0), (v, H), 255, 1)
                pass
            pass
        cv2.imwrite("ranks/%s_%s.png" % (gene, probe), image)
        index.write("<h2>%s - %s</h2>\n" % (gene, probe))
        index.write("<img src='%s_%s.png'></img><br/><hr/>\n" % (gene,probe))
        pass
    index.write("</BODY></HTML>\n")
    pass



