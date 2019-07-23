# Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
# Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy

# Specification of data:
# cicle: [x, y, r]
# circlesegment: [x, y, r, plow, phigh], plow/phigh counter-clockwise
# line: [x1, y1, x2, y2 ]

# Glyph data
def glyph(char):
    pi = numpy.pi;
    circles = []
    segments = []
    lines = []
    if char == 'a':
        lines = [ [ -.5,-.8,.2,.75], [.2,.75,.5,-.7], [-.25,.15,.25,-.15] ]
    elif char == 'b':
        lines =  [ [-.75,.75,.4,.6], [.4,.6,.6,.2], [.6,.2,.2,0.], [.2,0.,.5,-.5],[.5,-.5,.15,-.8],[.15,-.8,.75,-.6],[-.6,-.7,-.4,.7] ]
    elif char == 'c':
        lines = [ [.5,.5,.0,.8], [.0,.8,-.7,.3], [-.7,.3,-.5,-.6], [-.5,-.6,.2,-.7], [.2,-.7,.6,-.6] ]
    elif char == 'd':
        lines = [ [-.75,.75,.4,.6], [.4,.6,.7,-.1], [.7,-.1,.15,-.8], [.15,-.8,-.75,-.6],[-.6,-.7,-.4,.7] ]
    elif char == 'e':
        lines = [ [.4,.7,-.7,.6],[-.7,.6,-.4,-.8],[-.4,-.8,.6,-.6],[-.65,0.,.5,0.] ]
    elif char == 'f':
        lines = [ [.4,.7,-.7,.6],[-.7,.6,-.4,-.8],[-.65,0.,.5,0.] ]
    elif char == 'g':
        lines = [ [.5,.5,.0,.8], [.0,.8,-.7,.3], [-.7,.3,-.5,-.6], [-.5,-.6,.2,-.7], [.2,-.7,.6,-.6], [.6,-.6,.5,0.], [.5,0.,0.,0.] ]
    elif char == 'h':
        lines = [ [-.4,.6,-.6,-.7],[.5,.6,.25,-.8],[-.6,.15,.6,-.15] ]
    elif char == 'i':
        lines = [ [-.15,.75,.15,-.75] ]
    elif char == 'j':
        lines = [ [-.5,.75,.4,.5],[.4,.5,.3,-.6],[.3,-.6,-.3,-.8],[-.3,-.8,-.5,-.4] ]
    elif char == 'k':
        lines = [ [-.6,.6,-.2,-.7],[-.5,-.2,.4,.6],[-.2,.15,.8,-.7] ]
    elif char == 'l':
        lines = [ [-.3,.8,-.7,-.6],[-.7,-.6,.4,-.9] ]
    elif char == 'm':
        lines = [ [-.8,-.8,-.4,.8],[-.7,.7,0.,-.25],[-.2,-.35,.4,.7],[.25,.7,.8,-.9] ]
    elif char == 'n':
        lines = [ [-.8,-.8,-.4,.8],[-.7,.7,.5,-.8],[.3,-.7,.9,.9] ]
    elif char == 'o':
        lines = [ [.5,.5,.0,.8], [.0,.8,-.7,.3], [-.7,.3,-.5,-.6], [-.5,-.6,.2,-.7], [.2,-.7,.6,-.6], [.5,.5,.6,-.6] ]
    elif char == 'p':
        lines = [ [-.75,.75,.4,.6], [.4,.6,.6,.2], [.6,.2,.2,0.], [-.6,-.7,-.4,.7], [.2,0.,-.8,.1] ]
    elif char == 'q':
        lines = [ [.5,.5,.0,.8], [.0,.8,-.7,.3], [-.7,.3,-.5,-.6], [-.5,-.6,.2,-.7], [.2,-.7,.6,-.6], [.5,.5,.6,-.6], [.3,-.3,.4,-.9] ]
    elif char == 'r':
        lines = [ [-.75,.75,.4,.6], [.4,.6,.6,.2], [.6,.2,.2,0.], [-.6,-.7,-.4,.7], [.2,0.,-.8,.1], [-.2,.2,.7,-.8] ]
    elif char == 's':
        lines = [ [.4,.7,-.2,.7],[-.2,.7,-.5,.25],[-.5,.25,-.15,0.],[-.15,0.,.6,-.2],[.6,-.2,.7,-.75],[.7,-.75,0.,-.9],[0.,-.9,-.8,-.7] ]
    elif char == 't':
        lines = [ [ -.7,.65,.75,.85 ], [0.,.85,.2,-.8] ]
    elif char == 'u':
        lines = [ [-.3,.6,-.6,-.5],[-.6,-.5,.15,-.85],[.15,-.85,.7,-.5],[.7,-.5,.4,.6] ]
    elif char == 'v':
        lines = [ [-.8,.5,0.,-.7],[0.,-.7,.6,.9] ]
    elif char == 'w':
        lines = [ [-.8,.8,-.4,-.8],[-.7,-.7,0.,.25],[-.2,.35,.4,-.7],[.25,-.7,.8,.9] ]
    elif char == 'x':
        lines = [ [-.6,.6,.5,-.7], [-.5,-.8,.9,.75] ]
    elif char == 'y':
        lines = [ [-.6,-.9,.6,.6],[-.5,.75,.1,-.2] ]
    elif char == 'z':
        lines = [ [ -.65,.7,.5,.4 ],[.5,.4,-.65,-.7],[-.65,-.7,.75,-.6] ]
    elif char == '0':
        lines = [ [.5,.5,.0,.8], [.0,.8,-.7,.3], [-.7,.3,-.5,-.6], [-.5,-.6,.2,-.7], [.2,-.7,.6,-.6], [.5,.5,.6,-.6] ]
    elif char == '1':
        lines = [ [-.5,.5,.6,.8], [.6,.8,.3,-.8] ]        
    elif char == '2':
        lines = [ [-.4,.7,.2,.7],[.2,.7,.5,.25],[.5,.25,.15,0.],[.15,0.,-.6,-.2],[-.6,-.2,-.7,-.75],[-.7,-.75,0.,-.9],[0.,-.9,-.8,-.7] ]
    elif char == '3':
        lines =  [ [-.75,.75,.4,.6], [.4,.6,.6,.2], [.6,.2,.2,0.], [.2,0.,.5,-.5],[.5,-.5,.15,-.8],[.15,-.8,.75,-.6] ]
    elif char == '4':
        lines =  [ [-.4,.65,-.6,-.5], [-.6,-.5,.7,-.35],[.4,-.8,-.3,.8] ]
    elif char == '5':
        lines =  [ [.5,.8,-.5,.65],[-.5,.65,-.6,-.3],[-.6,.3,.2,0.],[.2,0.,.65,-.4],[.65,-.4,.4,-.9],[.4,-.9,-.6,-.8] ]
    elif char == '6':
        lines = [ [ .5,.8,-.5,.4],[-.5,.4,-.6,-.8], [-.6,.3,.2,0.],[.2,0.,.65,-.4],[.65,-.4,.4,-.9],[.4,-.9,-.6,-.8]]
    elif char == '7':
        lines = [ [-.8,.9,.8,.6],[.8,.6,.2,-.9],[0.,.15,.9,0.] ]
    elif char == '8':
        lines = [ [-.4,.7,.4,.6],[.4,.6,.6,.3],[.6,.3,.2,0.],[.2,0.,.6,-.3],[.6,-.3,.55,-.8],[.55,-.8,0.,-.9],[0.,-.9,-.8,-.8],[-.8,-.8,-.9,-.3],[-.9,-.3,0.,0.],[0.,0.,-.6,.4],[-.6,.4,-.4,.7] ]
    elif char == '9':
        lines = [ [ -.5,-.8,.5,-.4],[.5,-.4,.6,.8], [.6,-.3,-.2,-0.],[-.2,-0.,-.65,.4],[-.65,.4,-.4,.9],[-.4,.9,.6,.8]]
    elif char == '.':
        circles = [ [ 0., -1., 0. ] ]
    elif char == ':':
        circles = [ [ -1., -.5, 0. ], [-1., .5, 0.] ]
    elif char == ',':
        lines = [ [ -1., -1., -2./3., -2./3. ] ]
    elif char == '-':
        lines = [ [ -2./3., 0., 2./3., 0. ] ]
    elif char == '+':
        lines = [ [ -2./3., 0., 2./3., 0. ], [ 0., -2./3., 0., 2./3.] ]
    elif char == '?':
        lines = [ [-.75,.5,0.,.8],[0.,.8,.6,.6],[.6,.6,.7,0.],[.7,0.,-.15,-.25],[-.15,-.25,0.,-.75], [.15,-1.,.15,-.95] ]
    elif char == '!':
        lines = [ [-.1,-.4,.15,.5], [-.1,-1.,-.1,-.95] ]
    elif char == '/':
        lines = [ [ -1., -1., 1., 1. ] ]
        
    return [ circles, segments, lines ]
        
# Length of glyph data in texture (for building index with offsets)
def pack_length(char):
    gly = glyph(char)
    circles = gly[0]
    segments = gly[1]
    lines = gly[2]    
    return 3 + len(circles)*3 + len(segments)*5 + len(lines)*4
