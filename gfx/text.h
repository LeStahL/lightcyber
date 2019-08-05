/* Generated with shader-compressor by NR4/Team210. */
#ifndef TEXT_H
#define TEXT_H
const char * text_frag =
"/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" *\n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" *\n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" *\n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
" */\n"
" \n"
"#version 130\n"
"\n"
"uniform float iFontWidth, iTime;\n"
"uniform vec2 iResolution;\n"
"uniform sampler2D iChannel0, iFont;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"float a; // Aspect ratio\n"
"\n"
"// Hash function\n"
"void rand(in vec2 x, out float num)\n"
"{\n"
"    num = fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"// Arbitrary-frequency 2D noise\n"
"void lfnoise(in vec2 t, out float num)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise\n"
"    t = smoothstep(c.yy, c.xx, t); // TODO: add this for faster value noise\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    num = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\n"
"// Read short value from texture at index off\n"
"void rshort(in float off, out float val)\n"
"{\n"
"    // Parity of offset determines which byte is required.\n"
"    float hilo = mod(off, 2.);\n"
"    // Find the pixel offset your data is in (2 unsigned shorts per pixel).\n"
"    off *= .5;\n"
"    // - Determine texture coordinates.\n"
"    //     offset = i*iFontWidth+j for (i,j) in [0,iFontWidth]^2\n"
"    //     floor(offset/iFontWidth) = floor((i*iFontwidth+j)/iFontwidth)\n"
"    //                              = floor(i)+floor(j/iFontWidth) = i\n"
"    //     mod(offset, iFontWidth) = mod(i*iFontWidth + j, iFontWidth) = j\n"
"    // - For texture coordinates (i,j) has to be rescaled to [0,1].\n"
"    // - Also we need to add an extra small offset to the texture coordinate\n"
"    //   in order to always \"hit\" the right pixel. Pixel width is\n"
"    //     1./iFontWidth.\n"
"    //   Half of it is in the center of the pixel.\n"
"    vec2 ind = (vec2(mod(off, iFontWidth), floor(off/iFontWidth))+.05)/iFontWidth;\n"
"    // Get 4 bytes of data from the texture\n"
"    vec4 block = texture(iFont, ind);\n"
"    // Select the appropriate word\n"
"    vec2 data = mix(block.rg, block.ba, hilo);\n"
"    // Convert bytes to unsigned short. The lower bytes operate on 255,\n"
"    // the higher bytes operate on 65280, which is the maximum range \n"
"    // of 65535 minus the lower 255.\n"
"    val = round(dot(vec2(255., 65280.), data));\n"
"}\n"
"\n"
"// Read float value from texture at index off\n"
"void rfloat(in float off, out float val)\n"
"{\n"
"    // Convert the bytes to unsigned short as first step.\n"
"    float d;\n"
"    rshort(off, d);\n"
"    \n"
"    // Convert bytes to IEEE 754 float16. That is\n"
"    // 1 sign bit, 5 bit exponent, 11 bit mantissa.\n"
"    // Also it has a weird conversion rule that is not evident at all.\n"
"    float sign = floor(d/32768.),\n"
"        exponent = floor(d/1024.-sign*32.),\n"
"        significand = d-sign*32768.-exponent*1024.;\n"
"\n"
"    // Return full float16\n"
"    if(exponent == 0.)\n"
"    {\n"
"        val = mix(1., -1., sign) * 5.960464477539063e-08 * significand;\n"
"    }\n"
"    else\n"
"    {\n"
"        val = mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);\n"
"    }\n"
"}\n"
"\n"
"\n"
"// 2D box\n"
"void box(in vec2 x, in vec2 b, out float dst)\n"
"{\n"
"    vec2 d = abs(x) - b;\n"
"    dst = length(max(d,c.yy)) + min(max(d.x,d.y),0.);\n"
"}\n"
"\n"
"// Distance to circle\n"
"void circle(in vec2 x, out float d)\n"
"{\n"
"    d = abs(length(x)-1.0);\n"
"}\n"
"\n"
"// Distance to line segment\n"
"void lineseg(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\n"
"// 2D rhomboid\n"
"void rhomboid(in vec2 x, in vec2 b, in float tilt, out float dst)\n"
"{\n"
"    x.x -= tilt/2./b.y*x.y;\n"
"    box(x,b,dst);\n"
"}\n"
"\n"
"// Distance to hexagon pattern\n"
"void dhexagonpattern(in vec2 p, out float d, out vec2 ind) \n"
"{\n"
"    vec2 q = vec2( p.x*1.2, p.y + p.x*0.6 );\n"
"    \n"
"    vec2 pi = floor(q);\n"
"    vec2 pf = fract(q);\n"
"\n"
"    float v = mod(pi.x + pi.y, 3.0);\n"
"\n"
"    float ca = step(1.,v);\n"
"    float cb = step(2.,v);\n"
"    vec2  ma = step(pf.xy,pf.yx);\n"
"    \n"
"    d = dot( ma, 1.0-pf.yx + ca*(pf.x+pf.y-1.0) + cb*(pf.yx-2.0*pf.xy) );\n"
"    ind = pi + ca - cb*ma;\n"
"    ind = vec2(ind.x/1.2, ind.y);\n"
"    ind = vec2(ind.x, ind.y-ind.x*.6);\n"
"}\n"
"\n"
"// Distance to circle segment\n"
"void circlesegment(in vec2 x, in float r, in float p0, in float p1, out float d)\n"
"{\n"
"    float p = atan(x.y, x.x);\n"
"    vec2 philo = vec2(max(p0, p1), min(p0, p1));\n"
"    if((p < philo.x && p > philo.y) || (p+2.*pi < philo.x && p+2.*pi > philo.y) || (p-2.*pi < philo.x && p-2.*pi > philo.y))\n"
"    {\n"
"        d = abs(length(x)-r);\n"
"        return;\n"
"    }\n"
"    d = min(\n"
"        length(x-r*vec2(cos(p0), sin(p0))),\n"
"        length(x-r*vec2(cos(p1), sin(p1)))\n"
"        );\n"
"}\n"
"\n"
"// Compute distance to stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0) - s;\n"
"}\n"
"\n"
"// Get glyph data from texture\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst)\n"
"{\n"
"    float dis;\n"
"    box(x, 2.*size*c.xx, dis);\n"
"    if(dis > 0.)\n"
"    {\n"
"        dst = dis+.5*size;\n"
"        return;\n"
"    }\n"
"\n"
"    // Find glyph offset in glyph index\n"
"    float nglyphs, offset = 0;\n"
"    rfloat(1., nglyphs);\n"
"        \n"
"    for(float i=0.; i<nglyphs; i+=1.)\n"
"    {\n"
"        float ord;\n"
"        rfloat(2.+2.*i, ord);\n"
"        ord = floor(ord);\n"
"        \n"
"        if(ord == ordinal)\n"
"        {\n"
"            rfloat(2.+2.*i+1., offset);\n"
"            offset = floor(offset);\n"
"            break;\n"
"        }\n"
"    }\n"
"    \n"
"    if(offset == 0.) \n"
"    {\n"
"        dst = 1.;\n"
"        return;\n"
"    }\n"
"    \n"
"    // Get distance from glyph data\n"
"    float d = 1.;\n"
"    \n"
"    // Lines\n"
"    float nlines;\n"
"    rfloat(offset, nlines);\n"
"    nlines = floor(nlines);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<nlines; i+=1.)\n"
"    {\n"
"        float x1;\n"
"        rfloat(offset, x1);\n"
"        offset += 1.;\n"
"        float y1;\n"
"        rfloat(offset, y1);\n"
"        offset += 1.;\n"
"        float x2;\n"
"        rfloat(offset, x2);\n"
"        offset += 1.;\n"
"        float y2;\n"
"        rfloat(offset, y2);\n"
"        offset += 1.;\n"
"        float da;\n"
"        lineseg(x, size*vec2(x1,y1), size*vec2(x2, y2), da);\n"
"        d = min(d,da);\n"
"    }\n"
"    \n"
"    // Circles\n"
"    float ncircles;\n"
"    rfloat(offset, ncircles);\n"
"    ncircles = floor(ncircles);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<ncircles; i+=1.)\n"
"    {\n"
"        float xc;\n"
"        rfloat(offset, xc);\n"
"        offset += 1.;\n"
"        float yc;\n"
"        rfloat(offset, yc);\n"
"        offset += 1.;\n"
"        float r;\n"
"        rfloat(offset, r);\n"
"        offset += 1.;\n"
"        float da;\n"
"        circle( (x-size*vec2(xc, yc))/size/r,da);\n"
"        d = min(d, da*size*r);\n"
"    }\n"
"    \n"
"    // Circle segments\n"
"    float nsegments;\n"
"    rfloat(offset, nsegments);\n"
"    nsegments = floor(nsegments);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<nsegments; i+=1.)\n"
"    {\n"
"        float xc;\n"
"        rfloat(offset, xc);\n"
"        offset += 1.;\n"
"        float yc;\n"
"        rfloat(offset, yc);\n"
"        offset += 1.;\n"
"        float r;\n"
"        rfloat(offset, r);\n"
"        offset += 1.;\n"
"        float phi0;\n"
"        rfloat(offset, phi0);\n"
"        offset += 1.;\n"
"        float phi1;\n"
"        rfloat(offset, phi1);\n"
"        offset += 1.;\n"
"        float da;\n"
"        circlesegment(x-size*vec2(xc,yc), size*r, phi0, phi1, da);\n"
"        d = min(d, da);\n"
"    }\n"
"    \n"
"    if(nlines+ncircles+nsegments == 0.)\n"
"        dst = dis;\n"
"    else dst = d;\n"
"}\n"
"\n"
"// Get distance to string from database\n"
"void dstring(in vec2 x, in float ordinal, in float size, out float dst)\n"
"{\n"
"    // Get string database offset\n"
"    float stroff0;\n"
"    rfloat(0., stroff0);\n"
"    stroff0 = floor(stroff0);\n"
"    \n"
"    // Return 1 if wrong ordinal is supplied\n"
"    float nstrings;\n"
"    rfloat(stroff0, nstrings);\n"
"    nstrings = floor(nstrings);\n"
"    if(ordinal >= nstrings)\n"
"    {\n"
"        dst = 1.;\n"
"        return;\n"
"    }\n"
"    \n"
"    // Get offset and length of string from string database index\n"
"    float stroff;\n"
"    rfloat(stroff0+1.+2.*ordinal, stroff);\n"
"    stroff = floor(stroff);\n"
"    float len;\n"
"    rfloat(stroff0+2.+2.*ordinal, len);\n"
"    len = floor(len);\n"
"    \n"
"    // Draw glyphs\n"
"    vec2 dx = mod(x-size, 2.*size)-size, \n"
"        ind = ceil((x-dx+size)/2./size);\n"
"    \n"
"    // Bounding box\n"
"    float bound;\n"
"    box(x-size*(len-3.)*c.xy, vec2(size*len, 1.*size), bound);\n"
"    if(bound > 0.)\n"
"    {\n"
"        dst = bound+.5*size;\n"
"        return;\n"
"    }\n"
"    \n"
"    float da;\n"
"    rfloat(stroff+ind.x, da);\n"
"    da = floor(da);\n"
"    dglyph(dx, da, .7*size, dst);\n"
"}\n"
"\n"
"// distance to a floating point number string\n"
"// for debugging stuff while shader is loaded\n"
"void dfloat(in vec2 x, in float num, in float size, out float dst)\n"
"{\n"
"    float d = 1., index = 0.;\n"
"    \n"
"    // Determine sign and output it if present\n"
"    float sign = sign(num), exp = 0.;\n"
"    if(sign<0.)\n"
"    {\n"
"        float da;\n"
"        dglyph(x, 45., .7*size, da);\n"
"        d = min(d, da);\n"
"        index += 1.;\n"
"        num *= -1.;\n"
"    }\n"
"    \n"
"    // The first power of ten that floors num to anything not zero is the exponent\n"
"    for(exp = -15.; exp < 15.; exp += 1.)\n"
"        if(floor(num*pow(10.,exp)) != 0.)\n"
"            break;\n"
"    exp *= -1.;\n"
"    // Determine the significand and output it\n"
"    for(float i = exp; i >= max(exp-5.,-33); i -= 1.)\n"
"    {\n"
"        float po = pow(10.,i);\n"
"        float ca = floor(num/po);\n"
"        num -= ca*po;\n"
"        \n"
"        float da;\n"
"        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"        d = min(d, da);\n"
"        index += 1.;\n"
"        if(i == exp) // decimal point\n"
"        {\n"
"            dglyph(x-2.*index*size*c.xy, 46., .7*size, da);\n"
"            d = min(d, da);\n"
"            index += 1.;\n"
"        }\n"
"    }\n"
"    \n"
"    // Output the exponent\n"
"    float db;\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 101., .7*size, db);\n"
"    d = min(d, db);\n"
"    index += 1.;\n"
"    if(exp < 0.) // Sign\n"
"    {\n"
"        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 45., .7*size,db);\n"
"        d = min(d, db);\n"
"        index += 1.;\n"
"        exp *= -1.;\n"
"    }\n"
"    float ca = floor(exp/10.);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, db);\n"
"    d = min(d, db);\n"
"    index += 1.;\n"
"    ca = floor(exp-10.*ca);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, db);\n"
"    d = min(d, db);\n"
"    index += 1.;\n"
"    \n"
"    dst = d;\n"
"}\n"
"\n"
"// Add scene contents\n"
"void add(in vec4 src1, in vec4 src2, out vec4 dst)\n"
"{\n"
"    dst = mix(src1, src2, smoothstep(0., 1.5/iResolution.y, -src2.x));\n"
"}\n"
"\n"
"void blendadd(in vec4 src1, in vec4 src2, in float tlo, in float thi, out vec4 dst)\n"
"{\n"
"    vec4 added;\n"
"    add(src1, src2, added);\n"
"    dst = mix(src1, added, smoothstep(tlo-.5,tlo+.5,iTime)*(1.-smoothstep(thi-.5,thi+.5,iTime)));\n"
"}\n"
"\n"
"// UI Window Control\n"
"void window(in vec2 x, in vec2 size, in vec3 bg, in float title_index, out vec4 col)\n"
"{\n"
"    size.x *= .5;\n"
"    col = vec4(1., bg);\n"
"    \n"
"    const float cellsize = .015, bordersize = .005;\n"
"    vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),\n"
"        bordercolor = vec3(1.00,0.71,0.02);\n"
"    vec4 c2 = vec4(1., titlecolor);\n"
"    \n"
"    float dhx, dhy;\n"
"    vec2 ind;\n"
"    dhexagonpattern(72.*x,  dhx, ind);\n"
"    stroke(dhx, .1, dhx);\n"
"    lfnoise(ind-iTime, dhy);\n"
"    \n"
"    // Window background\n"
"    box(x+.5*size*c.yx,size*vec2(1.,.5),c2.x);\n"
"    c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/size.y), .5+.5*dhy*step(0.,dhx));\n"
"    add(col, c2, col);\n"
"    \n"
"    // Title bar\n"
"    c2.gba = titlecolor;\n"
"    rhomboid(x+.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);\n"
"   	add(col, c2, col);\n"
"    rhomboid(x, vec2(.65*size.x,cellsize), cellsize, c2.x);\n"
"   	add(col, c2, col);\n"
"    rhomboid(x-.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);\n"
"   	add(col, c2, col);\n"
"    \n"
"    // Border of title bar\n"
"    c2 = vec4(1., bordercolor);\n"
"    stroke(col.x,bordersize,c2.x);\n"
"    add(col,c2,col);\n"
"    \n"
"    // Window Border\n"
"    lineseg(x, -.9*size.x*c.xy, -size.x*c.xy, c2.x);\n"
"    float d;\n"
"    lineseg(x, -size.x*c.xy, -size, d);\n"
"    c2.x = min(c2.x, d);\n"
"    lineseg(x, -size, size*c.xz, d);\n"
"    c2.x = min(c2.x, d);\n"
"    lineseg(x, size*c.xz, size*c.xy, d);\n"
"    c2.x = min(c2.x, d);\n"
"    lineseg(x, .9*size.x*c.xy, size.x*c.xy, d);\n"
"    c2.x = min(c2.x, d);\n"
"    stroke(c2.x,.25*bordersize,c2.x);\n"
"    add(col, c2, col);\n"
"}\n"
"\n"
"void progressbar(in vec2 x, in float width, in float progress, out vec4 col)\n"
"{\n"
"    const float cellsize = .015, bordersize = .005;\n"
"    vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),\n"
"        bordercolor = vec3(1.00,0.71,0.02), bg = c.yyy;\n"
"    vec4 c2 = vec4(1., titlecolor);\n"
"    \n"
"    // Window background\n"
"    box(x+.5*width*c.yx,width*c.xy,c2.x);\n"
"    c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/cellsize), .5);\n"
"    add(col, c2, col);\n"
"    \n"
"    // Bar background\n"
"    c2.gba = titlecolor;\n"
"    rhomboid(x, vec2(.5*width,cellsize), cellsize, c2.x);\n"
"   	add(col, c2, col);\n"
"    \n"
"    // Border\n"
"    c2.gba = bordercolor;\n"
"    stroke(c2.x,.5*bordersize,c2.x);\n"
"    add(col, c2, col);\n"
"    \n"
"    // Progress\n"
"    float wc = width/cellsize;\n"
"    x.x -= .5*x.y;\n"
"    vec2 y = vec2(mod(x.x, 1.2*cellsize)-.6*cellsize, x.y),\n"
"        index = (x-y)/.6/cellsize;\n"
"    if(abs(index.x) < .8*wc && -index.x > .8*wc*(1.-2.*progress))\n"
"    {\n"
"        box(y, vec2(.5*cellsize, .8*cellsize), c2.x);\n"
"        add(col, c2, col);\n"
"    }\n"
"}\n"
"\n"
"// Revision logo of width 1.\n"
"void drevision(in vec2 x, in float r, out float dst)\n"
"{\n"
"    float l = length(x),\n"
"        p = atan(x.y,x.x),\n"
"	    d = abs(l-r*.07)-.02, \n"
"        k1 = abs(l-r*.16)-.03,\n"
"        k2 = abs(l-r*.21)-.02, \n"
"        k3 = abs(l-r*.35)-.03,\n"
"        k4 = abs(l-r*.45)-.02;\n"
"    vec4 n1;\n"
"    lfnoise(1.*c.xx-3.*iTime, n1.x);\n"
"    lfnoise(2.*c.xx-2.4*iTime, n1.y);\n"
"    lfnoise(3.*c.xx-2.9*iTime, n1.z);\n"
"    lfnoise(4.*c.xx-3.1*iTime, n1.w);\n"
"    n1 = mix(n1,c.yyyy, clamp((iTime-24.)/2.,0.,1.));\n"
"    d = min(d, mix(d, abs(l-.11)-.03, step(p, -1.71)*step(-2.73, p)));\n"
"    d = min(d, mix(d, k1, step(p+n1.x, 3.08)*step(2.82,p)));\n"
"    d = min(d, mix(d, k1, step(p+n1.x, 1.47)*step(.81,p)));\n"
"    d = min(d, mix(d, k1, step(p+n1.x, -.43)*step(-1.19,p)));\n"
"    d = min(d, mix(d, k2, step(p+n1.y, -2.88)*step(-pi,p)));\n"
"    d = min(d, mix(d, k2, step(p+n1.y, pi)*step(2.38,p)));\n"
"    d = min(d, mix(d, k2, step(p+n1.y, 2.1)*step(.51,p)));\n"
"    d = min(d, mix(d, k2, step(p+n1.y, .3)*step(-1.6,p)));\n"
"    d = min(d, abs(l-.24)-.02);\n"
"    d = min(d, mix(d, k3, step(p+n1.z, -2.18)*step(-pi, p)));\n"
"    d = min(d, mix(d, k3, step(p+n1.z, -1.23)*step(-1.7, p)));\n"
"    d = min(d, mix(d, k3, step(p+n1.z, -.58)*step(-.78, p)));\n"
"    d = min(d, mix(d, k3, step(p+n1.z, 0.)*step(-.29, p)));\n"
"    d = min(d, mix(d, k3, step(p+n1.z, 1.25)*step(1.06, p)));\n"
"    d = min(d, mix(d, k3, step(p+n1.z, 1.99)*step(.5*pi, p)));\n"
"    d = min(d, abs(l-.41)-.03);\n"
"    d = min(d, mix(d, k4, step(p+n1.w, 1.04)*step(.04, p)));\n"
"    d = min(d, mix(d, k4, step(p+n1.w, -2.2)*step(-2.34, p)));\n"
"    \n"
"    dst = d-.005;\n"
"}\n"
"\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z)\n"
"{\n"
"    vec2 y = floor(x);\n"
"       float ret = 1.;\n"
"    vec2 pf=c.yy, p;\n"
"    float df=10.;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            d = length(x-p);\n"
"            \n"
"            if(d < df)\n"
"            {\n"
"                df = d;\n"
"                pf = p;\n"
"            }\n"
"        }\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            vec2 o = p - pf;\n"
"            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);\n"
"            ret = min(ret, d);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    z = pf;\n"
"}\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    vec3 c1;\n"
"    vec2 ind,\n"
"        xv,\n"
"        xi;\n"
"    float d,\n"
"        vs = 16.,\n"
"        n,\n"
"        size = .1,\n"
"        xix = mod(x.x, size)-.5*size,\n"
"        xixj = (x.x - xix),\n"
"        ri,\n"
"        rim1,\n"
"        rip1,\n"
"        lines = 8.,\n"
"        da,\n"
"        op,\n"
"        s;\n"
"    \n"
"    // Background blending\n"
"    s = smoothstep(0.,.5,.5-abs(x.y));\n"
"    col = mix(1.e-4*c.xxx, vec3(0.04,0.18,0.24), s);\n"
"    \n"
"    // Background circles\n"
"    dvoronoi(vs*x, d, ind);\n"
"    xv = ind/vs-x;\n"
"    lfnoise(vec2(3.,33.)*ind/vs-3.*iTime*c.xy,n);\n"
"    n = .5+.5*n;\n"
"    d = length(xv)-mix(.0,.35,n)/vs;\n"
"    col = mix(col, n*.5*vec3(1.00,0.40,0.39), sm(d));\n"
"    d = abs(d-.005) -.002;\n"
"    col = mix(col, (1.-n)*vec3(0.49,0.71,0.78), sm(d));\n"
"    \n"
"    for(float i = 1.; i < 9.; i += 1.)\n"
"    {\n"
"        rand((9.-i)*c.xx, op);\n"
"        op = .5+.5*round(16.*op)/16.;\n"
"        x += -.1+.2*op;\n"
"        \n"
"        xix = mod(x.x, size)-.5*size;\n"
"        xixj = (x.x - xix);\n"
"        \n"
"        // Edges\n"
"        lfnoise(2.e0*xixj*c.xx+14.*i, ri);\n"
"        lfnoise(2.e0*(xixj+size)*c.xx+14.*i, rip1);\n"
"        lfnoise(2.e0*(xixj-size)*c.xx+14.*i, rim1);\n"
"\n"
"        float h = .2;\n"
"        \n"
"        ri = h*round(lines*ri)/lines;\n"
"        rip1 = h*round(lines*rip1)/lines;\n"
"        rim1 = h*round(lines*rim1)/lines;\n"
"\n"
"        //if(ri < 0.)\n"
"        {\n"
"            dlinesegment(vec2(xix, x.y), vec2(-.5*size, mix(ri,rim1,.5)), vec2(-.25*size, ri), d);\n"
"            dlinesegment(vec2(xix, x.y), vec2(-.25*size, ri), vec2(.25*size, ri), da);\n"
"            d = min(d, da);\n"
"            dlinesegment(vec2(xix, x.y), vec2(.25*size, ri), vec2(.5*size, mix(ri,rip1,.5)), da);\n"
"            d = min(d, da);\n"
"            stroke(d, .002+.002*op, d);\n"
"            col = mix(col, op*(1.-n)*vec3(1.00,0.40,0.39), sm(d));\n"
"\n"
"            // Dots\n"
"            lfnoise(8.*xixj*c.xx-3.*iTime*c.xy+14.*i, n);\n"
"            n = .5+.5*n;\n"
"            d = length(vec2(xix, x.y-ri))-mix(.0,.35,n)/vs;\n"
"            c1 = mix(vec3(1.00,0.40,0.39), vec3(0.85,0.87,0.89), n);\n"
"            col = mix(col, op*(1.-n)*c1, sm(d));\n"
"            stroke(d - .009, (1.-n)*.005, d);\n"
"            c1 *= 2.4;\n"
"            col = mix(col, op*(1.-n)*c1, sm(d));\n"
"        }\n"
"        \n"
"        x -= -.1+.2*op;\n"
"    }\n"
"    \n"
"    //mix to blackish\n"
"    lfnoise(3.*x.xy-vec2(1.,.1)*iTime, n);\n"
"    stroke(n, .3, n);\n"
"    col = mix(col, 1.e-4*c.xxx, n);\n"
"    col = mix(col, .1*col, 1.-s);\n"
"    \n"
"    col = mix(col, mix(col, vec3(1.00,0.40,0.39), mix(.4,.8,.5+.5*x.y/.1)), sm(abs(x.y)-.1));\n"
"    col = mix(col, c.xxx, sm(abs(abs(x.y)-.11)-.001));\n"
"    \n"
"    col = mix(col, col*col, clamp(-x.y/.1,0.,1.));\n"
"    col *= col;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    \n"
"    float d;\n"
"\n"
"    vec4 old = vec4(-1.,texture(iChannel0, fragCoord/iResolution.xy).rgb), \n"
"    new = old; // Scene\n"
"    \n"
"    // Add overlay\n"
"    colorize(2.*(c.xz*uv-.45*vec2(-a,1.)-12.*c.xy), new.gba);\n"
"    new.gba = mix(old.gba, mix(old.gba, new.gba,.4), step(5.e-2,length(new.gba)));\n"
"    \n"
"    // Add Static text\n"
"    dstring(uv-.45*vec2(-.85*a,1.), 3., .02, d); // Team210\n"
"    \n"
"    stroke(d, .003, d);\n"
"\n"
"    new.gba = mix(new.gba, mix(new.gba, c.xxx, .6), sm(d));\n"
"    \n"
"    stroke(d-.002, .001, d);\n"
"    new.gba = mix(new.gba, vec3(1.00,0.40,0.39), sm(d));\n"
"\n"
"//     \n"
"    // Display time\n"
"    /*\n"
"    vec4 b = vec4(1., vec3(0.99,0.64,0.02)), bd = vec4(1., .5*vec3(0.99,0.64,0.02));\n"
"    box(uv-vec2(-.48,.45)-.03*sin(iTime)*c.xy, vec2(.2,.02), b.x);    \n"
"    stroke(b.x, .001, bd.x);\n"
"    add(b, bd, b);\n"
"    box(uv-vec2(-.08,.45)-.03*sin(iTime)*c.xy, vec2(.2,.02), bd.x);\n"
"    bd.gba = vec3(0.44,0.07,0.66);\n"
"    add(b, bd, b);\n"
"    stroke(bd.x, .001, bd.x);\n"
"    add(b, bd, b);\n"
"    dfloat(uv-vec2(-.63,.45)-.03*sin(iTime)*c.xy, iTime, .018, bd.x);\n"
"    stroke(bd.x, .004, bd.x);\n"
"    add(b, bd, b);\n"
"//     dfloat(uv-vec2(-.23,.45)-.03*sin(iTime)*c.xy, iExecutableSize, .018, bd.x);\n"
"    dstring(uv-vec2(-.225,.45)-.03*sin(iTime)*c.xy, 4., .018, bd.x);\n"
"    stroke(bd.x, .004, bd.x);\n"
"    bd.gba = vec3(0.99,0.64,0.02);\n"
"    add(b, bd, b);\n"
"    b.gba = mix(old.gba, b.gba, .8);\n"
"    \n"
"    blendadd(old, b, 5., 999., old);\n"
"    \n"
"    vec4 bda;\n"
"    \n"
"    if(iTime < 49.655)\n"
"    {\n"
"        dstring(uv+.6*c.xy, 1., .05, d); //Meme210 present\n"
"        stroke(d, .01, d);\n"
"        new = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(old,new,2.,8.,new);\n"
"        \n"
"        dstring(uv+.6*c.xy+.1*c.yx, 2., .03, d); // no partycoding this time!\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,4.,8.,new);\n"
"        \n"
"        dstring(uv+.6*c.xy+.15*c.yx, 3., .02, d); // well, shit.\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,6.,8.,new);\n"
"        \n"
"        dstring(uv+.6*c.xy, 5., .03, d); // Anyway, we present\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,9.,15.,new);\n"
"    \n"
"        dstring(uv+.6*c.xy-.15*c.yx, 6., .07, d); // GROSS GLOSS\n"
"        stroke(d, .015, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,11.,15.,new);\n"
"\n"
"        dstring(uv+.6*c.xy-.35*c.yx, 7., .03, d); // QM\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,30.,37.,new);\n"
"\n"
"        dstring(uv+.3*c.xy-.35*c.yx, 8., .03, d); // NR4\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,31.,38.,new);\n"
"\n"
"        dstring(uv+.0*c.xy-.35*c.yx, 9., .03, d); // MIC\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,32.,39.,new);\n"
"\n"
"        dstring(uv-.3*c.xy-.35*c.yx, 10., .03, d); // grenzdevil\n"
"        stroke(d, .005, d);\n"
"        old = vec4(d, mix(old.gba, c.xxx, .6));\n"
"        blendadd(new,old,33.,40.,new);\n"
"\n"
"    }\n"
"    else \n"
"    {\n"
"        new = old;\n"
"    }\n"
"    */\n"
"    \n"
"    fragColor = vec4(new.gba, 1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\n"
;
#endif
