//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

extern float progress;int scale_handle, dsmoothvoronoi_handle, rand_handle, hash31_handle, lfnoise_handle, mfnoise_handle, dbox_handle, dlinesegment3_handle, stroke_handle, zextrude_handle, add_handle, smoothmin_handle, dspline3_handle, dvoronoi_handle, normal_handle, dbox3_handle, rot3_handle, dtriangle_handle, dlinesegment_handle, dpolygon_handle, rot_handle, dcircle_handle, dschnappsgirls_handle, dspacepigs_handle, dkewlers_handle, dfarbrausch_handle, dhaujobb_handle, dmercury_handle, rshort_handle, rfloat_handle, drhomboid_handle, dcirclesegment_handle, dglyph_handle, dstring_handle, dfloat_handle, dint_handle, dtime_handle, window_handle, progressbar_handle, hash13_handle;
const int nsymbols = 40;
const char *scale_source = "#version 130\n\n"
"uniform float iTime;\n"
"void scale(out float s)\n"
"{\n"
"    if(iTime >=  0.0  && iTime <  4.705882 )\n"
"    {\n"
"        s = mod(iTime+.3- 0.0 , 0.588225 )- 0.2941125 ;\n"
"        s = smoothstep( -0.04901875 ,0.,s)*(1.-smoothstep(0., 0.14705625 ,s));\n"
"    }\n"
"    if(iTime >=  4.705882  && iTime <  18.552036 )\n"
"    {\n"
"        s = mod(iTime+.3- 4.705882 , 0.576975 )- 0.2884875 ;\n"
"        s = smoothstep( -0.04808125 ,0.,s)*(1.-smoothstep(0., 0.14424375 ,s));\n"
"    }\n"
"    if(iTime >=  18.552036  && iTime <  22.996481 )\n"
"    {\n"
"        s = mod(iTime+.3- 18.552036 , 0.55555 )- 0.277775 ;\n"
"        s = smoothstep( -0.046295833333333335 ,0.,s)*(1.-smoothstep(0., 0.1388875 ,s));\n"
"    }\n"
"    if(iTime >=  22.996481  && iTime <  25.139338 )\n"
"    {\n"
"        s = mod(iTime+.3- 22.996481 , 0.535675 )- 0.2678375 ;\n"
"        s = smoothstep( -0.04463958333333334 ,0.,s)*(1.-smoothstep(0., 0.13391875 ,s));\n"
"    }\n"
"    if(iTime >=  25.139338  && iTime <  27.208303 )\n"
"    {\n"
"        s = mod(iTime+.3- 25.139338 , 0.517275 )- 0.2586375 ;\n"
"        s = smoothstep( -0.043106250000000006 ,0.,s)*(1.-smoothstep(0., 0.12931875 ,s));\n"
"    }\n"
"    if(iTime >=  27.208303  && iTime <  65.208303 )\n"
"    {\n"
"        s = mod(iTime+.3- 27.208303 , 0.5 )- 0.25 ;\n"
"        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));\n"
"    }\n"
"    if(iTime >=  65.208303  && iTime <  71.109943 )\n"
"    {\n"
"        s = mod(iTime+.3- 65.208303 , 0.491825 )- 0.2459125 ;\n"
"        s = smoothstep( -0.04098541666666667 ,0.,s)*(1.-smoothstep(0., 0.12295625 ,s));\n"
"    }\n"
"    if(iTime >=  71.109943  && iTime <  82.722846 )\n"
"    {\n"
"        s = mod(iTime+.3- 71.109943 , 0.48385 )- 0.241925 ;\n"
"        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));\n"
"    }\n"
"    if(iTime >=  82.722846  && iTime <  86.722846 )\n"
"    {\n"
"        s = mod(iTime+.3- 82.722846 , 0.5 )- 0.25 ;\n"
"        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));\n"
"    }\n"
"    if(iTime >=  86.722846  && iTime <  94.998708 )\n"
"    {\n"
"        s = mod(iTime+.3- 86.722846 , 0.517275 )- 0.2586375 ;\n"
"        s = smoothstep( -0.043106250000000006 ,0.,s)*(1.-smoothstep(0., 0.12931875 ,s));\n"
"    }\n"
"    if(iTime >=  94.998708  && iTime <  103.134301 )\n"
"    {\n"
"        s = mod(iTime+.3- 94.998708 , 0.50845 )- 0.254225 ;\n"
"        s = smoothstep( -0.04237083333333333 ,0.,s)*(1.-smoothstep(0., 0.1271125 ,s));\n"
"    }\n"
"    if(iTime >=  103.134301  && iTime <  105.134301 )\n"
"    {\n"
"        s = mod(iTime+.3- 103.134301 , 0.5 )- 0.25 ;\n"
"        s = smoothstep( -0.041666666666666664 ,0.,s)*(1.-smoothstep(0., 0.125 ,s));\n"
"    }\n"
"    if(iTime >=  105.134301  && iTime <  107.069785 )\n"
"    {\n"
"        s = mod(iTime+.3- 105.134301 , 0.48385 )- 0.241925 ;\n"
"        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));\n"
"    }\n"
"    if(iTime >=  107.069785  && iTime <  129.569785 )\n"
"    {\n"
"        s = mod(iTime+.3- 107.069785 , 0.468775 )- 0.2343875 ;\n"
"        s = smoothstep( -0.03906458333333333 ,0.,s)*(1.-smoothstep(0., 0.11719375 ,s));\n"
"    }\n"
"    if(iTime >=  129.569785  && iTime <  131.505269 )\n"
"    {\n"
"        s = mod(iTime+.3- 129.569785 , 0.48385 )- 0.241925 ;\n"
"        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));\n"
"    }\n"
"    if(iTime >=  131.505269  && iTime <  141.341334 )\n"
"    {\n"
"        s = mod(iTime+.3- 131.505269 , 0.491825 )- 0.2459125 ;\n"
"        s = smoothstep( -0.04098541666666667 ,0.,s)*(1.-smoothstep(0., 0.12295625 ,s));\n"
"    }\n"
"    if(iTime >=  141.341334  && iTime <  143.276818 )\n"
"    {\n"
"        s = mod(iTime+.3- 141.341334 , 0.48385 )- 0.241925 ;\n"
"        s = smoothstep( -0.040320833333333334 ,0.,s)*(1.-smoothstep(0., 0.1209625 ,s));\n"
"    }\n"
"    if(iTime >=  143.276818  && iTime <  145.151818 )\n"
"    {\n"
"        s = mod(iTime+.3- 143.276818 , 0.468775 )- 0.2343875 ;\n"
"        s = smoothstep( -0.03906458333333333 ,0.,s)*(1.-smoothstep(0., 0.11719375 ,s));\n"
"    }\n"
"    if(iTime >=  145.151818  && iTime <  186.97 )\n"
"    {\n"
"        s = mod(iTime+.3- 145.151818 , 0.45455 )- 0.227275 ;\n"
"        s = smoothstep( -0.037879166666666665 ,0.,s)*(1.-smoothstep(0., 0.1136375 ,s));\n"
"    }\n"
"}\n"
"\0";
const char *dsmoothvoronoi_source = "#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform float iFader0;\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void lfnoise(in vec2 t, out float n);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void rand(in vec2 x, out float d);\n"
"void dsmoothvoronoi(in vec2 x, out float d, out vec2 z)\n"
"{\n"
"    float n;\n"
"//     lfnoise(x-iTime*c.xy, n);\n"
"    \n"
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
"            smoothmin(ret, d, .05, ret);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    z = pf;\n"
"}\n"
"\0";
const char *rand_source = "#version 130\n\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\0";
const char *hash31_source = "// Creative Commons Attribution-ShareAlike 4.0 International Public License\n"
"// Created by David Hoskins.\n"
"// See https://www.shadertoy.com/view/4djSRW\n"
"void hash31(in float p, out vec3 d)\n"
"{\n"
"   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));\n"
"   p3 += dot(p3, p3.yzx+33.33);\n"
"   d = fract((p3.xxy+p3.yzz)*p3.zyx); \n"
"}\n"
"\0";
const char *lfnoise_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void rand(in vec2 x, out float d);\n"
"void lfnoise(in vec2 t, out float n)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    t = smoothstep(c.yy, c.xx, t);\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    n = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\0";
const char *mfnoise_source = "#version 130\n\n"
"// const vec3 c = vec3(1.,0.,-1.);\n"
"void lfnoise(in vec2 x, out float d);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)\n"
"{\n"
"    n = 0.;\n"
"    float a = 1., nf = 0., buf;\n"
"    for(float f = d; f<b; f *= 2.)\n"
"    {\n"
"        lfnoise(f*x, buf);\n"
"        n += a*buf;\n"
"        a *= e;\n"
"        nf += 1.;\n"
"    }\n"
"    n *= (1.-e)/(1.-pow(e, nf));\n"
"}\n"
"\0";
const char *dbox_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"    vec2 da = abs(x)-b;\n"
"    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\0";
const char *dlinesegment3_source = "#version 130\n\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)\n"
"{\n"
"    vec3 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *stroke_source = "// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\0";
const char *zextrude_source = "// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\0";
const char *add_source = "void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = mix(sda, sdb, step(sdb.x, sda.x));\n"
"}\n"
"\0";
const char *smoothmin_source = "// iq's smooth minimum\n"
"void smoothmin(in float a, in float b, in float k, out float dst)\n"
"{\n"
"    float h = max( k-abs(a-b), 0.0 )/k;\n"
"    dst = min( a, b ) - h*h*h*k*(1.0/6.0);\n"
"}\n"
"\0";
const char *dspline3_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"\n"
"//distance to spline with parameter t\n"
"float dist3(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)\n"
"{\n"
"    t = clamp(t, 0., 1.);\n"
"    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);\n"
"}\n"
"\n"
"//minimum dist3ance to spline\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds)\n"
"{\n"
"    //coefficients for 0 = t^3 + a * t^2 + b * t + c\n"
"    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,\n"
"    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);\n"
"\n"
"	//discriminant and helpers\n"
"    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;\n"
"    \n"
"    //triple real root\n"
"    if(dis > 0.) \n"
"    {\n"
"        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);\n"
"        ds = dist3(p0,p1,p2,x,ui.x+ui.y-tau);\n"
"        return;\n"
"    }\n"
"    \n"
"    //three dist3inct real roots\n"
"    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;\n"
"    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;\n"
"    ds = min(\n"
"        dist3(p0,p1,p2,x, t.x),\n"
"        min(\n"
"            dist3(p0,p1,p2,x,t.y),\n"
"            dist3(p0,p1,p2,x,t.z)\n"
"        )\n"
"    );\n"
"}\n"
"\0";
const char *dvoronoi_source = "#version 130\n\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"void rand(in vec2 x, out float d);\n"
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
"\0";
const char *normal_source = "const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"void scene(in vec3 x, out vec2 s);\n"
"void normal(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\0";
const char *dbox3_source = "#version 130\n\n"
"void dbox3(in vec3 x, in vec3 b, out float d)\n"
"{\n"
"  vec3 da = abs(x) - b;\n"
"  d = length(max(da,0.0))\n"
"         + min(max(da.x,max(da.y,da.z)),0.0);\n"
"}\n"
"\0";
const char *rot3_source = "const vec3 c = vec3(1.,0.,-1.);\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\0";
const char *dtriangle_source = "// Adapted from iq, https://www.shadertoy.com/view/XsXSz4\n"
"void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst)\n"
"{\n"
"	vec2 e0 = p1 - p0;\n"
"	vec2 e1 = p2 - p1;\n"
"	vec2 e2 = p0 - p2;\n"
"\n"
"	vec2 v0 = p - p0;\n"
"	vec2 v1 = p - p1;\n"
"	vec2 v2 = p - p2;\n"
"\n"
"	vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );\n"
"	vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );\n"
"	vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );\n"
"    \n"
"    float s = sign( e0.x*e2.y - e0.y*e2.x );\n"
"    vec2 d = min( min( vec2( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) ),\n"
"                       vec2( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) )),\n"
"                       vec2( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) ));\n"
"\n"
"	dst = -sqrt(d.x)*sign(d.y);\n"
"}\n"
"\0";
const char *dlinesegment_source = "#version 130\n\n"
"\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\0";
const char *dpolygon_source = "#version 130\n\n"
"// compute distance to regular polygon\n"
"const float pi = acos(-1.);\n"
"void dpolygon(in vec2 x, in float N, out float d)\n"
"{\n"
"    d = 2.0*pi/N;\n"
"    float t = mod(acos(x.x/length(x)), d)-0.5*d;\n"
"    d = -0.5+length(x)*cos(t)/cos(0.5*d);\n"
"}\n"
"\0";
const char *rot_source = "#version 130\n\n"
"void rot(in float phi, out mat2 m)\n"
"{\n"
"    vec2 cs = vec2(cos(phi), sin(phi));\n"
"    m = mat2(cs.x, -cs.y, cs.y, cs.x);\n"
"}\n"
"\0";
const char *dcircle_source = "#version 130\n\n"
"\n"
"void dcircle(in vec2 x, out float d)\n"
"{\n"
"    d = length(x)-1.0;\n"
"}\n"
"\0";
const char *dschnappsgirls_source = "#version 130\n\n"
"\n"
"void dtriangle(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dcircle(in vec2 x, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"\n"
"void dschnappsgirls(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    // Dress\n"
"    dtriangle(x, vec2(-.1,-.3), vec2(.5,-.3), vec2(.2, .6), d0);\n"
"    dlinesegment(x, vec2(-.1,.325), vec2(.5,.325), da);\n"
"    stroke(da,.06,da);\n"
"    d0 = max(d0,-da);\n"
"    \n"
"    // Head\n"
"    dcircle(7.*(x-vec2(.2,.5)), da);\n"
"    d0 = max(d0, -da+.5);\n"
"    d0 = min(d0, da/7.);\n"
"    \n"
"    // Legs\n"
"    dlinesegment(x, vec2(.125,-.3), vec2(.125,-.6), da);\n"
"    stroke(da, .06, da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x, vec2(.275,-.3), vec2(.275,-.6), da);\n"
"    stroke(da, .06, da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Shoulders\n"
"    dlinesegment(x, vec2(0.05,.25), vec2(.35,.25), da);\n"
"    stroke(da, .085, da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Arms\n"
"    dlinesegment(x, vec2(.385,.25), vec2(.5, -.1), da);\n"
"    stroke(da, .055, da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x, vec2(.017,.25), vec2(-.1, -.1), da);\n"
"    stroke(da, .055, da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Glass\n"
"    dtriangle(x, vec2(-.6,.3), vec2(-.4,.1), vec2(-.2,.3), da);\n"
"    stroke(da, .0125, da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x, vec2(-.4,.15), vec2(-.4,-.1), da);\n"
"    stroke(da, .0125, da);\n"
"    d0 = min(d0, da);\n"
"    dtriangle(x, vec2(-.5,-.15), vec2(-.3,-.15), vec2(-.4,-.1), da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Liquid\n"
"    dtriangle(x, vec2(-.55,.25), vec2(-.4,.1), vec2(-.25,.25), da);\n"
"    d0 = min(d0, da);\n"
"    \n"
"    // Salad\n"
"    dlinesegment(x, vec2(-.4,.1), vec2(-.2,.5), da);\n"
"    stroke(da, .01, da);\n"
"    d0 = min(d0, da);\n"
"    dcircle(24.*(x-vec2(-.3,.3)), da);\n"
"    d0 = min(d0, da/24.);\n"
"    dcircle(24.*(x-vec2(-.25,.4)), da);\n"
"    d0 = min(d0, da/24.);\n"
"    \n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dspacepigs_source = "#version 130\n\n"
"\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"void dcircle(in vec2 x, out float d);\n"
"\n"
"void dear(in vec2 x, out float d)\n"
"{\n"
"    d = abs(2.*x.y)\n"
"        -.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))\n"
"        -.5*min(-abs(x.x),.01);\n"
"}\n"
"\n"
"void dspacepigs(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    // Head\n"
"    dcircle(2.5*x,d0);\n"
"    d0 /= 2.5;\n"
"    \n"
"    // Ears\n"
"    dear(vec2(2.,5.)*x-vec2(.8,1.3), da);\n"
"    d0 = min(d0,da/10.);\n"
"    dear(vec2(2.,5.)*x+vec2(.8,-1.3), da);\n"
"    d0 = min(d0,da/10.);\n"
"    \n"
"    // Nose\n"
"    dcircle(6.*x-vec2(0.,-.5),da);\n"
"    d0 = max(d0,-da/6.);\n"
"    dcircle(24.*x-vec2(-1.5,-2.),da);\n"
"    d0 = min(d0,da/24.);\n"
"    dcircle(24.*x-vec2(1.5,-2.),da);\n"
"    d0 = min(d0,da/24.);\n"
"    \n"
"    // Eyes\n"
"    dcircle(16.*x-vec2(-3.5,2.5),da);\n"
"    d0 = max(d0,-da/16.);\n"
"    dcircle(16.*x-vec2(3.5,2.5),da);\n"
"    d0 = max(d0,-da/16.);\n"
"    dcircle(24.*x-vec2(-5.,3.5),da);\n"
"    d0 = min(d0,da/24.);\n"
"    dcircle(24.*x-vec2(5.,3.5),da);\n"
"    d0 = min(d0,da/24.);\n"
"    \n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dkewlers_source = "#version 130\n\n"
"\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"\n"
"void dkewlers(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    x *= 1.2;\n"
"    \n"
"    dbox(x-vec2(0.,-.3),vec2(.6,.1),d0);\n"
"    dbox(x-vec2(-.5,-.0),vec2(.1,.25),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x-vec2(-.5+1./3.,.25),vec2(.1,.5),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x-vec2(-.5+2./3.,-.0),vec2(.1,.25),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x-vec2(.5,-.0),vec2(.1,.25),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dfarbrausch_source = "#version 130\n\n"
"\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"\n"
"void dfarbrausch(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    \n"
"    x += vec2(.1,0.);\n"
"    x *= 1.2;\n"
"    \n"
"    dlinesegment(x,vec2(-.65,.05),vec2(-.5,.05),d0);\n"
"    dlinesegment(x,vec2(-.5,.05),vec2(-.2,-.49),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(-.2,-.49),vec2(-.0,-.49),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(-.0,-.49),vec2(-.27,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(-.07,0.),vec2(-.27,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.2,-.49),vec2(-.07,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.4,-.49),vec2(.13,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.4,-.49),vec2(.2,-.49),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.33,0.),vec2(.13,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.33,0.),vec2(.51,-.33),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.6,-.15),vec2(.51,-.33),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.53,0.),vec2(.6,-.15),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.7,0.),vec2(.53,.0),da);\n"
"    d0 = min(d0, da);\n"
"    dlinesegment(x,vec2(.7,0.),vec2(.68,-.04),da);\n"
"    d0 = min(d0, da);\n"
"    dpolygon(5.*(x+vec2(.3,.65)),6.,da);\n"
"    d0 = min(d0, da/5.);\n"
"    dpolygon(5.*(x+vec2(-.5,.65)),6.,da);\n"
"    d0 = min(d0, da/5.);\n"
"    \n"
"    stroke(d0,.035, d0);\n"
"    d = max(d, -d0);\n"
"}\n"
"\0";
const char *dhaujobb_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"void rot(in float phi, out mat2 m);\n"
"void dcircle(in vec2 x, out float d);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"\n"
"void dhaujobb(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da, d0;\n"
"    mat2 m;\n"
"	rot(.3,m);\n"
"    x = 1.1*x*m;\n"
"    x.x *= 1.1;\n"
"        \n"
"    x += vec2(-.05,.2);\n"
"    \n"
"    // Left leg\n"
"    dbox(x+.35*c.xx,vec2(.1,.05),d0);\n"
"    dbox(x+vec2(.3,.25),vec2(.05,.15),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x+vec2(.2,.15),vec2(.1,.05),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(x+vec2(.15,.05),vec2(.05,.15),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Right leg\n"
"    dbox(x-vec2(.65,.35),vec2(.05,.15),da);\n"
"    d0 = min(d0,da);\n"
"\n"
"    // Torso\n"
"    rot(.2, m);\n"
"    dbox(m*(x-vec2(.25,.15)),vec2(.45,.05),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(m*(x-vec2(-.15,.35)),vec2(.45,.05),da);\n"
"    d0 = min(d0,da);\n"
"    rot(pi/8.,m);\n"
"    dbox(m*(x-vec2(.0,.25)),vec2(.1,.15),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Penis\n"
"    dbox(m*(x-vec2(.1,-.0)),vec2(.025,.1),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Left hand\n"
"    rot(.3,m);\n"
"    dbox(m*(x-vec2(.235,.535)),vec2(.035,.15),da);\n"
"    d0 = min(d0,da);\n"
"    dbox(m*(x-vec2(.225,.7)),vec2(.075,.025),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Right hand\n"
"    rot(-.2,m);\n"
"    dbox(m*(x+vec2(.585,-.2)),vec2(.0375,.1),da);\n"
"    d0 = min(d0,da);\n"
"    \n"
"    // Head\n"
"    dcircle(6.*(x-vec2(-.15,.58)),da);\n"
"    d0 = min(d0,da/6.);\n"
"    \n"
"    d0 -= .05*(abs(x.x)+abs(x.y)-.2);\n"
"    d = max(d,-d0);\n"
"}\n"
"\0";
const char *dmercury_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dpolygon(in vec2 x, in float N, out float d);\n"
"\n"
"void dmercury(in vec2 x, out float d)\n"
"{\n"
"    dpolygon(.5*x,6.0,d);\n"
"    float da;\n"
"\n"
"    x += .1*c.yx;\n"
"\n"
"    // Upper part\n"
"    dbox(x-.35*c.yx,vec2(.4,.35), da);\n"
"    d = max(d, -da);\n"
"    dbox(x-.7*c.yx, vec2(.2,.2), da);\n"
"    d = min(d,da);\n"
"    dbox(x-.25*c.yx,vec2(.2,.05),da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Lower part\n"
"    dbox(x+.2*c.yx,vec2(.1,.4),da);\n"
"    d = max(d, -da);\n"
"    dbox(x+.2*c.yx, vec2(.4,.1),da);\n"
"    d = max(d, -da);\n"
"}\n"
"\0";
const char *rshort_source = "#version 130\n\n"
"\n"
"uniform float iFontWidth;\n"
"uniform sampler2D iFont;\n"
"\n"
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
"\0";
const char *rfloat_source = "#version 130\n\n"
"\n"
"void rshort(in float off, out float val);\n"
"\n"
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
"\0";
const char *drhomboid_source = "#version 130\n\n"
"\n"
"void dbox(in vec2 x, in vec2 b, out float dst);\n"
"\n"
"void drhomboid(in vec2 x, in vec2 b, in float tilt, out float dst)\n"
"{\n"
"    x.x -= tilt/2./b.y*x.y;\n"
"    dbox(x,b,dst);\n"
"}\n"
"\0";
const char *dcirclesegment_source = "#version 130\n\n"
"const float pi = acos(-1.);\n"
"void dcirclesegment(in vec2 x, in float R, in float p0, in float p1, out float d)\n"
"{\n"
"    float p = atan(x.y, x.x);\n"
"    vec2 philo = vec2(max(p0, p1), min(p0, p1));\n"
"    if((p < philo.x && p > philo.y) || (p+2.0*pi < philo.x && p+2.0*pi > philo.y) || (p-2.0*pi < philo.x && p-2.0*pi > philo.y))\n"
"        d = abs(length(x)-R);\n"
"    else d = min(\n"
"        length(x-vec2(cos(p0), sin(p0))),\n"
"        length(x-vec2(cos(p1), sin(p1)))\n"
"        );\n"
"}\n"
"\0";
const char *dglyph_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void rfloat(in float off, out float val);\n"
"void dbox(in vec2 x, in vec2 b, out float dst);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void dcircle(in vec2 x, out float d);\n"
"void dcirclesegment(in vec2 x, in float r, in float p0, in float p1, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst)\n"
"{\n"
"    float dis;\n"
"    dbox(x, 2.*size*c.xx, dis);\n"
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
"    float d = 1., da = 1.;\n"
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
"        dlinesegment(x, size*vec2(x1,y1), size*vec2(x2, y2), da);\n"
"        d = min(d,da);\n"
"    }\n"
"    \n"
"    stroke(d,.2*size,d);\n"
"    \n"
"    // Smooth lines\n"
"    float nsmoothlines, db = 1.;\n"
"    da = 1.;\n"
"    rfloat(offset, nsmoothlines);\n"
"    nsmoothlines = floor(nsmoothlines);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<nsmoothlines; i+=1.)\n"
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
"        dlinesegment(x, size*vec2(x1,y1), size*vec2(x2, y2), db);\n"
"        da = min(da, db);\n"
"    }\n"
"    stroke(da,.2*size, da);\n"
"    smoothmin(d,da,.1*size,d);\n"
"    \n"
"    \n"
"//     if(nlines+nsmoothlines== 0.)\n"
"//         dst = dis;\n"
"//     else dst = d;\n"
"    dst = d;\n"
"}\n"
"\0";
const char *dstring_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void rfloat(in float off, out float val);\n"
"void dbox(in vec2 x, in vec2 b, out float dst);\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst);\n"
"\n"
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
"    dbox(x-size*(len-3.)*c.xy, vec2(size*len, 1.*size), bound);\n"
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
"\0";
const char *dfloat_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst);\n"
"\n"
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
"\0";
const char *dint_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst);\n"
"\n"
"void dint(in vec2 x, in float num, in float size, in float ndigits, out float dst)\n"
"{\n"
"    float d = 1., index = 0.;\n"
"    \n"
"    if(num == 0.)\n"
"    {\n"
"        index = ndigits;\n"
"        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48., .7*size, dst);\n"
"        return;\n"
"    } \n"
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
"    \n"
"    int hit = 0;\n"
"    \n"
"    // Determine the significand and output it\n"
"    for(float i = ndigits; i >= 0.; i -= 1.)\n"
"    {\n"
"        float po = pow(10.,i);\n"
"        float ca = floor(num/po);\n"
"        if(ca == 0.) \n"
"        {\n"
"            if(hit == 0)\n"
"            {\n"
"                index += 1.;\n"
"                continue;\n"
"            }\n"
"            \n"
"        }\n"
"        else hit = 1;\n"
"        num -= ca*po;\n"
"        \n"
"        float da;\n"
"        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"        d = min(d, da);\n"
"        index += 1.;\n"
"    }\n"
"    \n"
"    \n"
"    \n"
"    dst = d;\n"
"}\n"
"\0";
const char *dtime_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst);\n"
"\n"
"// Time in format: 00:00\n"
"void dtime(in vec2 x, in float num, in float size, out float dst)\n"
"{\n"
"    float d = 1., index = 0.;\n"
"    \n"
"    num = floor(num);\n"
"\n"
"    // 10 minutes\n"
"    float ca = floor(num/600.), da = 1.;\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"    d = min(d, da);\n"
"    index += 1.;\n"
"    num -= ca*600.;\n"
"\n"
"    // minutes\n"
"    ca = floor(num/60.);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"    d = min(d, da);\n"
"    index += 1.;\n"
"    num -= ca*60.;\n"
"\n"
"    // :\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 58., .7*size, da);\n"
"    d = min(d, da);\n"
"    index += 1.;\n"
"    \n"
"    // 10 seconds\n"
"    ca = floor(num/10.);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"    d = min(d, da);\n"
"    index += 1.;\n"
"    num -= ca*10.;\n"
"    \n"
"    // 1 seconds\n"
"    ca = floor(num/1.);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"    d = min(d, da);\n"
"    \n"
"    dst = d;\n"
"}\n"
"\0";
const char *window_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"uniform float iTime;\n"
"\n"
"void dhexagonpattern(in vec2 p, out float d, out vec2 ind);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void lfnoise(in vec2 t, out float num);\n"
"void box(in vec2 x, in vec2 b, out float dst);\n"
"void drhomboid(in vec2 x, in vec2 b, in float tilt, out float dst);\n"
"\n"
"// Fixme: Add\n"
"\n"
"void window(in vec2 x, in vec2 size, in vec3 bg, in float title_index, out vec4 col)\n"
"{\n"
"//     size.x *= .5;\n"
"//     col = vec4(1., bg);\n"
"//     \n"
"//     const float cellsize = .015, bordersize = .005;\n"
"//     vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),\n"
"//         bordercolor = vec3(1.00,0.71,0.02);\n"
"//     vec4 c2 = vec4(1., titlecolor);\n"
"//     \n"
"//     float dhx, dhy;\n"
"//     vec2 ind;\n"
"//     dhexagonpattern(72.*x,  dhx, ind);\n"
"//     stroke(dhx, .1, dhx);\n"
"//     lfnoise(ind-iTime, dhy);\n"
"//     \n"
"//     // Window background\n"
"//     box(x+.5*size*c.yx,size*vec2(1.,.5),c2.x);\n"
"//     c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/size.y), .5+.5*dhy*step(0.,dhx));\n"
"//     add(col, c2, col);\n"
"//     \n"
"//     // Title bar\n"
"//     c2.gba = titlecolor;\n"
"//     drhomboid(x+.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);\n"
"//    	add(col, c2, col);\n"
"//     drhomboid(x, vec2(.65*size.x,cellsize), cellsize, c2.x);\n"
"//    	add(col, c2, col);\n"
"//     drhomboid(x-.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);\n"
"//    	add(col, c2, col);\n"
"//     \n"
"//     // Border of title bar\n"
"//     c2 = vec4(1., bordercolor);\n"
"//     stroke(col.x,bordersize,c2.x);\n"
"//     add(col,c2,col);\n"
"//     \n"
"//     // Window Border\n"
"//     dlinesegment(x, -.9*size.x*c.xy, -size.x*c.xy, c2.x);\n"
"//     float d;\n"
"//     dlinesegment(x, -size.x*c.xy, -size, d);\n"
"//     c2.x = min(c2.x, d);\n"
"//     dlinesegment(x, -size, size*c.xz, d);\n"
"//     c2.x = min(c2.x, d);\n"
"//     dlinesegment(x, size*c.xz, size*c.xy, d);\n"
"//     c2.x = min(c2.x, d);\n"
"//     dlinesegment(x, .9*size.x*c.xy, size.x*c.xy, d);\n"
"//     c2.x = min(c2.x, d);\n"
"//     stroke(c2.x,.25*bordersize,c2.x);\n"
"//     add(col, c2, col);\n"
"}\n"
"\0";
const char *progressbar_source = "#version 130\n\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void progressbar(in vec2 x, in float width, in float progress, out vec4 col)\n"
"{\n"
"//     const float cellsize = .015, bordersize = .005;\n"
"//     vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),\n"
"//         bordercolor = vec3(1.00,0.71,0.02), bg = c.yyy;\n"
"//     vec4 c2 = vec4(1., titlecolor);\n"
"//     \n"
"//     // Window background\n"
"//     box(x+.5*width*c.yx,width*c.xy,c2.x);\n"
"//     c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/cellsize), .5);\n"
"//     add(col, c2, col);\n"
"//     \n"
"//     // Bar background\n"
"//     c2.gba = titlecolor;\n"
"//     drhomboid(x, vec2(.5*width,cellsize), cellsize, c2.x);\n"
"//    	add(col, c2, col);\n"
"//     \n"
"//     // Border\n"
"//     c2.gba = bordercolor;\n"
"//     stroke(c2.x,.5*bordersize,c2.x);\n"
"//     add(col, c2, col);\n"
"//     \n"
"//     // Progress\n"
"//     float wc = width/cellsize;\n"
"//     x.x -= .5*x.y;\n"
"//     vec2 y = vec2(mod(x.x, 1.2*cellsize)-.6*cellsize, x.y),\n"
"//         index = (x-y)/.6/cellsize;\n"
"//     if(abs(index.x) < .8*wc && -index.x > .8*wc*(1.-2.*progress))\n"
"//     {\n"
"//         box(y, vec2(.5*cellsize, .8*cellsize), c2.x);\n"
"//         add(col, c2, col);\n"
"//     }\n"
"}\n"
"\0";
const char *hash13_source = "// Creative Commons Attribution-ShareAlike 4.0 International Public License\n"
"// Created by David Hoskins.\n"
"// See https://www.shadertoy.com/view/4djSRW\n"
"void hash13(in vec3 p3, out float d)\n"
"{\n"
"	p3  = fract(p3 * .1031);\n"
"    p3 += dot(p3, p3.yzx + 33.33);\n"
"    d = fract((p3.x + p3.y) * p3.z);\n"
"}\n"
"\0";
const char *voronoidesign_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"void scale(out float s);\n"
"void dsmoothvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void rand(in vec2 x, out float n);\n"
"void hash31(in float p, out vec3 d);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"\n"
"vec2 vind,vind2;\n"
"float v, fn, r1, fb;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.y += mix(.2,-.2,step(150., iTime))*iTime;\n"
"    \n"
"    dvoronoi(1.5*x.xy, v, vind);\n"
"    \n"
"    vec3 y = vec3(vind/1.5-x.xy,x.z);\n"
"    \n"
"    float n, n2;\n"
"    \n"
"    lfnoise(c.xx-.3*iTime+vind*3., n);\n"
"    lfnoise(5.*x.z*c.xx-iTime-vind*4., n2);\n"
"    n2 *= .2;\n"
"    \n"
"    mat2 RR = mat2(cos(n2), sin(n2), -sin(n2), cos(n2));\n"
"    vec2 a = x.xy;\n"
"    x.xy = RR * x.xy;\n"
"    rand(vind, r1);\n"
"    float r2;\n"
"    rand(vind -1336., r2);\n"
"    \n"
"    float phi = atan(y.y, y.x),\n"
"        dp = pi/24.,\n"
"        phii = mod(phi, dp)-.5*dp,\n"
"        pa = phi - phii, \n"
"        R1 = .05,\n"
"        R2 = mix(.4,.25,1.-r2);\n"
"    \n"
"    R2 = mix(R1, R2, .5+.5*n2);\n"
"    \n"
"    float r0;\n"
"    rand(pa*c.xx, r0);\n"
"    r0 = mix(r0,.5+.5*n,.5);\n"
"    \n"
"    dspline3(y, vec3(1.4*R1*cos(pa), 1.4*R1*sin(pa), -.5), vec3(R1*cos(pa), R1*sin(pa), .1*r1), vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), sdf.x);\n"
"    float da;\n"
"    dspline3(y, vec3(mix(R1,R2,.5)*cos(pa), mix(R1,R2,.5)*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1*r1), vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0), da);\n"
"    sdf.x = min(sdf.x, da);\n"
"    stroke(sdf.x, .25*mix(.02,.05, .5+.5*n2), sdf.x);\n"
"    sdf.y = 2.;\n"
"    \n"
"    add(sdf, vec2(length(y-vec3(R2*cos(pa), R2*sin(pa), .1-.4*r0))-.01, 3.), sdf);\n"
"    \n"
"    float fa;\n"
"    lfnoise(4.*a,  fa);\n"
"    dvoronoi(a,fn, vind2); \n"
"    fa = x.z+.4+.1*mix((v+fn),fa,.5);\n"
"    add(sdf, vec2(fa,4.), sdf);\n"
"    smoothmin(sdf.x, fa, .1, sdf.x);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"float nan;\n"
"void vs(in vec3 x, out vec2 sdf)\n"
"{\n"
"    vec2 vi;\n"
"    dsmoothvoronoi(3.*(x.xy+.2*iTime*c.yx), sdf.x, vi);\n"
"    sdf.x = x.z-.1-.2*sdf.x;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, \n"
"        s;\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    uv *= 2.;\n"
"    \n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx), \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 400,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = -(o.z-.3)/dir.z;\n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        vs(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        d += s.x;\n"
"    }\n"
"    float v1, rar, dx = 1.e-3;\n"
"    vec2 vi1, na;\n"
"    vec3 cv, l;\n"
"    x = o + d * dir;\n"
"    dsmoothvoronoi(3.*(x.xy+.2*iTime*c.yx), v1, vi1);\n"
"\n"
"    rand(vi1, rar);\n"
"    cv = mix(c.yyy,vec3(.23,.23,.23), rar);\n"
"    v1 = abs(v1)-.01;\n"
"    cv = mix(cv, c.yyy, sm(v1));\n"
"    v1 = abs(v1-.01)-.005;\n"
"    cv = mix(cv, c.xxx, sm(v1));\n"
"    vs(x,s);\n"
"    vs(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    vs(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    vs(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"    l = normalize(x+c.yyx);\n"
"    cv = .2*cv\n"
"        +.2*cv*abs(dot(l,n))\n"
"        +.4*cv*pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"    cv = mix(cv, 1.5*vec3(0.76,0.20,0.13), smoothstep(0.858, 1.02, dot(n, c.yyx)));\n"
"    dir = refract(dir, n, .98);\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"\n"
"        d += min(s.x,3.e-2);\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 1.e-4);\n"
"       \n"
"        if(s.y == 3.)\n"
"        {\n"
"            l = normalize(x+c.yyx);\n"
"            float r;\n"
"            \n"
"            col = mix(c.xxx, vec3(0.76,0.20,0.13), .8);\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .8*col * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"        }\n"
"        else if(s.y == 4.)\n"
"        {\n"
"            l = normalize(x+c.yyx);\n"
"            float r;\n"
"            rand(vind+vind2,r);\n"
"            \n"
"            col = mix(.023*c.xxx, vec3(0.76,0.40,0.23), r);\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            stroke(v, .01, v);\n"
"            stroke(fn, .01, fn);\n"
"            col = mix(col, c.yyy, sm(v));\n"
"            col = mix(col, col*col, sm(fn));\n"
"        }\n"
"		if(s.y == 2.)\n"
"        {\n"
"            l = normalize(x+c.yyx);\n"
"            float r;\n"
"            rand(vind, r);\n"
"            \n"
"            col = mix(mix(vec3(0.76,0.20,0.23),vec3(.18,.32,.13), r), vec3(0.23,0.23,0.23),  clamp((x.z)/r1/.1,0.,1.));\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            col = mix(col, 5.*col, .25*n.x*n.x);\n"
"        }\n"
"        \n"
"    }\n"
"    \n"
"    col *=3.6;\n"
"    col *= col;\n"
"\n"
"    if(s.y != 3.)\n"
"    {\n"
"        col = mix(length(col)/sqrt(3.)*c.xxx, col,.3);\n"
"    }\n"
"    \n"
"    col = mix(col,cv,.8);\n"
"    col = mix(col, c.yyy, smoothstep(1.,5.,d));\n"
"    \n"
"    col *= mix(1.,15.,mix(.28,.88, 0.*iScale));\n"
"    col *= col;\n"
"    \n"
"    col = mix(col, .01*col, smoothstep(-.6,-1.,uv.y));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}	\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *groundboxes_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
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
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"float iScale, nbeats;\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dspline3(in vec3 x, in vec3 p0, in vec3 p1, in vec3 p2, out float ds);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void scale(out float s);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"\n"
"vec2 ind=c.yy;\n"
"\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{    \n"
"    \n"
"    float d,\n"
"        s = .1;\n"
"\n"
"    sdf = c.xy;\n"
"    \n"
"    for(float size = .0; size <= .5; size += .025)\n"
"    {\n"
"        dbox3(x, size*c.xxx, d);\n"
"        stroke(d, .001, d);\n"
"        vec2 sda = vec2(d,3.+size);\n"
"\n"
"        float n;\n"
"        vec3 y = mod(x,.125*size)-.5*.125*size,\n"
"            yi = (x-y)/size;\n"
"		ind = yi.xy+yi.yz+yi.xz;\n"
"        lfnoise(3.6*ind+15.*size-1.1*iTime, n); // TODO: use shifts to 13.6 from 3.6\n"
"        \n"
"        if(n>-.3)\n"
"        {\n"
"            \n"
"            dbox3(y, .25*size*c.xxx, d);\n"
"            sda.x = max(sda.x, -d);\n"
"        }\n"
"        \n"
"        add(sdf, sda, sdf);\n"
"    }   \n"
"}  \n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    col = .5*c.xxx;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    mat3 RR;\n"
"    rot3(.2*iTime*vec3(1.1,1.4,1.6), RR);\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = RR * c.yyx,\n"
"        r = normalize(c.xyy), \n"
"        u = normalize(c.yxy),\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x,\n"
"        size = .301*c.xxx;\n"
"        \n"
"    int N = 150,\n"
"        i = 0;\n"
"    t = uv.x * r + uv.y * u;\n"
"    t = RR * t;\n"
"    dir = normalize(t-o);\n"
"    float d = 0., ra, ss, inside = 0., dd;\n"
"    \n"
"    for(ss = .5; ss >= .0; ss -= .025)\n"
"    {\n"
"        size = ss*c.xxx+1.e-4;\n"
"\n"
"        vec3 tlo = min((size-o)/dir,(-size-o)/dir); // Select 3 visible planes\n"
"        vec2 abxlo = abs(o.yz + tlo.x*dir.yz),\n"
"            abylo = abs(o.xz + tlo.y*dir.xz),\n"
"            abzlo = abs(o.xy + tlo.z*dir.xy);\n"
"\n"
"        vec4 dn = 100.*c.xyyy;\n"
"        \n"
"        dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));\n"
"        dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));\n"
"        dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));\n"
"\n"
"        inside += .05;\n"
"        \n"
"        if(ss == 3.)dd = dn.r;\n"
"        d = dn.r;\n"
"        \n"
"        x = o + d * dir;\n"
"        scene(x, s);\n"
"        if(s.x < 1.e-4)break;\n"
"        \n"
"                normal(x,n, 5.e-4);\n"
"        \n"
"        lfnoise(x.xy*vec2(3.,8.)-iTime, ra);\n"
"        \n"
"        vec3 f;\n"
"        float dd = 5.e-1;\n"
"        \n"
"        r = RR*c.xyy;\n"
"        f = RR*c.yzy;\n"
"        u = RR*c.yyx;\n"
"        \n"
"        vec3 dp = abs(vec3(dot(n,r), dot(n,f), dot(n,u)));\n"
"        \n"
"        if(dp.y < dd && dp.z < dd) n = r;\n"
"        else if(dp.x < dd && dp.z < dd) n = f;\n"
"        else if(dp.x < dd && dp.y < dd) n = u;\n"
"        s.y = -1.;\n"
"\n"
"        vec3 l = normalize(x+.03*normalize(x-o));\n"
"        \n"
"        vec3 c1 = col;\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            colorize(x.xy, c1);\n"
"            c1 = .1*c1\n"
"                + 1.*c1 * abs(dot(l,n))\n"
"                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"            \n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            \n"
"            c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));\n"
"            vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5*sin(2.*iScale*ra*x));\n"
"            c1 = mix(c1, c1, .5+.5*ra);\n"
"            c1 = .3*c.xxx;\n"
"            c1 = .1*c1\n"
"                + .4*c1 * abs(dot(l,n))\n"
"                + .3*c1 * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"        }\n"
"        else if(s.y >= 3.)\n"
"        {\n"
"            lfnoise(.1*ind+3.*s.y*c.xx-iTime, c1.x);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, c1.y);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, c1.z);\n"
"            float na;\n"
"            ind = x.xy + x.yz + x.zx;\n"
"            ind = mod(ind, .01)-.005;\n"
"            ind = x.xy + x.yz + x.zx-ind;\n"
"            rand(ind, na);\n"
"            na = 0.;\n"
"            c1 = .8+.2*c1;\n"
"            c1 *= na;\n"
"            c1 = mix(.1,mix(.2,.4,iScale),step(na,.05))*c1\n"
"                + mix(.1,.2,step(na,.95))*c1 * abs(dot(l,n))\n"
"                + .5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"        }\n"
"        else if(s.y == -1.)\n"
"        {\n"
"            lfnoise(.1*ind+3.*s.y*c.xx-iTime, c1.x);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime, c1.y);\n"
"            lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime, c1.z);\n"
"            c1 = .8+.2*c1;\n"
"            c1 = .1*c1\n"
"                + 1.*c1 * abs(dot(l,n))\n"
"                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"                \n"
"            vec3 dc;\n"
"\n"
"            vec3 zz = mod(x, .025)-.5*.025, zi = x-zz;\n"
"            rand(zi.xy+zi.yz+zi.zx, dc.x);\n"
"            rand(zi.xy+zi.yz+zi.zx+1337., dc.y);\n"
"            rand(zi.xy+zi.yz+zi.zx+2337., dc.z);\n"
"            float da;\n"
"            dbox3(zz, .01*c.xxx, da);\n"
"            stroke(da, .001, da);\n"
"            c1 = mix(c1, 1.2*c1*c1+.2*dc, sm(da));\n"
"            stroke(da-.002, .001, da);\n"
"            c1 = mix(c1, 1.6*c1*c1, sm(da));\n"
"        }\n"
"        else if(s.y == -2.)\n"
"        {\n"
"            float s = .05;\n"
"            c1 = .5*c.xxx;\n"
"            vec2 dd = mod(x.xz, s)-.5*s;\n"
"            stroke(dd.x, .005, dd.x);\n"
"            stroke(dd.y, .005, dd.y);\n"
"            c1 = mix(c1, c.xxx, sm(min(dd.x, dd.y)));\n"
"        }\n"
"        col = mix(col, c1, mix(.2,.6,iScale));\n"
"        col = .9*col;\n"
"    }\n"
"    \n"
"    if(s.x > 1.e-4 && uv.y<-2.e-4)\n"
"    {\n"
"        d = -(o.y+.375)/dir.y;\n"
"        x = o + d * dir;\n"
"        scene(x, s);\n"
"        i=N;\n"
"    }\n"
"    else x = o + d * dir;\n"
"    \n"
"    if(s.y == -1.)\n"
"    {\n"
"        vec3 c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));\n"
"    	\n"
"        col = mix(col, c1, .3+.5*inside);\n"
"    }\n"
"    \n"
"    col = mix(col, c.yyy, smoothstep(2., 3., d));\n"
"    \n"
"    col *= 1.2;\n"
"    col *= col;\n"
"    \n"
"    col *= col*col;\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *graffiti_source = "/* Lightcyber by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"float iScale, nbeats;\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n);\n"
"void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void rot3(in vec3 phi, out mat3 R); \n"
"void scale(out float s);\n"
"\n"
"void graf(in vec2 x, out float d)\n"
"{\n"
"    x.y *= .7;\n"
"    float size = .4,\n"
"        n,\n"
"        da;\n"
"    vec2 y = vec2(mod(x.x, size)-.5*size, x.y),\n"
"        yi = (x-y)/size,\n"
"        x1,\n"
"        x2,\n"
"        x3;\n"
"    \n"
"    dbox(y,vec2(.75,.75)*size, d);\n"
"    \n"
"    // lines\n"
"    rand(yi, n);\n"
"    x1 = vec2(-.5+.02+n*.96, .75)*size,\n"
"    x2 = vec2(.5-.02-n*.96, -.75)*size;\n"
"    x1.x = floor(5.*x1.x)/5.;\n"
"    x2.x = floor(5.*x2.x)/5.;\n"
"    x1.x = max(x1.x,-.4*size);\n"
"    x1.x = min(x1.x,.4*size);\n"
"    x2.x = max(x2.x,-.4*size);\n"
"    x2.x = min(x2.x,.4*size);\n"
"    dlinesegment(y, x1, x2, da);\n"
"    stroke(da, .02, da);\n"
"    d = max(d, -da);\n"
"    \n"
"    // upper triangles\n"
"    rand(yi+1337., n);\n"
"	x1 = vec2(-.55+n,.75)*size*1.05,\n"
"    x2 = vec2(.5-n,.75)*size*1.05,\n"
"    x3 = .75*(.8-n)*size*1.05*c.yx-.01*c.yx;\n"
"    x1 = round(15.*x1)/15.;\n"
"    x2 = round(15.*x2)/15.;\n"
"    x3 = round(15.*x3)/15.;\n"
"    dtriangle(y, x1, x2, x3, da);\n"
"    d = max(d, -da);\n"
"    \n"
"    // lower triangles\n"
"    rand(yi+2337., n);\n"
"    x1 = vec2(-.5+n,-.75)*size*1.05,\n"
"    x2 = vec2(.55-n,-.75)*size*1.05,\n"
"    x3 = -.75*(.8-n)*size*1.05*c.yx+.01*c.yx;\n"
"    x1 = round(15.*x1)/15.;\n"
"    x2 = round(15.*x2)/15.;\n"
"    x3 = round(15.*x3)/15.;\n"
"    dtriangle(y, x1, x2, x3, da);\n"
"    d = max(d, -da);\n"
"}\n"
"\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    \n"
"    x.x += .3*iTime;\n"
"    x *= 2.;\n"
"    \n"
"    vec3 n;\n"
"    lfnoise(x.x*c.xx-iTime, n.x);\n"
"    lfnoise(2.*x.x*c.xx-iTime-1337., n.y);\n"
"    lfnoise(x.x*c.xx+2.*iTime-2337., n.z);\n"
"\n"
"    x.yz += .1*vec2(cos(x.x), sin(x.x))*n.xy;\n"
"    \n"
"    mat3 RR;\n"
"    rot3(1.3*mix(.2,1.5, .5+.5*n.x)*n.z * c.xyy, RR);\n"
"    x = RR * x;\n"
"    \n"
"    x.z = abs(x.z);\n"
"    \n"
"    float d, da, db;\n"
"    \n"
"    graf(x.xy, d);\n"
"    stroke(d+mix(.01,.04, iScale), mix(.01,.04, iScale), da);\n"
"    //stroke(d,.01,da);\n"
"    \n"
"    float v;\n"
"    vec2 ind;\n"
"    dvoronoi(12.*x.xy, v, ind); // 12.\n"
"    \n"
"    zextrude(x.z, -d, .1-.1*v, d);\n"
"    \n"
"	sdf = vec2(d,1.);\n"
"    float modsize = .025,\n"
"		y = mod(d-.3-.02*iTime,modsize)-.5*modsize,\n"
"        yi = (d-y)/modsize;\n"
"    \n"
"    float na;\n"
"    lfnoise(2.*yi*c.xx-.3*iTime, na);\n"
"\n"
"    zextrude(x.z-.05*na, -y, mix(0.,.05+.05*na,iScale), d);\n"
"    stroke(d,.035,d);\n"
"        \n"
"    \n"
"    \n"
"    zextrude(x.z, -da, .25, da);\n"
"    \n"
"	add(sdf, vec2(da, 1.), sdf);\n"
"	\n"
"	lfnoise(5.*x.xy, da);\n"
"	    mfnoise(x.xy, 32., 422., .45, db);\n"
"        da = .5*(db+da);\n"
"		sdf.x -= .001*da; // .001\n"
"        stroke(da, .1, da);\n"
"        sdf.x -=.005*da; // .005\n"
"    add(sdf, vec2(d, 1.), sdf);\n"
"    add(sdf, vec2(x.z+.25,1.), sdf);\n"
"    \n"
"    // xa(t)\n"
"    float xa = mix(x.x+3.*a,x.x-3.*a,clamp(iTime/3.,0.,1.));\n"
"    xa = mix(xa,-xa, clamp((iTime-3.)/6., 0.,1.));\n"
"    sdf.x += mix(0., 2., smoothstep(-.5, .5, xa));\n"
"    \n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void colorize(in vec2 x, out vec3 col)\n"
"{\n"
"    x.x += .3*iTime;\n"
"    x *= 2.;\n"
"\n"
"    float n;\n"
"    lfnoise(x.x*c.xx-iTime, n);\n"
"    x.y += .3*cos(x.x)*n;\n"
" \n"
"    float d;\n"
"    graf(x, d);\n"
"    col = mix(col, mix(mix(vec3(0.85,0.87,0.89), c.xxx, step(50., iTime)), mix(vec3(0.04,0.18,0.24),vec3(0.00,0.20,0.36),step(50.,iTime)), clamp(abs(x.y/2.),0.,1.)), sm(d-.2));\n"
"    col = mix(col, mix(vec3(1.00,0.40,0.39), vec3(0.00,0.67,0.91), step(50., iTime)), sm(d));\n"
"    float da = d;\n"
"    stroke(d+mix(.01,.03, iScale), mix(.01,.04,iScale), d);\n"
"    //stroke(d,.01, d);\n"
"    col = mix(col, 1.4*col, sm(d));\n"
"    stroke(d, .001, d);\n"
"    col = mix(col, 1.3*col, sm(d));\n"
"    \n"
"    if(da < .02 && da > -.02)\n"
"    {\n"
"        lfnoise(5.*x, da);\n"
"	    mfnoise(x, 32., 422., .45, d);\n"
"        d = .5*(d+da);\n"
"		col = mix(col, vec3(0.27,0.27,0.27), sm(d));\n"
"        stroke(d, .1, d);\n"
"        col = mix(col, 1.5*col, sm(d));\n"
"    }\n"
"    \n"
"    col *= mix(1., 1.6, iScale);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"        \n"
"    float sc2 = 0.,\n"
"        sc3 = 0.;\n"
"        \n"
"    mat3 R;\n"
"    \n"
"    vec3 col = c.yyy, \n"
"        o = mix(1.,.5,smoothstep(0.,5.,clamp(iTime-71.,0.,5.)))*c.yzx,\n"
"        r = c.xyy,\n"
"        u = normalize(c.yxx),\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 250,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    vec3 c1;\n"
"    float d = -(o.z-.35)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(x.z<-.15)\n"
"        {\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,5.e-2);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 1.e-2);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            vec3 l = normalize(x+.5*c.yzx);\n"
"            colorize(x.xy, c1);\n"
"            c1 = .1*c1\n"
"                + 1.*c1 * abs(dot(l,n))\n"
"                + 1.5 * c1 * abs(pow(dot(reflect(x-l,n),dir),2.));\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            vec3 l = normalize(x+c.xzx);\n"
"            float r;\n"
"            lfnoise(x.xy, r);\n"
"            c1 = mix(vec3(0.99,0.43,0.15),vec3(0.44,0.07,0.66),sin(2.*iScale*r*x));\n"
"            c1 = .1*c1\n"
"                + .8*c1 * abs(dot(l,n))\n"
"                + 6.5*c1 * abs(pow(dot(reflect(x-l,n),dir),3.));\n"
"        }\n"
"        col = c1;\n"
"    }\n"
"    \n"
"    \n"
"    //col += col;\n"
"    \n"
"    col *= col*col;\n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    \n"
"    // Background geraffel\n"
"    if(length(col) < .001)\n"
"    {\n"
"        float v, ra, v2;\n"
"        vec2 ind, ind2;\n"
"        lfnoise(iTime*c.xx, ind2.x);\n"
"        lfnoise(iTime*c.xx-1337., ind2.y);\n"
"        dvoronoi(12.*(uv-.03*ind2), v, ind);\n"
"//         dvoronoi(43.*uv, v2, ind2);\n"
"//         ind += .1*ind2;\n"
"        rand(ind, ra);\n"
"        stroke(-v, .05, v);\n"
"        v = -v;\n"
"        col = mix(col, .3*ra*mix( .5*vec3(1.00,0.40,0.39), .05*c.xxx, clamp(tanh(1.5*length(uv)),0.,1.)), sm(v));\n"
"        col *= mix(13.,1.,smoothstep(0.,.5,clamp((iTime-6.),0.,1.)));\n"
"    }\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *greet_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
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
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"float iScale;\n"
"\n"
"// void dpolygon(in vec2 x, in float N, out float d);\n"
"// void rot(in float phi, out mat2 m);\n"
"// void dcircle(in vec2 x, out float d);\n"
"// void dbox(in vec2 x, in vec2 b, out float d);\n"
"// void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"// void dtriangle(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, out float dst);\n"
"void scale(out float s);\n"
"void rand(in vec2 x, out float n);\n"
"void hash31(in float p, out vec3 d);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void dschnappsgirls(in vec2 x, out float d);\n"
"void dspacepigs(in vec2 x, out float d);\n"
"void dkewlers(in vec2 x, out float d);\n"
"void dfarbrausch(in vec2 x, out float d);\n"
"void dhaujobb(in vec2 x, out float d);\n"
"void dmercury(in vec2 x, out float d);\n"
"\n"
"vec2 ind, indc;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    float d;\n"
"\n"
"    d = mix(mix(mix(mix(mix(mix(mix(0.,0.14173228346456693,smoothstep(0.,0.341334,iTime)),\n"
"            0.14173228346456693+.25,smoothstep(0.341334,2.276818,iTime)),\n"
"            0.14173228346456693+.5,smoothstep(2.276818,4.151818,iTime)),\n"
"            0.14173228346456693+.75,smoothstep(4.151818,4.151818+1*1.8182,iTime)),\n"
"            0.14173228346456693+1.,smoothstep(4.151818+1*1.8182,4.151818+2*1.8182,iTime)),\n"
"            0.14173228346456693+1.25,smoothstep(4.151818+2*1.8182,4.151818+3*1.8182,iTime)),\n"
"            0.,smoothstep(4.151818+3*1.8182,4.151818+5*1.8182,iTime));\n"
"//     x.z -= mix(0.,1.5,iFader0);\n"
"    x.z -= d;\n"
"    \n"
"    // Corridor\n"
"    dbox3(x, vec3(.1,.1,1.e3), d);\n"
"    sdf = vec2(-d, 2.);\n"
"    \n"
"    // Wall tiles\n"
"    float distortion;\n"
"	lfnoise(5.2e2*x.yz, distortion);\n"
"    float tsize = .005,\n"
"    	dy = mod(x.y, tsize)-.5*tsize,\n"
"        yi = (x.y-dy)/tsize,\n"
"        zpar = x.z+mix(0., .5*tsize, mod(yi,2.)),\n"
"        dz = mod(zpar, tsize)-.5*tsize,\n"
"        zi = (zpar-dz)/tsize;\n"
"    dbox3(vec3(abs(x.x)-.1, dy, dz), vec3(.0005+.00001*distortion, .39*tsize*c.xx), d);\n"
"    add(sdf, vec2(d, 3.), sdf);\n"
"    smoothmin(sdf.x, d, .001, sdf.x);\n"
"    ind = vec2(yi, zi);\n"
"    \n"
"    // Ceiling \n"
"    tsize = .025;\n"
"    dz = mod(x.z, tsize)-.5*tsize;\n"
"    float dx = mod(x.x, tsize)-.5*tsize;\n"
"    zi = (x.z-dz)/tsize;\n"
"    float xi = (x.x-dx)/tsize;\n"
"    dbox3(vec3(dx, abs(x.y)-.1, dz), vec3(.48*tsize, .0005, .48*tsize), d);\n"
"    add(sdf, vec2(d, 4.), sdf);\n"
"    smoothmin(sdf.x, d, .002, sdf.x);\n"
"    indc = vec2(xi, zi);\n"
"\n"
"    // Logos\n"
"    tsize = .25;\n"
"    float tw = .0005;\n"
"    dz = mod(x.z-.5*tsize, tsize)-.5*tsize;\n"
"    zi = round((x.z-dz)/tsize);\n"
"    zi = mod(zi, 6.);\n"
"    \n"
"    if(zi < .5)dmercury(20.*x.xy, d);\n"
"    else if(zi < 1.5)dhaujobb(20.*x.xy, d);\n"
"    else if(zi < 2.5)dfarbrausch(20.*x.xy, d);\n"
"    else if(zi < 3.5)dkewlers(20.*x.xy, d);\n"
"	else if(zi < 4.5)dspacepigs(20.*x.xy, d);\n"
"	else if(zi < 5.5)dschnappsgirls(20.*x.xy, d);\n"
"    stroke(d/20.,tw, d);\n"
"    zextrude(dz, -d, .005, d);\n"
"    add(sdf, vec2(d, 5.), sdf);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{    \n"
"    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, \n"
"        s;\n"
" \n"
"    scale(iScale);\n"
"    //uv *= 2.;\n"
"    \n"
"    vec3 col = c.yyy, \n"
"        o = c.yyx,\n"
"        r = c.xyy, \n"
"        u = c.yxy, \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 400,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = 0.;//-(o.z-.1)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        d += min(s.x,3.e-2);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 5.e-4);\n"
"        vec3 l = mix(1.5,.2,abs(pow(sin(2.*2.*pi*(x.z-.05*iTime)), 2.)))*n;//normalize(x-.1*c.yxy);\n"
"       \n"
"		if(s.y == 2.)\n"
"        {\n"
"            col = .23*c.xxx;\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"        }\n"
"        else if(s.y == 3.)\n"
"        {\n"
"            float r;\n"
"            rand(ind, r);\n"
"            \n"
"            col = .2*vec3(0.02,0.11,0.24)\n"
"                + .2*vec3(0.25,0.75,0.85) * abs(dot(l,n))\n"
"                + mix(.5,.1,r)*vec3(0.45,0.69,0.76) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            \n"
"            // Reflections\n"
"            d = 1.e-2;\n"
"            o = x;\n"
"            dir = reflect(dir, n);\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                d += min(s.x,3.e-1);\n"
"                //d += s.x;\n"
"            }\n"
"            \n"
"            normal(x,n, 5.e-4);\n"
"        	vec3 l = mix(1.5,.2,abs(pow(sin(4.*2.*pi*(x.z-.05*iTime)), 2.)))*n;//normalize(x-.1*c.yxy);\n"
"       \n"
"            vec3 c1;\n"
"            if(s.y == 2.)\n"
"            {\n"
"                c1 = .23*c.xxx;\n"
"                c1 = .2*c1\n"
"                    + .2*c1 * abs(dot(l,n))\n"
"                    + .6*c1 * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            }\n"
"            else if(s.y == 3.)\n"
"            {\n"
"                float r;\n"
"                rand(ind, r);\n"
"\n"
"                c1 = .2*vec3(0.02,0.11,0.24)\n"
"                    + .2*vec3(0.25,0.75,0.85) * abs(dot(l,n))\n"
"                    + mix(.5,.1,r)*vec3(0.45,0.69,0.76) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            }\n"
"            else if(s.y == 4.)\n"
"            {\n"
"                float r;\n"
"                rand(indc, r);\n"
"\n"
"                c1 = .2*.2*c.xxx\n"
"                    + .2*.5*mix(c.xxx, vec3(0.02,0.11,0.24), step(0.,-x.y)) * abs(dot(l,n))\n"
"                    + mix(.5,.1,r)*.8*c.xxx * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            }\n"
"            else if(s.y == 5.)\n"
"            {\n"
"                c1 = .2*.2*c.xyy\n"
"                    + .2*.5*mix(c.xyy, vec3(0.24,0.11,0.024), step(0.,-x.y)) * abs(dot(l,n))\n"
"                    + .8*vec3(.8,.3,.2) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"                c1 = mix(c1, c.xxx, .1);\n"
"            }\n"
"            \n"
"            col = mix(col, c1, .5);\n"
"        }\n"
"        else if(s.y == 4.)\n"
"        {\n"
"            float r;\n"
"            rand(indc, r);\n"
"            \n"
"            col = .2*.2*c.xxx\n"
"                + .2*.5*mix(c.xxx, vec3(0.02,0.11,0.24), step(0.,-x.y)) * abs(dot(l,n))\n"
"                + mix(.5,.1,r)*.8*c.xxx * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"        }\n"
"        else if(s.y == 5.)\n"
"        {\n"
"            col = .2*.2*c.xyy\n"
"                + .2*.5*mix(c.xyy, vec3(0.24,0.11,0.024), step(0.,-x.y)) * abs(dot(l,n))\n"
"                + .8*vec3(.8,.3,.2) * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            col = mix(col, c.xxx, .1);\n"
"        }\n"
"\n"
"    }\n"
"    \n"
"    col = mix(col, 0.*.23*c.xxx*vec3(0.76,0.20,0.13),smoothstep(1.,5.,d));\n"
"    col *=3.6;\n"
"    //col = mix(col, 2.*col, iScale);\n"
"    \n"
"    //col = atan(col);\n"
"    col *= col;\n"
"    col = clamp(col, 0., 1.);\n"
"    col = mix(col, c.yyy, smoothstep(4.151818+3*1.8182,13.,iTime));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}	\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *evoke_source = "/* Lightcyber by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
"* Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"float iScale, nbeats;\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void dbox(in vec2 x, in vec2 b, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"void scale(out float s);\n"
"\n"
"void devoke(in vec2 x, out float d)\n"
"{\n"
"    x.x += .225;\n"
"    x *= 1.1;\n"
"    \n"
"    // o\n"
"    d = length(x+.35*c.xy)-.1;\n"
"    stroke(d,.06,d);\n"
"    \n"
"    // I\n"
"    float da;\n"
"    dbox(x+.1*c.xy, vec2(.05, .25), da);\n"
"    d = min(d, da);\n"
"    \n"
"    x = 2.*x - vec2(.4,-.2);\n"
"    // Mercury\n"
"    // Upper part\n"
"    dbox(x-.35*c.yx,vec2(.4,.35), da);\n"
"    d = min(d,da);\n"
"    dbox(x-.7*c.yx, vec2(.2,.2), da);\n"
"    d = max(d,-da);\n"
"    dbox(x-.25*c.yx,vec2(.2,.05),da);\n"
"    d = max(d,-da);\n"
"    \n"
"    // Lower part\n"
"    dbox(x+.1*c.yx,vec2(.1,.2),da);\n"
"    d = min(d,da);\n"
"    dbox(x+.2*c.yx, vec2(.4,.1),da);\n"
"    d = min(d,da);\n"
"    \n"
"    x = .5*(x + vec2(.4,-.2));\n"
"    \n"
"    // E\n"
"    // Right\n"
"    dbox(x-.9*c.xy, vec2(.05, .25), da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Top/bot\n"
"    dbox(vec2(x.x-.7, abs(x.y)-.2), vec2(.2, .05), da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Middle\n"
"    dbox(x-.7*c.xy, vec2(.2, .05), da);\n"
"    d = min(d,da);\n"
"    \n"
"    // Appendix\n"
"    dbox(vec2(x.x-.95,x.y+.2), vec2(.05,.05), da);\n"
"    d = min(d,da);\n"
"    \n"
"    stroke(d,.001,d);\n"
"}\n"
"\n"
"void dstripe(in vec2 x, out float d)\n"
"{\n"
"    dlinesegment(x-a*mix(-.4*c.xy, .4*c.xy, clamp(iTime/6.,0.,1.)), -.5*c.yx, .5*c.yx, d);\n"
"    d -= .005;\n"
"    float dd;\n"
"    vec2 vi;\n"
"    dvoronoi(5.*x, dd, vi);\n"
"    vi = x-vi/5.;\n"
"    dd = abs(length(vi)-.002)-.001;\n"
"    d = min(d,dd);\n"
"    stroke(d,.001,d);\n"
"    d = mix(1.,d, clamp(iTime, 0., 1.));\n"
"}\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), \n"
"        s;\n"
"    vec3 col = vec3(0.20,0.01,0.14), \n"
"        o = c.yyx,\n"
"        r = c.xyy, \n"
"        u = c.yxy,\n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    float d, i, ra;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    for(i=1.4; i>=0.; i -= .01)\n"
"    {\n"
"        lfnoise(102.*i*c.xx-mix(102.,0.,smoothstep(0.,1.,clamp(iTime-6.,0.,1.)))*iTime, ra);\n"
"        ra = .5+.5*ra;\n"
"        \n"
"		d = -(o.z-.2+i)/dir.z;\n"
"        x = o + d * dir;\n"
"        \n"
"        float da;\n"
"        dstripe(x.xy, da);\n"
"        devoke(x.xy, s.x);\n"
"        s.x = mix(da, s.x, smoothstep(0.,1., clamp(iTime-5.5,0.,1.)));\n"
"        s.x -= .01*iScale;\n"
"        \n"
"        if(ra < .5)\n"
"        {\n"
"            vec3 c1 = mix(mix(vec3(0.75,0.24,0.31), vec3(1.00,0.87,0.57), smoothstep(1.25,1.4,1.4-i)),vec3(0.20,0.01,0.14),i/1.4);\n"
"	        if(iTime > 6.)col = mix(col, c1, sm(s.x));\n"
"            col = mix(col, mix(1.1,1.,clamp(-iScale+smoothstep(0.,1.,clamp(iTime-6.,0.,1.)),0.,1.))*mix(col,vec3(.7,.45,.3), mix(.02,.1,iScale)), sm(s.x/64.));\n"
"        }\n"
"    }\n"
"    \n"
"    col = mix(col, c.yyy, clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));\n"
"    \n"
"    // Blend to voronoi background\n"
"    vec3 c1 = c.yyy;\n"
"    float v, v2;\n"
"    vec2 ind, ind2;\n"
"    lfnoise((iTime-12.)*c.xx, ind2.x);\n"
"    lfnoise((iTime-12.)*c.xx-1337., ind2.y);\n"
"    dvoronoi(12.*(uv-.03*ind2), v, ind);\n"
"    rand(ind, ra);\n"
"    stroke(-v, .05, v);\n"
"    v = -v;\n"
"    c1 = mix(c1, .3*ra*mix( .5*vec3(1.00,0.40,0.39), .05*c.xxx, clamp(tanh(1.5*length(uv)),0.,1.)), sm(v));\n"
"    \n"
"    c1 *= mix(1.,13., smoothstep(0.,1., clamp((iTime-11.), 0., 1.)));\n"
"    \n"
"    col = mix(col, c1, smoothstep(0.,1., clamp((iTime-11.),0.,1.)));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *canal_source = "/* Gross Gloss by Team210 - 64k intro by Team210 at Solskogen 2k19\n"
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
"\n"
"// Update1: Changes implementing FabriceNeyret2's comments.\n"
"\n"
"#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"const float pi = acos(-1.);\n"
"\n"
"void scale(out float s);\n"
"void rand(in vec2 x, out float n);\n"
"void hash31(in float p, out vec3 d);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void zextrude(in float z, in float d2d, in float h, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dsmoothvoronoi(in vec2 x, out float d, out vec2 z);\n"
"\n"
"vec2 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.z -= mix(1.3,-1.3,step(156., iTime))*iTime;\n"
"    \n"
"    float dx,\n"
"        d, v;\n"
"    \n"
"    lfnoise(.5*x.z*c.xx, dx);\n"
"    x.xy-=.4*dx*c.xy;\n"
"    \n"
"    // Voronoi\n"
"    float phi = atan(x.y,x.x);\n"
"    dsmoothvoronoi(2.*vec2(mod(phi+pi/4., 2.*pi), x.z), v, ind);\n"
"    stroke(v, .01, v);\n"
"    d = length(x.xy) - mix(1.,1.1, smoothstep(.0,.2,v));\n"
"    \n"
"    zextrude(length(x.xy)-1.0,d,.05, d);\n"
"    d -= .05;\n"
"    sdf = vec2(d,1.);\n"
"    \n"
"    // Smaller voronoi\n"
"    dsmoothvoronoi(8.*vec2(mod(phi+pi/4., 2.*pi), x.z), v, ind);\n"
"    stroke(v, .02, v);\n"
"    d = length(x.xy) - mix(1.1,1.2, smoothstep(.0,.2,v));\n"
"    \n"
"    zextrude(length(x.xy)-1.11,d,.01, d);\n"
"    d -= .1;\n"
"    add(sdf, vec2(d,1.), sdf);\n"
"    smoothmin(d, sdf.x , .1, sdf.x);\n"
"    \n"
"    float dy = 0.;\n"
"    //lfnoise(33.*x.x*c.xx-iTime,dy); \n"
"    \n"
"    for(float i=1.; i<=15.; i+=1.)\n"
"    {\n"
"        float f, a;\n"
"        vec2 dir;\n"
"        \n"
"        float n;\n"
"        \n"
"        lfnoise(i*c.xx-n, f);\n"
"        f = .5+.5*f;\n"
"        lfnoise(i*c.xx+1337.-n, a);\n"
"        a = .5+.5*a;\n"
"        lfnoise(i*c.xx+2337.-2.*n, dir.x);\n"
"        lfnoise(i*c.xx+3337.-3.*n, dir.y);\n"
"        dir = mix(c.yx-.2*c.xy, x.yx+.2*c.xy, 2.*dir);\n"
"        dir = normalize(dir);\n"
"        \n"
"        \n"
"        float dya = pow(1.01,f)* a * sin(-2.e-3*2.*pi*pow(1.95,abs(f+.01*a))*(1.*f-.01*a)*iTime-2.e-4*2.*pi*pow(1.99,abs(i-.1*a))*dot(dir,vec2(.5,4.)*(2.*(x.xz+1.3*(iTime))*c.yx)));\n"
"    	dy += 2.*pow((dya+1.)/2.,4.)-1.;\n"
"    }\n"
"    dy = .4*dy;\n"
"    \n"
"    add(sdf, vec2(x.y+.4+.001*dy, 2.), sdf);\n"
"    //smoothmin(x.y+.4-.02*dy, sdf.x , .5, sdf.x);\n"
"    \n"
"    \n"
"\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    vec2 uv = ( fragCoord -.5* iResolution.xy) / iResolution.y, \n"
"        s;\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    float dx, dx2, d0;\n"
"    lfnoise(-.5*1.3*iTime*c.xx, dx);\n"
"    lfnoise(-.5*1.3*(iTime+1.e-3)*c.xx, dx2);\n"
"\n"
"    vec3 col = c.yyy, \n"
"        o = c.yyx+.1*c.yxy+.4*dx*c.xyy,\n"
"        r = c.xyy, \n"
"        u = c.yxy, \n"
"        t = c.yyy+.4*dx2*c.xyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 400,\n"
"        i, a = 0;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"\n"
"    float d = .5/length(dir.xy);// -(o.z-.2)/dir.z;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"     	x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        if(length(x.xy-.4*dx*c.xy)>1.5)\n"
"        {\n"
"            col = c.yyy;\n"
"            i = N;\n"
"            break;\n"
"        }\n"
"        d += min(s.x,2.e-2);\n"
"        //d += s.x;\n"
"    }\n"
"    \n"
"    if(i < N)\n"
"    {\n"
"        normal(x,n, 5.e-4);\n"
"        vec3 l = normalize(x+.5*n);\n"
"       \n"
"		if(s.y == 1.)\n"
"        {\n"
"            col = mix(vec3(0.76,0.20,0.23), vec3(0.07,0.64,0.29), step(166., iTime));\n"
"            \n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            \n"
"            vec3 c1 = 2.*col;\n"
"            col = mix(col, c1, smoothstep(0.658, 1.02, 1.-abs(dot(n, c.yyz))));\n"
"            \n"
"            vec3 c2 = vec3(0.96,0.7,0.423);\n"
"        	col = mix(col, c2, smoothstep(0.658, 1.02, abs(dot(n, c.yyz))));\n"
"        }\n"
"        else if(s.y == 2.) // Mirror material\n"
"        {\n"
"            col = .3*c.xxx;\n"
"            col = .2*col\n"
"                + .2*col * abs(dot(l,n))\n"
"                + .3*col * pow(abs(dot(reflect(-l,n),dir)),2.);\n"
"            \n"
"            \n"
"            N = 50;\n"
"            o = x;\n"
"            dir = reflect(dir, n);\n"
"            d0 = d;\n"
"            d = 1.e-2;\n"
"            vec3 c1 = c.yyy;\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                if(length(x.xy)>1.5)\n"
"                {\n"
"                    c1 = c.yyy;\n"
"                    i = N;\n"
"                    break;\n"
"                }\n"
"                //d += min(s.x,5.e-3);\n"
"                d += s.x;\n"
"            }\n"
"            \n"
"            if(i < N)\n"
"            {\n"
"                normal(x,n, 5.e-4);\n"
"                vec3 l = normalize(x+.5*n);\n"
"\n"
"                if(s.y == 1.)\n"
"                {\n"
"                    c1 = vec3(0.76,0.20,0.23);\n"
"                    c1 = .2*c1\n"
"                        + .2*c1 * abs(dot(l,n))\n"
"                        + .6*c1 * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"                    c1 = mix(c1, 2.*vec3(0.76,0.20,0.13), smoothstep(0.658, 1.02, clamp(1.-abs(dot(n, c.yyz)),0.,1.)));\n"
"        			c1 = mix(c1, vec3(0.96,0.7,0.423), smoothstep(0.658, 1.02, abs(dot(n, c.yyz))));\n"
"                }\n"
"            }\n"
"            //c1 = mix(c1, c.yyy, smoothstep(3.,6.,d));\n"
"            \n"
"            col = mix(col, c1, .35);\n"
"            col = mix(col, .1*c.yxx, .3);\n"
"            //a = 1;\n"
"        }\n"
"\n"
"    }\n"
"    \n"
"    col = mix(col, c.yyy, smoothstep(2.,22.,d+d0));\n"
"    //col = mix(col, 2.*col, iScale);\n"
"    \n"
"    float nn;\n"
"    lfnoise(12.*(x.z-1.3*iTime)*c.xx, nn);\n"
"    \n"
"    col *=mix(1.1,2.6,mix(.5+.5*nn,1.,0.*iScale));\n"
"    col *= col;\n"
"    col = clamp(col, 0., 1.);\n"
"    if(col == c.xxx) col = c.yyy;\n"
"    \n"
"//     col = mix(col,c.yyy, .1);\n"
"    \n"
"//     col = mix(col, mix(col, length(col)/sqrt(3.)*c.xxx, .7), iScale);\n"
"//     col = mix(col, 1.7*1.7*col*col, iScale);\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *text_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
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
"#version 130\n\n"
"\n"
"uniform float iFontWidth, iTime;\n"
"uniform vec2 iResolution;\n"
"uniform sampler2D iChannel0, iFont;\n"
"uniform float iFSAA;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"float a; // Aspect ratio\n"
"\n"
"void rand(in vec2 x, out float num);\n"
"void lfnoise(in vec2 t, out float num);\n"
"void rshort(in float off, out float val);\n"
"void rfloat(in float off, out float val);\n"
"void dbox(in vec2 x, in vec2 b, out float dst);\n"
"void dcircle(in vec2 x, out float d);\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d);\n"
"void drhomboid(in vec2 x, in vec2 b, in float tilt, out float dst);\n"
"void dcirclesegment(in vec2 x, in float r, in float p0, in float p1, out float d);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst);\n"
"void dstring(in vec2 x, in float ordinal, in float size, out float dst);\n"
"void dfloat(in vec2 x, in float num, in float size, out float dst);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dint(in vec2 x, in float num, in float size, in float ndigits, out float dst);\n"
"void dtime(in vec2 x, in float num, in float size, out float dst);\n"
"\n"
"// Fixme: remove vec4 technique in favor of separate distance\n"
"// void blendadd(in vec4 src1, in vec4 src2, in float tlo, in float thi, out vec4 dst)\n"
"// {\n"
"//     vec4 added;\n"
"//     add(src1, src2, added);\n"
"//     dst = mix(src1, added, smoothstep(tlo-.5,tlo+.5,iTime)*(1.-smoothstep(thi-.5,thi+.5,iTime)));\n"
"// }\n"
"\n"
"void window(in vec2 x, in vec2 size, in vec3 bg, in float title_index, out vec4 col);\n"
"void progressbar(in vec2 x, in float width, in float progress, out vec4 col);\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
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
"    vec4 old = c.yyyy, \n"
"        new = c.yyyy;\n"
"    \n"
"    float bound = sqrt(iFSAA)-1.;\n"
"\n"
"    for(float i = -.5*bound; i<=.5*bound; i+=1.)\n"
"        for(float j=-.5*bound; j<=.5*bound; j+=1.)\n"
"        {\n"
"            old.gba += texture(iChannel0, (fragCoord+vec2(i,j)*3./max(bound, 1.))/iResolution.xy).xyz;\n"
"        }\n"
"    old.gba /= iFSAA;\n"
"    \n"
"    new = old;\n"
"    \n"
"//     if(uv.y < -.3 && iTime > 12.)\n"
"//     {\n"
"//         // Add overlay\n"
"//         colorize(2.*(c.xz*uv-.45*vec2(-a,1.)-12.*c.xy), new.gba);\n"
"//         new.gba = mix(old.gba, mix(old.gba, new.gba,.4), smoothstep(3.e-2, 5.e-2,length(new.gba)));\n"
"//     }\n"
"    \n"
"    if(uv.y > .4)\n"
"    {\n"
"        float da;\n"
"        dstring((uv-.45*vec2(-.55*a,1.+4.*.008)), 9., .004, d);\n"
"        dstring((uv-.45*vec2(-.55*a,1.+2.*.008)), 10., .004, da);\n"
"        d = min(d,da);\n"
"        dstring((uv-.45*vec2(-.55*a,1.)), 11., .004, da);\n"
"        d = min(d,da);\n"
"        dstring((uv-.45*vec2(-.55*a,1.-2.*.008)), 12., .004, da);\n"
"        d = min(d,da);\n"
"        dstring((uv-.45*vec2(-.55*a,1.-4.*.008)), 13., .004, da);\n"
"        d = min(d,da);\n"
"        new.gba = mix(new.gba, mix(new.gba, c.xxx, .5), sm(d));\n"
"        \n"
"        // Add Static text\n"
"        dstring((uv-.45*vec2(-.85*a,1.)), 3., .02, d); // Team210\n"
"        \n"
"        stroke(d-.002, .001, d);\n"
"        new.gba = mix(new.gba, vec3(1.00,0.40,0.39), sm(d));\n"
"\n"
"        // Add time overlay\n"
"        dtime((uv-.45*vec2(.975*a,1.05)), iTime+11., .01, d);\n"
"        new.gba = mix(new.gba, c.xxx, sm(d));\n"
"        \n"
"        // Add exact millisecond\n"
"        dint(uv-.45*vec2(.975*a,1.0), floor(1.e3*fract(iTime)), .01, 4., d);\n"
"//         new.gba = mix(new.gba, vec3(1.00,0.40,0.39), sm(d));\n"
"        stroke(d-.001, .0005, d);\n"
"        new.gba = mix(new.gba, c.xxx, sm(d));\n"
"    }\n"
"    \n"
"    if(iTime < 0.) \n"
"    {\n"
"        new.gba = old.gba;\n"
"        \n"
"        float sc = smoothstep(0.,1.,clamp(iTime+3.,0.,1.))*(1.-smoothstep(0.,1.,clamp(iTime+1.,0.,1.)));\n"
"        \n"
"        dstring((uv-vec2(-.085,-.3)), 3., .02, d); // Team210\n"
"        float da;\n"
"        dstring((uv-vec2(-.08,-.35)), 40., .02, da); // present\n"
"        d = min(d,da);\n"
"\n"
"        new.gba = mix(new.gba, mix(new.gba,c.yyy,sc), sm(d));\n"
"    }\n"
"    else if(iTime < 6.)\n"
"    {\n"
"        vec2 dx = (.25*a+.3*c.xy)*c.xy;\n"
"        if(iTime < 3.)\n"
"        {\n"
"            float ind = mix(100000., 2., clamp(iTime/3.,0.,1)), da;\n"
"            dint(uv+dx*c.xy, ind, .02, 6., d);\n"
"            \n"
"            dstring(uv+dx-2.*9.*.02*c.xy, 4., .02, da);\n"
"            d = min(d, da);\n"
"        }\n"
"        else if(iTime < 4.)\n"
"        {\n"
"            dint(uv+dx, 2., .02, 6., d);\n"
"            \n"
"            float da;\n"
"            dstring(uv+dx-2.*9.*.02*c.xy, 4., .02, da);\n"
"            d = min(d, da);\n"
"        }\n"
"        else if(iTime < 5.)\n"
"        {\n"
"            dint(uv+dx+.04*c.yx, 1., .02, 6., d);\n"
"         \n"
"            float da;\n"
"            dint(uv+dx, 2., .02, 6., da);\n"
"            d = min(d, da);\n"
"            \n"
"            dstring(uv+dx-2.*9.*.02*c.xy+.04*c.yx, 4., .02, da);\n"
"            d = min(d, da);\n"
"        }\n"
"        else if(iTime < 6.)\n"
"        {\n"
"            dint(uv+dx+.08*c.yx, 0., .02, 6., d);\n"
"            \n"
"            float da;\n"
"            dint(uv+dx+.04*c.yx, 1., .02, 6., da);\n"
"            d = min(d, da);\n"
"            \n"
"            dint(uv+dx, 2., .02, 6., da);\n"
"            d = min(d, da);\n"
"            \n"
"            dstring(uv+dx-2.*9.*.02*c.xy+.08*c.yx, 4., .02, da);\n"
"            d = min(d, da);\n"
"        }\n"
"        \n"
"            \n"
"        new.gba = mix(new.gba, mix(new.gba, vec3(1.00,0.87,0.57), .7), sm(d));\n"
"        stroke(d-.002, .001, d);\n"
"        new.gba = mix(new.gba, c.xxx, sm(d));\n"
"    }\n"
"    \n"
"    else if(iTime < 12. && iTime > 7.)\n"
"    {\n"
"        // EVK\n"
"        float da, db;\n"
"        dbox(vec2(uv.x+.75,uv.y-.35), vec2(.013,.035), da);\n"
"        stroke(da, .002, da);\n"
"        \n"
"        // E\n"
"        dglyph(vec2(uv.x+.75,uv.y-.35-.02).yx*c.zx, 101., .01, db);\n"
"        da = min(da, db);\n"
"        \n"
"        // V\n"
"        dglyph(vec2(uv.x+.75,uv.y-.35), 118., .01, db);\n"
"        da = min(da, db);\n"
"        \n"
"        // K\n"
"        dglyph(vec2(uv.x+.75,uv.y-.35+.02).yx*c.zx, 107., .01, db);\n"
"        da = min(da, db);\n"
"        \n"
"        // 333 block\n"
"        vec2 b = vec2(uv.x+.75,uv.y-.35+.02)-.01*c.xx-.02*c.xy,\n"
"            b1 = mod(b, .02)-.01,\n"
"            b1i = floor((b-b1)/.02);\n"
"        \n"
"        if(abs(b1i.y) <= 1. && b1i.x >= 0. && b1i.x <= 10.)\n"
"        {\n"
"            // random letter\n"
"            lfnoise(b1i-12.*iTime, db);\n"
"            db = 97.+mod(floor(26.*(.5+.5*db)),26.);\n"
"            dglyph(b1, db, .008, db);\n"
"            da = min(da, db);\n"
"        }\n"
"        \n"
"        dlinesegment(vec2(uv.x+.75,uv.y-.35+.06), -.015*c.xy, .25*c.xy, db);\n"
"        stroke(db, .001, db);\n"
"        da = min(da, db);\n"
"        \n"
"        // No more partycoding this time\n"
"        dstring(vec2(uv.x+.75,uv.y+.35), 5., .015, db);\n"
"        da = min(da, db);\n"
"        \n"
"        // Yeah. sure.\n"
"        dstring(vec2(uv.x-.2,uv.y+.35), 6., .015, db);\n"
"        float dc;\n"
"        dbox(vec2(uv.x-.2-.12,uv.y+.35), vec2(.165, .015), dc);\n"
"        db = max(dc,-db);\n"
"        da = min(da, db);\n"
"        \n"
"        // well, that worked.\n"
"        dstring(vec2(uv.x+.75,uv.y+.4),7., .015, db);\n"
"        da = min(da, db);\n"
"        \n"
"        // not\n"
"        dstring(vec2(uv.x-.2,uv.y+.4), 8., .015, db);\n"
"        dbox(vec2(uv.x-.2-.12,uv.y+.4), vec2(.165, .015), dc);\n"
"        db = max(dc,-db);\n"
"        da = min(da, db);\n"
"        \n"
"        \n"
"        new.gba = mix(new.gba, vec3(0.75,0.24,0.30), sm(da));\n"
"    }\n"
"    else if(iTime < 28.)\n"
"    {\n"
"        float da = length(uv)-.45, db;\n"
"        \n"
"        // Lightcyber\n"
"        dstring((uv+.3*c.xy), 2., .0415, db);\n"
"        db -= .001;\n"
"        da = max(da, -db);\n"
"        da = mix(1., da, smoothstep(0.,.5,clamp(iTime-18.5, 0., 1.))*(1.-smoothstep(0.,.5,clamp(iTime-22.,0.,1.))));\n"
"        new.gba = mix(new.gba, vec3(1.00,0.40,0.39), sm(da));\n"
"        \n"
"        // Team210\n"
"        da = length(uv - .3*c.xx)-.2, db;\n"
"        dstring(2.*(uv+.075*c.xy-.3*c.xx), 3., .0415, db);\n"
"        db -= .001;\n"
"//         da = max(da, -db);\n"
"        da = mix(1., da, smoothstep(0.,.5,clamp(iTime-19.5, 0., 1.))*(1.-smoothstep(0.,.5,clamp(iTime-22.,0.,1.))));\n"
"        db = mix(1., db, smoothstep(0.,.5,clamp(iTime-19.5, 0., 1.))*(1.-smoothstep(0.,.5,clamp(iTime-22.,0.,1.))));\n"
"        new.gba = mix(new.gba, vec3(1.00,0.40,0.39)*vec3(1.00,0.40,0.39), sm(da));\n"
"        new.gba = mix(new.gba, c.yyy, sm(db));\n"
"        \n"
"        // Graffiti tricks.\n"
"        dstring((uv-vec2(-.75,-.35)).yx*c.xz, 18., .045, db);\n"
"        dstring((uv-vec2(-.65,-.35)).yx*c.xz, 19., .045, da);\n"
"        db = min(db,da);\n"
"        db = mix(1., db, smoothstep(0.,.5,clamp(iTime-24.5, 0., 1.))*(1.-smoothstep(0.,.5,clamp(iTime-28.,0.,1.))));\n"
"        new.gba = mix(new.gba, mix(new.gba, c.xxx, .8), sm(db));\n"
"        \n"
"        stroke(db-.005, .0005, db);\n"
"        new.gba = mix(new.gba, mix(new.gba, vec3(1.00,0.40,0.39), .8), sm(db));\n"
"        \n"
"        // Nice!\n"
"        da = length((uv-vec2(-.6,-.325)).yx*c.xz)-.1;\n"
"        dstring((uv-vec2(-.6,-.35)).yx*c.xz, 20., .015, db);\n"
"        da = max(da, -db);\n"
"        da = mix(1., da, smoothstep(0.,.5,clamp(iTime-25.5, 0., 1.))*(1.-smoothstep(0.,.5,clamp(iTime-28.,0.,1.))));\n"
"        new.gba = mix(new.gba, mix(new.gba, c.xxx, .6), sm(da));\n"
"    }\n"
"    \n"
"    // \n"
"    float dc;\n"
"    dbox(uv, .5*vec2(a,1.), dc);\n"
"    stroke(dc, .005, dc);\n"
"    new.gba = mix(new.gba, c.yyy, sm(dc));\n"
"    \n"
"    fragColor = vec4(new.gba, 1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *post_source = "/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
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
"#version 130\n\n"
"\n"
"uniform float iFSAA;\n"
"uniform vec2 iResolution;\n"
"uniform sampler2D iChannel0;\n"
"uniform float iTime;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"float a = 1.0;\n"
"\n"
"float lscale, rscale;\n"
"float size;\n"
"\n"
"float nbeats;\n"
"float iScale;\n"
"\n"
"\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void dvoronoi(in vec2 x, out float d, out vec2 z);\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"float dot2( in vec3 v ) { return dot(v,v); }\n"
"\n"
"// Adapted from https://www.shadertoy.com/view/4sXXRN\n"
"void dtriangle3(in vec3 p,  in vec3 v1, in vec3 v2, in vec3 v3, out float dst)\n"
"{\n"
"    vec3 v21 = v2 - v1; vec3 p1 = p - v1;\n"
"    vec3 v32 = v3 - v2; vec3 p2 = p - v2;\n"
"    vec3 v13 = v1 - v3; vec3 p3 = p - v3;\n"
"    vec3 nor = cross( v21, v13 );\n"
"\n"
"    dst = sqrt( (sign(dot(cross(v21,nor),p1)) + \n"
"                  sign(dot(cross(v32,nor),p2)) + \n"
"                  sign(dot(cross(v13,nor),p3))<2.0) \n"
"                  ?\n"
"                  min( min( \n"
"                  dot2(v21*clamp(dot(v21,p1)/dot2(v21),0.0,1.0)-p1), \n"
"                  dot2(v32*clamp(dot(v32,p2)/dot2(v32),0.0,1.0)-p2) ), \n"
"                  dot2(v13*clamp(dot(v13,p3)/dot2(v13),0.0,1.0)-p3) )\n"
"                  :\n"
"                  dot(nor,p1)*dot(nor,p1)/dot2(nor) );\n"
"}\n"
"\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"mat3 R;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    float d;\n"
"    \n"
"	// Big red box    \n"
"    dbox3(x, .2*c.xxx, sdf.x);\n"
"    sdf.y = 1.;\n"
"    \n"
"    // Holes\n"
"    \n"
"    // 2 upper bar\n"
"    dbox3(x-.1*c.xyy, vec3(.02,.3,.12), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 right bar\n"
"    dbox3(x-.05*c.xyy-.1*c.yyx, vec3(.07,.3,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 mid bar\n"
"    dbox3(x, vec3(.02,.3,.1), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 left bar\n"
"    dbox3(x+.05*c.xyy+.1*c.yyx, vec3(.07,.3,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 dot\n"
"    dbox3(x+.1*c.xyy-.1*c.yyx, vec3(.02,.3,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 1 bar\n"
"    dbox3(x+.04*c.yyx, vec3(.3,.02,.08), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 1 dot\n"
"    dbox3(x-.1*c.yyx, vec3(.3,.02,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 0 big stripes\n"
"    vec3 y = vec3(x.x, abs(x.y), x.z);\n"
"    dbox3(y-.05*c.yxy, vec3(.1,.03,.3), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"\n"
"	// 0 small stripes\n"
"    dbox3(y-.1*c.yxy-.06*c.xyy, vec3(.08,.021,.3), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"\n"
"    // 0 upper/lower stripes\n"
"    vec3 z = vec3(abs(x.x), x.yz);\n"
"	dbox3(z-.119*c.xyy, vec3(.021,.08,.3), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"}\n"
"\n"
"void scene2(in vec3 x, out vec2 sdf)\n"
"{\n"
"    float v = 0.;\n"
"    vec2 vi = c.yy;\n"
"    dvoronoi(x.xy/size, v, vi);\n"
"    vec3 y = vec3(x.xy-vi*size, x.z);\n"
"    vec2 yi = vi*size;\n"
"    \n"
"    float n = 0.;\n"
"    lfnoise(4.*(yi-.5*iTime), n);\n"
"    lfnoise(12.*vec2(n,1.)*yi-(.8+.2*n)*c.xy, n);\n"
"    n *= iScale;\n"
"    //sdf = vec2(length(y-.05*n*c.yyx)-.5*size, 1.);\n"
"    sdf = vec2(length(y-.05*n*c.yyx)-mix(.05,1.,length(texture(iChannel0, yi/vec2(a,1.)).rgb)/sqrt(3.))*size, 1.);\n"
"}\n"
"\n"
"void normal2(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene2(x,s);\n"
"    scene2(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene2(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene2(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"void scene3(in vec3 x, out vec2 sdf)\n"
"{\n"
"    vec3 y = vec3(mod(x.xy,2.*size)-size, x.z);\n"
"    vec2 yi = x.xy-y.xy;\n"
"    float ss = mix(.0,.05,size/.01);\n"
"    \n"
"    vec2 p0 = .8*size*c.xx,\n"
"        p1 = .8*size*c.zx,\n"
"        p2 = .8*size*c.xz;\n"
"    \n"
"    vec2 ind;\n"
"    \n"
"    float y0, y1, y2;\n"
"    lfnoise(4.e1*(yi+p0-.5e-4*iTime), y0);\n"
"    lfnoise(12.e1*vec2(y0,1.)*(yi+p0)-1.e-4*(.8+.2*y0)*iTime*c.xy, y0);\n"
"    lfnoise(4.e1*(yi+p1-.5e-4*iTime), y1);\n"
"    lfnoise(12.e1*vec2(y1,1.)*(yi+p1)-1.e-4*(.8+.2*y1)*iTime*c.xy, y1);\n"
"    lfnoise(4.e1*(yi+p2-.5e-4*iTime), y2);\n"
"    lfnoise(12.e1*vec2(y2,1.)*(yi+p2)-1.e-4*(.8+.2*y2)*iTime*c.xy, y2);\n"
"    y0 *= ss;\n"
"    y1 *= ss;\n"
"    y2 *= ss;\n"
"    \n"
"    dtriangle3(y, vec3(p0,y0), vec3(p1,y1), vec3(p2,y2), sdf.x);\n"
"    \n"
"    float d;\n"
"    vec2 p3 = .8*size*c.zz,\n"
"        p4 = .8*size*c.xz,\n"
"        p5 = .8*size*c.zx;\n"
"    \n"
"    float y3, y4, y5;\n"
"    lfnoise(4.e1*(yi+p3-.5e-4*iTime), y3);\n"
"    lfnoise(12.e1*vec2(y3,1.)*(yi+p3)-1.e-4*(.8+.2*y3)*iTime*c.xy, y3);\n"
"    lfnoise(4.e1*(yi+p4-.5e-4*iTime), y4);\n"
"    lfnoise(12.e1*vec2(y4,1.)*(yi+p4)-1.e-4*(.8+.2*y4)*iTime*c.xy, y4);\n"
"    lfnoise(4.e1*(yi+p5-.5e-4*iTime), y5);\n"
"    lfnoise(12.e1*vec2(y5,1.)*(yi+p5)-1.e-4*(.8+.2*y5)*iTime*c.xy, y5);\n"
"    y3 *= ss;\n"
"    y4 *= ss;\n"
"    y5 *= ss;\n"
"    \n"
"    dtriangle3(y, vec3(p3,y3), vec3(p4,y4), vec3(p5,y5), d);\n"
"    sdf.x = min(sdf.x, d);\n"
"\n"
"    stroke(sdf.x, .1*size, sdf.x);\n"
"    sdf.y = 1.;\n"
"}\n"
"\n"
"void normal3(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene3(x,s);\n"
"    scene3(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene3(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene3(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord_ )\n"
"{\n"
"    vec2 fragCoord = fragCoord_;\n"
"    float a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    \n"
"    nbeats = mod(iTime, 60./29.);\n"
"    iScale = nbeats-30./29.;\n"
"    iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0., 15./29., iScale));\n"
"    \n"
"    vec3 col = texture(iChannel0, fragCoord_/iResolution).rgb;\n"
"    float delta = 0.;\n"
"//     vec2 n = c.yy;\n"
"    \n"
"    // Box\n"
"    rot3(vec3(-2.*pi/8.,2.*pi/8.,2.*pi/4.)-iTime*vec3(1.1,1.3,1.5), R);\n"
"    \n"
"    float d;\n"
"    vec2 s;\n"
"    vec3 o, r, u, t, ssize, dir, x, n;\n"
"    vec2 uv2 = 10.*(uv-vec2(-.45*a,.45));\n"
"    o = R * c.yyx;\n"
"	r = c.xyy; \n"
"	u = c.yxy;\n"
"	t = c.yyy; \n"
"    int N = 250,\n"
"        i;\n"
"    t = uv2.x * r + uv2.y * u;\n"
"    t = R * t;\n"
"    dir = normalize(t-o);\n"
"\n"
"    ssize = .2*c.xxx;\n"
"\n"
"	vec3 tlo = min((ssize-o)/dir,(-ssize-o)/dir); // Select 3 visible planes\n"
"    vec2 abxlo = abs(o.yz + tlo.x*dir.yz),\n"
"        abylo = abs(o.xz + tlo.y*dir.xz),\n"
"        abzlo = abs(o.xy + tlo.z*dir.xy);\n"
"\n"
"    vec4 dn = 100.*c.xyyy;\n"
"\n"
"    dn = mix(dn, vec4(tlo.x,c.xyy), float(all(lessThan(abxlo,ssize.yz)))*step(tlo.x,dn.x));\n"
"    dn = mix(dn, vec4(tlo.y,c.yxy), float(all(lessThan(abylo,ssize.xz)))*step(tlo.y,dn.x));\n"
"    dn = mix(dn, vec4(tlo.z,c.yyx), float(all(lessThan(abzlo,ssize.xy)))*step(tlo.z,dn.x));\n"
"    \n"
"    uv = (fragCoord)/iResolution.xy*vec2(a,1.);\n"
"    \n"
"    d = dn.r;\n"
"    \n"
"    float nan;\n"
"    lfnoise(iTime*c.xx, nan);\n"
"    nan += .5;\n"
"    if(nan > 0.) d = 3.;\n"
"    \n"
"    if(d<=2.)\n"
"    {\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        \n"
"        if(s.x > 1.e-4)\n"
"        {\n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                d += s.x;\n"
"            }\n"
"        }\n"
"        \n"
"        if(i<N)\n"
"        {\n"
"            normal(x,n, 5.e-4);\n"
"            \n"
"            if(s.y == 1.)\n"
"            {\n"
"                vec3 l = normalize(x+c.zzx*vec3(1.3,.9,1.2));\n"
"                col = vec3(0.81,0.15,0.18);\n"
"                col = .3*col\n"
"                    + .4*col * abs(dot(l,n))\n"
"                    + .6 * col * abs(pow(dot(reflect(-l,n),dir),2.));\n"
"            }\n"
"            else if(s.y == 2.)\n"
"            {\n"
"                vec3 l = normalize(x+c.zzx*vec3(1.3,.9,1.2));\n"
"                col = .7*c.xxx;\n"
"                col = .5*col\n"
"                    + .4*col * abs(dot(l,n))\n"
"                    + .8 * col * abs(pow(dot(reflect(-l,n),dir),2.));\n"
"            }\n"
"        }\n"
"        \n"
"        if(iTime < 0.) col = texture(iChannel0, fragCoord_/iResolution).rgb;\n"
"    }\n"
"    else\n"
"    {\n"
"        iScale = nbeats-30./29.;\n"
"        iScale = smoothstep(-5./29., 0., iScale)*(1.-smoothstep(0./29., 35./29., iScale));\n"
"//         lscale = iScale;\n"
"        lscale = smoothstep(0.,.5,clamp((iTime-10.),0.,1.))*(1.-smoothstep(0.,.5,clamp((iTime-18.),0.,1.)));\n"
"//         lscale += smoothstep(0.,.5,clamp((iTime-10.),0.,1.))*(1.-smoothstep(0.,.5,clamp((iTime-18.),0.,1.)));\n"
"        \n"
"        \n"
"//         rscale = iScale;\n"
"        rscale = 0.;\n"
"        size = mix(.005, .01, rscale);\n"
"        size = mix(0., size, max(rscale, lscale));\n"
"     \n"
"        if(lscale > 0.)\n"
"        {\n"
"            col = c.yyy;\n"
"            \n"
"            o = c.yyx+.5*vec3(cos(iTime), sin(iTime),0.);\n"
"            r = c.xyy;\n"
"            u = c.yxy;\n"
"            t = c.yyy;\n"
"            dir = c.yyy;\n"
"            n = c.yyy;\n"
"            x = c.yyy;\n"
"            N = 200;\n"
"            t = uv.x * r + uv.y * u;\n"
"            dir = normalize(t-o);\n"
"\n"
"            d = -(o.z-.05-.5*size)/dir.z;\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene2(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                \n"
"                if(x.z<-.05-.5*size)\n"
"                {\n"
"                    col = c.yyy;\n"
"                    i = N;\n"
"                    break;\n"
"                }\n"
"                d += min(s.x,1.e-3);\n"
"                //d += s.x;\n"
"            }\n"
"            \n"
"            if(i < N)\n"
"            {\n"
"                normal2(x,n, 5.e-4);\n"
"                vec3 l = normalize(x+.5*n);\n"
"            \n"
"                if(s.y == 1.)\n"
"                {\n"
"                    float v;\n"
"                    vec2 vi;\n"
"                    dvoronoi(x.xy/size, v, vi);\n"
"                    vec3 y = vec3(x.xy-vi*size, x.z);\n"
"                    vec2 yi = vi*size;\n"
"                    \n"
"                    float bound = sqrt(iFSAA)-1.;\n"
"\n"
"                    for(float i = -.5*bound; i<=.5*bound; i+=1.)\n"
"                        for(float j=-.5*bound; j<=.5*bound; j+=1.)\n"
"                        {\n"
"                            col += texture(iChannel0, yi/vec2(a,1.)+vec2(i,j)*3./max(bound, 1.)/iResolution.xy).xyz;\n"
"                        }\n"
"                    col /= iFSAA;   \n"
"                    \n"
"                    col = .4*col\n"
"                        + .9*col * abs(dot(l,n))\n"
"                        + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"                }\n"
"            }\n"
"            else col = c.yyy;\n"
"        }\n"
"        else if(rscale > 0.)\n"
"        {\n"
"            col = c.yyy;\n"
"            \n"
"            o = c.yyx+.5*vec3(-1., -1.,0.);\n"
"            r = c.xyy;\n"
"            u = c.yxy;\n"
"            t = c.yyy;\n"
"            dir = c.yyy;\n"
"            n = c.yyy;\n"
"            x = c.yyy;\n"
"            N = 300;\n"
"            t = uv.x * r + uv.y * u;\n"
"            dir = normalize(t-o);\n"
"\n"
"            d = -(o.z-.05-.5*size)/dir.z;\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene3(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                \n"
"                if(x.z<-.05-.5*size)\n"
"                {\n"
"                    col = c.yyy;\n"
"                    i = N;\n"
"                    break;\n"
"                }\n"
"                d += min(s.x,1.e-3);\n"
"                //d += s.x;\n"
"            }\n"
"            \n"
"            if(i < N)\n"
"            {\n"
"                normal3(x,n, 5.e-4);\n"
"                vec3 l = normalize(x+.5*n);\n"
"            \n"
"                if(s.y == 1.)\n"
"                {\n"
"                    vec3 y = vec3(mod(x.xy,size)-.5*size, x.z);\n"
"                    vec2 yi = x.xy-y.xy;\n"
"                    \n"
"                    col = texture(iChannel0, yi/vec2(a,1.)).rgb;\n"
"                    \n"
"//                     col = .7*c.xxy;\n"
"                    \n"
"                    col = .4*col\n"
"                        + .9*col * abs(dot(l,n))\n"
"                        + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"                    \n"
"                }\n"
"            }\n"
"            else col = c.yyy;\n"
"        }\n"
"    }\n"
"    \n"
"    // Scan lines\n"
"    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *logo210_source = "/* Fuer Elite - 64k Intro by Team210 at Underground Conference 9\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" * \n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" * \n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" * \n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
" \n"
"#version 130\n \n"
" \n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"float a = 1.0;\n"
"\n"
"float nbeats, iScale;\n"
"\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void dbox210(in vec3 x, in float size, out vec2 sdf)\n"
"{\n"
"    x /= size;\n"
"    \n"
"    float d = 1.;\n"
"    \n"
"    // Big red box    \n"
"    dbox3(x, .2*c.xxx, sdf.x);\n"
"    sdf.y = 1.;\n"
"    \n"
"    // Holes\n"
"    \n"
"    // 2 upper bar\n"
"    dbox3(x-.1*c.xyy, vec3(.02,.3,.12), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 right bar\n"
"    dbox3(x-.05*c.xyy-.1*c.yyx, vec3(.07,.3,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 mid bar\n"
"    dbox3(x, vec3(.02,.3,.1), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 left bar\n"
"    dbox3(x+.05*c.xyy+.1*c.yyx, vec3(.07,.3,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 2 dot\n"
"    dbox3(x+.1*c.xyy-.1*c.yyx, vec3(.02,.3,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 1 bar\n"
"    dbox3(x+.04*c.yyx, vec3(.3,.02,.08), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 1 dot\n"
"    dbox3(x-.1*c.yyx, vec3(.3,.02,.02), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    // 0 big stripes\n"
"    vec3 y = vec3(x.x, abs(x.y), x.z);\n"
"    dbox3(y-.05*c.yxy, vec3(.1,.03,.3), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"\n"
"	// 0 small stripes\n"
"    dbox3(y-.1*c.yxy-.06*c.xyy, vec3(.08,.021,.3), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"\n"
"    // 0 upper/lower stripes\n"
"    vec3 z = vec3(abs(x.x), x.yz);\n"
"	dbox3(z-.119*c.xyy, vec3(.021,.08,.3), d);\n"
"    sdf.x = max(-d, sdf.x);\n"
"    sdf.y = mix(sdf.y, 2., step(d, sdf.x));\n"
"    \n"
"    sdf.x *= size;\n"
"}\n"
"\n"
"mat3 R;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    sdf = c.xy;\n"
"    \n"
"    float d, da;\n"
"    \n"
"    rot3(vec3(-pi/2.,0.,pi/2.), R);\n"
"    x = R * x;\n"
"    \n"
"    vec2 sda = c.xy;\n"
"    \n"
"	dbox210(x+.1*c.xyy, .5, sdf);\n"
"	rot3(vec3(pi/2.,0.,pi/2.), R);\n"
"    x = R * x;\n"
"    //add(sdf,sda,sdf);\n"
"    \n"
"    dbox210(x,5.,sda);\n"
"    add(sdf,sda,sdf);\n"
"    \n"
"    rot3(vec3(pi/2.,-pi/2.,pi/2.), R);\n"
"    x = R * x;\n"
"    \n"
"    dbox210(x-2.*c.yxy,50.,sda);\n"
"    add(sdf,sda,sdf);\n"
"    \n"
"    stroke(sdf.x,.001, sdf.x);\n"
"    \n"
"    dbox3(x, 100.*c.xxx, sda.x);\n"
"    sda.y = 2.;\n"
"    \n"
"    add(sdf, sda*c.zx, sdf);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    float a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    \n"
"    vec3 col = c.xxx;\n"
"    \n"
"    float d = 0.;\n"
"    vec2 s;\n"
"    vec3 o, t, dir, x, n;\n"
"    \n"
"    mat3 Ra;\n"
"    rot3(mix(c.yyy,vec3(-5.7884463,2.4242211,0.3463173),clamp((iTime-6.)/1.5,0.,1.)), Ra);\n"
"    //vec3 a = vec3(uv,-1.);\n"
"	\n"
"    //rot3(mix(c.yyy,vec3(-3.*pi/4.,3.*pi/4.,-7.*pi/4.),clamp((iTime-6.)/1.5,0.,1.)), Ra);\n"
"    //Ra *= mix(1.,-1.,clamp((iTime-6.)/1.5,0.,1.));\n"
"       \n"
"    \n"
"    o = Ra * mix(mix(mix(c.yyy-.1*c.yxy,c.yyx,clamp(iTime/2.,0.,1.)),10.*c.yyx,clamp((iTime-2.)/2.,0.,1.)), 100.*c.yyx, clamp((iTime-4.)/2.,0.,1.));\n"
"	t = c.yyy;\n"
"    int N = 650,\n"
"        i;\n"
"    dir = Ra *normalize(vec3(uv,-1.));//normalize(t-o);\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        d += s.x;\n"
"    }\n"
"        \n"
"    if(s.x < 1.e-4)\n"
"    {\n"
"        normal(x,n, 5.e-4);\n"
"        vec3 l = normalize(x+.1*n);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            col = vec3(0.81,0.15,0.18);\n"
"            col = .3*col\n"
"                + .4*col * abs(dot(l,n))\n"
"                + .6 * col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            col = .7*c.xxx;\n"
"            col = .5*col\n"
"                + .4*col * abs(dot(l,n))\n"
"                + .8 * col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"            \n"
"            vec3 c1 = c.yyy;\n"
"            \n"
"            o = x;\n"
"            dir = reflect(dir,n);\n"
"            d = 1.e-1;\n"
"            \n"
"            N = 150;\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                d += s.x;\n"
"            }\n"
"                \n"
"            if(s.x < 1.e-4)\n"
"            {\n"
"                normal(x,n, 5.e-4);\n"
"                vec3 l = normalize(x+.1*n);\n"
"\n"
"                if(s.y == 1.)\n"
"                {\n"
"                    c1 = vec3(0.81,0.15,0.18);\n"
"                    c1 = .3*c1\n"
"                        + .4*c1 * abs(dot(l,n))\n"
"                        + .6 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                }\n"
"                else if(s.y == 2.)\n"
"                {\n"
"                    c1 = .7*c.xxx;\n"
"                    c1 = .5*c1\n"
"                        + .4*c1 * abs(dot(l,n))\n"
"                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                }\n"
"                c1 = clamp(c1, 0., 1.);\n"
"                \n"
"                col = mix(col, c1, .2);\n"
"            }\n"
"            \n"
"            col = clamp(col, 0., 1.);\n"
"        }\n"
"        \n"
"    }\n"
"    col = mix(col,vec3(0.20,0.01,0.14),smoothstep(0.,1.,iTime-10.));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *transbubbles_source = "#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"float a = 1.0;\n"
"\n"
"float nbeats, iScale;\n"
"\n"
"void scale(out float s);\n"
"void rand(in vec2 x, out float n);\n"
"void lfnoise(in vec2 t, out float n);\n"
"void dbox3(in vec3 x, in vec3 b, out float d);\n"
"void rot3(in vec3 p, out mat3 rot);\n"
"void stroke(in float d0, in float s, out float d);\n"
"void hash13(in vec3 p3, out float d);\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf);\n"
"void smoothmin(in float a, in float b, in float k, out float dst);\n"
"void dvoronoi3(in vec3 x, out float d, out vec3 z)\n"
"{\n"
"    vec3 y = floor(x);\n"
"    float ret = 1.;\n"
"    vec3 pf=c.yyy, p;\n"
"    float df=10.;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            for(int k=-1; k<=1; k+=1)\n"
"            {\n"
"                p = y + vec3(float(i), float(j), float(k));\n"
"                float pa;\n"
"                hash13(p, pa);\n"
"                p += pa;\n"
"\n"
"                d = length(x-p);\n"
"\n"
"                if(d < df)\n"
"                {\n"
"                    df = d;\n"
"                    pf = p;\n"
"                }\n"
"            }\n"
"        }\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            for(int k=-1; k<=1; k+=1)\n"
"            {\n"
"                p = y + vec3(float(i), float(j), float(k));\n"
"                float pa;\n"
"                hash13(p, pa);\n"
"                p += pa;\n"
"\n"
"                vec3 o = p - pf;\n"
"                d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);\n"
"                ret = min(ret, d);\n"
"            }\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    z = pf;\n"
"}\n"
"\n"
"mat3 R;\n"
"vec3 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x = R * x;\n"
"    x -= mix(0.,.1*iTime, step(150., iTime));\n"
"//     x.z -= .1*iTime;\n"
"    \n"
"    sdf = c.xy;\n"
"    \n"
"    float bsize = 2.;\n"
"    float d, da;\n"
"    \n"
"    dvoronoi3(bsize*x, d, ind);\n"
"    vec3 y = x-ind/bsize;\n"
"    \n"
"    float n;\n"
"    hash13(ind,n);\n"
"    \n"
"    //add(sdf, vec2(abs(length(y)-.3)-.001,2.), sdf);\n"
"	add(sdf, vec2(abs(length(y)-.3)-.001,2.), sdf);\n"
"	sdf.x = max(sdf.x, -length(x)+.5);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx);\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    float a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    rot3(.1*vec3(1.1,1.3,1.5)*iTime, R);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    float d = 0.;\n"
"    vec2 s;\n"
"    vec3 o, t, dir, x, n;\n"
"    \n"
"    //mat3 Ra;\n"
"    //rot3(mix(c.yyy,vec3(-3.*pi/4.,3.*pi/4.,-7.*pi/4.),clamp((iTime-6.)/1.5,0.,1.)), Ra);\n"
"    //Ra *= mix(1.,-1.,clamp((iTime-6.)/1.5,0.,1.));\n"
"       \n"
"    float mx = clamp((iTime-5.),0.,1.), \n"
"        my = clamp(iTime-10., 0., 1.);\n"
"    float nai;\n"
"    lfnoise(2.*iTime*c.xx, nai);\n"
"    nai = .5*(nai);\n"
"    //o = Ra * mix(mix(mix(c.yyy-.1*c.yxy,c.yyx,clamp(iTime/2.,0.,1.)),10.*c.yyx,clamp((iTime-2.)/2.,0.,1.)), 100.*c.yyx, clamp((iTime-4.)/2.,0.,1.));\n"
"	o = c.yyx;\n"
"    t = c.yyy;\n"
"    int N = 250,\n"
"        i;\n"
"    dir = normalize(vec3(uv,-1.));//normalize(t-o);\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        d += s.x;\n"
"        //d += min(s.x, .1);\n"
"    }\n"
"        \n"
"    if(s.x < 1.e-4)\n"
"    {\n"
"        normal(x,n, 1.e-4);\n"
"        vec3 l = normalize(x+.1*n);\n"
"        \n"
"        if(s.y == 2.)\n"
"        {\n"
"            col = mix(vec3(0.86,0.21,0.13), vec3(0.02,0.46,0.44), mx);\n"
"            col = .1*col\n"
"                            + .1*col * abs(dot(l,n))\n"
"                            + .5 * col * abs(pow(dot(reflect(-l,n),dir),2.));\n"
"            \n"
"            vec3 c1 = c.yyy;\n"
"            for(float fraction = 0.; fraction <= 3.; fraction += 1.)\n"
"    		{\n"
"//                 o = x;\n"
"//                 dir = refract(dir,n,.9);\n"
"//                 //dir = reflect(dir,n);\n"
"//                 d = 1.e-2;\n"
"                o = x;\n"
"                vec3 ddir = refract(dir,n,.95+.05*nai);\n"
"                vec3 dddir = refract(dir,n,.1);\n"
"                dir = refract(dir,n,.5+.05*nai);//reflect(dir,n);\n"
"                \n"
"                dir = mix(ddir, dir,mx);\n"
"                dir = mix(dir,dddir,my);\n"
"                dir = mix(dir, ddir, step(150.,iTime));\n"
"                \n"
"                d = 2.e-2;\n"
"\n"
"                for(i = 0; i<N; ++i)\n"
"                {\n"
"                    x = o + d * dir;\n"
"                    scene(x,s);\n"
"                    if(s.x < 1.e-4)break;\n"
"                    d += s.x;\n"
"                    //d += min(s.x, .01);\n"
"                }\n"
"                \n"
"//                 if((R*x).z<-1.)\n"
"//                 {\n"
"//                     fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"//                     return;\n"
"//                 }\n"
"                \n"
"                if(s.x < 1.e-4)\n"
"                {\n"
"                    normal(x,n, 1.e-4);\n"
"                    vec3 l = normalize(x+.1*n);\n"
"\n"
"                    if(s.y == 2.)\n"
"                    {\n"
"                        //c1 = .7*c.xxx;\n"
"                        c1 = (fraction == 0.)?vec3(0.86,0.21,0.13):\n"
"                        	(fraction == 1.)?vec3(0.85,0.80,0.62):\n"
"                        	(fraction == 2.)?vec3(0.22,0.25,0.25):\n"
"                        	(fraction == 3.)?vec3(0.16,0.17,0.17):\n"
"                        vec3(0.12,0.12,0.13); // sieht ganz ok aus\n"
"//                         vec3 c2 = (fraction == 4.)?vec3(0.20,0.15,0.20):\n"
"//                         	(fraction == 3.)?vec3(0.40,0.30,0.32):\n"
"//                         	(fraction == 2.)?vec3(0.97,0.48,0.33):\n"
"//                         	(fraction == 1.)?vec3(1.00,0.59,0.31):\n"
"//                         vec3(0.65,0.60,0.53);\n"
"                        vec3 c2 = (fraction == 0.)?vec3(0.99,0.33,0.05):\n"
"                        	(fraction == 1.)?vec3(0.94,0.94,0.94):\n"
"                        	(fraction == 2.)?vec3(0.75,0.82,0.88):\n"
"                        	(fraction == 3.)?vec3(0.25,0.34,0.39):\n"
"                        vec3(0.17,0.22,0.27);\n"
"                        vec3 c3 = (fraction == 0.)?vec3(1.00,0.55,0.03):\n"
"                        	(fraction == 1.)?vec3(0.84,0.20,0.18):\n"
"                        	(fraction == 2.)?vec3(0.13,0.55,0.57):\n"
"                        	(fraction == 3.)?vec3(0.29,0.22,0.30):\n"
"                        vec3(0.00,0.00,0.00);\n"
"                        vec3 c4 = (fraction == 0.)?vec3(0.09,0.00,0.13):\n"
"                        	(fraction == 1.)?vec3(0.87,0.38,0.61):\n"
"                        	(fraction == 2.)?vec3(0.97,0.45,0.50):\n"
"                        	(fraction == 3.)?vec3(0.97,0.70,0.72):\n"
"                        vec3(0.04,0.04,0.05);\n"
"                        vec3 c5 = (fraction == 0.)?vec3(0.94,0.24,0.27):\n"
"                        	(fraction == 1.)?vec3(0.48,0.49,0.55):\n"
"                        	(fraction == 2.)?vec3(0.81,0.76,0.90):\n"
"                        	(fraction == 3.)?vec3(0.44,0.33,0.60):\n"
"                        vec3(0.34,0.24,0.50);\n"
"                        vec3 c6 = (fraction == 0.)?vec3(0.11,0.32,0.23):\n"
"                        	(fraction == 1.)?vec3(0.25,0.89,0.79):\n"
"                        	(fraction == 2.)?vec3(0.31,0.59,0.56):\n"
"                        	(fraction == 3.)?vec3(0.00,0.34,0.39):\n"
"                        vec3(0.00,0.15,0.17);\n"
"                        c1 = mix(c1,c2, mx);\n"
"                        c1 = mix(c1,c3, my);\n"
"                        c1 = mix(c1, c4, step(150.,iTime));\n"
"                        c1 = mix(c1, c5, step(163.5,iTime));\n"
"                        c1 = mix(c1, c6, step(170.,iTime));\n"
"\n"
"                        c1 = .1*c1\n"
"                            + .4*c1 * abs(dot(l,n))\n"
"                            + mix(mix(5.8,3.79,mx),4.753,my) * c1 * abs(pow(dot(reflect(-l,n),dir),2.));\n"
"                        c1 = mix(c1, 2.*c1, smoothstep(mix(1.,.6,iScale), 1.02, 1.-abs(dot(n, c.xyy))));\n"
"                        c1 = mix(c1, 2.*c1, smoothstep(mix(1.,.6,iScale), 1.02, abs(dot(n, c.zyy))));\n"
"                    }//5.8\n"
"					//col = clamp(col, 0., 1.);\n"
"                    col = mix(col, c1, .15);\n"
"                }\n"
"\n"
"                col = clamp(col, 0., 1.);\n"
"            }\n"
"        }\n"
"        \n"
"    }\n"
"    col *= col;\n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *volclouds_source = "#version 130\n\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"float a = 1.0;\n"
"\n"
"float nbeats, iScale;\n"
"\n"
"void scale(out float s);\n"
"// Creative Commons Attribution-ShareAlike 4.0 International Public License\n"
"// Created by David Hoskins.\n"
"// See https://www.shadertoy.com/view/4djSRW\n"
"void hash13(in vec3 p3, out float d)\n"
"{\n"
"	p3  = fract(p3 * .1031);\n"
"    p3 += dot(p3, p3.yzx + 33.33);\n"
"    d = fract((p3.x + p3.y) * p3.z);\n"
"}\n"
"\n"
"// Arbitrary-frequency 2D noise\n"
"void lfnoise3(in vec3 t, out float num)\n"
"{\n"
"    t -= vec3(11.,13.,5.);\n"
"    vec3 i = floor(t);\n"
"    t = fract(t);\n"
"    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise\n"
"    t = smoothstep(c.yyy, c.xxx, t); // TODO: add this for faster value noise\n"
"    vec2 v1, v2, v3, v4;\n"
"    hash13(i, v1.x);\n"
"    hash13(i+c.xyy, v1.y);\n"
"    hash13(i+c.yxy, v2.x);\n"
"    hash13(i+c.xxy, v2.y);\n"
"    hash13(i+c.yyx, v3.x);\n"
"    hash13(i+c.xyx, v3.y);\n"
"    hash13(i+c.yxx, v4.x);\n"
"    hash13(i+c.xxx, v4.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    v3 = c.zz+2.*mix(v3, v4, t.y);\n"
"    v2.x = mix(v1.x, v1.y, t.x);\n"
"    v2.y = mix(v3.x, v3.y, t.x);\n"
"    num = mix(v2.x, v2.y, t.z);\n"
"}\n"
"\n"
"void mfnoise3(in vec3 x, in float d, in float b, in float e, out float n)\n"
"{\n"
"    n = 0.;\n"
"    float a = 1., nf = 0., buf;\n"
"    for(float f = d; f<b; f *= 2.)\n"
"    {\n"
"        lfnoise3(f*x-vec3(11.,13.,5.), buf);\n"
"        n += a*buf;\n"
"        a *= e;\n"
"        nf += 1.;\n"
"    }\n"
"    n *= (1.-e)/(1.-pow(e, nf));\n"
"}\n"
"\n"
"// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"mat3 R;\n"
"vec3 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x = R * x;\n"
"    \n"
"    float n;\n"
"    mfnoise3(x-.1*iTime*c.yyx,10.,400.,.25,n);\n"
"    n = .5+.5*n;\n"
"    sdf = vec2(-n, 2.);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"void palette1(in float scale, out vec3 col)\n"
"{\n"
"    const int N = 5;\n"
"   \n"
"    /*\n"
"    const vec3 colors[N] = vec3[N](\n"
"            vec3(0.82,0.27,0.13),\n"
"            vec3(0.85,0.77,0.68),\n"
"            vec3(0.65,0.59,0.55),\n"
"            vec3(0.45,0.29,0.24),\n"
"            vec3(0.85,0.27,0.15)\n"
"        );\n"
"    //*/\n"
"    \n"
"    /*\n"
"	const vec3 colors[N] = vec3[N](\n"
"       	vec3(0.86,0.21,0.13),\n"
"        vec3(0.85,0.80,0.62),\n"
"        vec3(0.22,0.25,0.25),\n"
"        vec3(0.16,0.17,0.17),\n"
"        vec3(0.12,0.12,0.13)\n"
"    );\n"
"    //*/\n"
"    \n"
"	//*\n"
"    vec3 colors[N];\n"
"    if(iTime < 150.)\n"
"        colors = vec3[N](\n"
"       	vec3(0.99,0.33,0.05),\n"
"        vec3(0.94,0.94,0.94),\n"
"        vec3(0.75,0.82,0.88),\n"
"        vec3(0.25,0.34,0.39),\n"
"        vec3(0.17,0.22,0.27));\n"
"    else if(iTime < 160.)\n"
"        colors = vec3[N](\n"
"            vec3(0.82,0.27,0.13),\n"
"            vec3(0.85,0.77,0.68),\n"
"            vec3(0.65,0.59,0.55),\n"
"            vec3(0.45,0.29,0.24),\n"
"            vec3(0.85,0.27,0.15)\n"
"        );\n"
"    else if(iTime < 165.)\n"
"        colors = vec3[N](\n"
"            vec3(0.11,0.32,0.23),\n"
"            vec3(0.25,0.89,0.79),\n"
"            vec3(0.31,0.59,0.56),\n"
"            vec3(0.00,0.34,0.39),\n"
"            vec3(0.00,0.15,0.17)\n"
"        );\n"
"    else colors = vec3[N](\n"
"       	vec3(0.86,0.21,0.13),\n"
"        vec3(0.85,0.80,0.62),\n"
"        vec3(0.22,0.25,0.25),\n"
"        vec3(0.16,0.17,0.17),\n"
"        vec3(0.12,0.12,0.13)\n"
"    );\n"
"    //*/\n"
"	float index = floor(scale*float(N)), \n"
"        remainder = scale*float(N)-index;\n"
"    col = mix(colors[int(index)],colors[int(index)+1], remainder);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    rot3(.01*vec3(1.1,1.3,1.5)*iTime, R);\n"
"    \n"
"    scale(iScale);\n"
"    \n"
"    float a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    float d = 0.;\n"
"    vec2 s;\n"
"    vec3 o, t, dir, x, n;\n"
"    \n"
"	o = c.yyx;\n"
"    t = c.yyy;\n"
"    int N = 50,\n"
"        i;\n"
"    dir = normalize(vec3(uv,-1.));//normalize(t-o);\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"        d += .5/float(N);\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        normal(x,n,5.e-4);\n"
"        vec3 l = normalize(x+.1*n);\n"
"        vec3 c1;\n"
"        palette1(-s.x, c1);\n"
"        c1 = .1*c1\n"
"                            + .1*c1 * abs(dot(l,n))\n"
"                            + 1.354 * c1 * abs(pow(dot(reflect(-l,n),dir),2.));\n"
"    	c1 = mix(c1, 2.*c1, smoothstep(mix(1.,.6,iScale), 1.02, 1.-abs(dot(n, c.xyy))));\n"
"        c1 = mix(c1, 2.*c1, smoothstep(mix(1.,.6,iScale), 1.02, abs(dot(n, c.zyy))));\n"
"    	col = mix(col, c1, d*d);\n"
"    	\n"
"    }\n"
"\n"
"    col *= col;\n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
const char *chart_source = "#version 130\n \n"
" \n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"uniform float iFader0;\n"
"uniform float iFader1;\n"
"uniform float iFader2;\n"
"uniform float iFader3;\n"
"uniform float iFader4;\n"
"uniform float iFader5;\n"
"uniform float iFader6;\n"
"uniform float iFader7;\n"
"\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"float a = 1.0;\n"
"\n"
"float nbeats, iScale;\n"
"\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"void lfnoise(in vec2 t, out float n)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    t = smoothstep(c.yy, c.xx, t);\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    n = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\n"
"void mfnoise(in vec2 x, in float d, in float b, in float e, out float n)\n"
"{\n"
"    n = 0.;\n"
"    float a = 1., nf = 0., buf;\n"
"    for(float f = d; f<b; f *= 2.)\n"
"    {\n"
"        lfnoise(f*x, buf);\n"
"        n += a*buf;\n"
"        a *= e;\n"
"        nf += 1.;\n"
"    }\n"
"    n *= (1.-e)/(1.-pow(e, nf));\n"
"}\n"
"\n"
"void dbox3(in vec3 x, in vec3 b, out float d)\n"
"{\n"
"  vec3 da = abs(x) - b;\n"
"  d = length(max(da,0.0))\n"
"         + min(max(da.x,max(da.y,da.z)),0.0);\n"
"}\n"
"\n"
"void rot3(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\n"
"// Add sdfs\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = sda.x<sdb.x?sda:sdb;\n"
"}\n"
"\n"
"void rotAB(in vec3 a, in vec3 b, out mat3 R)\n"
"{\n"
"    a = normalize(a);\n"
"    b = normalize(b);\n"
"    vec3 v = cross(a,b);\n"
"    float co = dot(a,b);\n"
"    R = mat3(0.,v.z,-v.y,-v.z,0.,v.x,v.y,-v.x,0.);\n"
"    R = R*R/(1.+co) + R;\n"
"    R += mat3(1.);\n"
"}\n"
"\n"
"mat3 R;\n"
"float ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    sdf = 3.*c.yx;\n"
"    \n"
"    float d, da;\n"
"    dbox3(x, 2.*c.xxx, sdf.x);\n"
"    sdf.x *= -1.;\n"
"	sdf.x = min(sdf.x, x.z);\n"
"    \n"
"    vec3 y = x;\n"
"    x = vec3(x.x, mod(x.y,.2)-.1, x.z);\n"
"    vec3 xi = (y-x)/.2;\n"
"    ind = xi.y;\n"
"    \n"
"    float n;\n"
"    mfnoise((x.x+.3*iTime)*c.xx-xi.y,12.,120.,.45, n);\n"
"    n *= clamp(.2-.01*xi.y,0.,1.);\n"
"    \n"
"    vec2 sda;\n"
"    dbox3(x-.2*xi.y*c.xyy, vec3(.5,.01,.3-n), sda.x);\n"
"    sda.y = 1.;\n"
"    add(sdf,sda,sdf);\n"
"    \n"
"    dbox3(x-(.3-n)*c.yyx-.2*xi.y*c.xyy, vec3(.5,.03,.01), sda.x);\n"
"    sda.y = 4.;\n"
"    add(sdf, sda, sdf);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n, in float dx)\n"
"{\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"float sm(float d)\n"
"{\n"
"    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);\n"
"}\n"
"\n"
"void palette1(in float scale, out vec3 col)\n"
"{\n"
"    const int N = 5;\n"
"   \n"
"    /*\n"
"    const vec3 colors[N] = vec3[N](\n"
"            vec3(0.82,0.27,0.13),\n"
"            vec3(0.85,0.77,0.68),\n"
"            vec3(0.65,0.59,0.55),\n"
"            vec3(0.45,0.29,0.24),\n"
"            vec3(0.85,0.27,0.15)\n"
"        );\n"
"    //*/\n"
"    \n"
"    /*\n"
"	const vec3 colors[N] = vec3[N](\n"
"       	vec3(0.86,0.21,0.13),\n"
"        vec3(0.85,0.80,0.62),\n"
"        vec3(0.22,0.25,0.25),\n"
"        vec3(0.16,0.17,0.17),\n"
"        vec3(0.12,0.12,0.13)\n"
"    );\n"
"    //*/\n"
"    \n"
"	//*\n"
"    const vec3 colors[N] = vec3[N](\n"
"       	vec3(0.99,0.33,0.05),\n"
"        vec3(0.94,0.94,0.94),\n"
"        vec3(0.75,0.82,0.88),\n"
"        vec3(0.25,0.34,0.39),\n"
"        vec3(0.17,0.22,0.27)\n"
"    );\n"
"    //*/\n"
"	float index = floor(scale*float(N)), \n"
"        remainder = scale*float(N)-index;\n"
"    col = mix(colors[int(index)],colors[int(index)+1], remainder);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    float a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"\n"
"    vec3 col = c.yyy, \n"
"        o = c.yzx,\n"
"        r = c.xyy, \n"
"        u = normalize(c.yxx), \n"
"        t = c.yyy, \n"
"        dir,\n"
"        n,\n"
"        x;\n"
"    int N = 400,\n"
"        i;\n"
"    t = uv.x * r + uv.y * u;\n"
"    dir = normalize(t-o);\n"
"    \n"
"    float d = 0.;\n"
"    vec2 s;\n"
"    \n"
"    for(i = 0; i<N; ++i)\n"
"    {\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4)break;\n"
"        //d += s.x;\n"
"        d += min(s.x, 8.e-3);\n"
"    }\n"
"        \n"
"    if(s.x < 1.e-4)\n"
"    {\n"
"        normal(x,n, 5.e-4);\n"
"        vec3 l = normalize(x+.1*n);\n"
"        \n"
"        if(s.y == 1.)\n"
"        {\n"
"            col = mix(vec3(0.81,0.15,0.18),vec3(0.62,0.27,0.35),clamp(1.-x.z/.3,0.,1.));\n"
"            col = .3*col\n"
"                + .4*col * abs(dot(l,n))\n"
"                + .9 * col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"        }\n"
"        else if(s.y == 4.)\n"
"        {\n"
"        	col = c.xxx;\n"
"            col = .5*col\n"
"            	+ .4*col * abs(dot(l,n))\n"
"                + .8 * col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"        }\n"
"        else if(s.y == 2. || s.y == 3.)\n"
"        {\n"
"            col = c.xxx;\n"
"            vec2 ma = abs(mod(x.xy,.05)-.025)-.0015;\n"
"            col = mix(col,.5*c.xxx, sm(min(ma.x,ma.y)));\n"
"            col = .5*col\n"
"                + .4*col * abs(dot(l,n))\n"
"                + .8 * col * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"            \n"
"            vec3 c1 = c.yyy;\n"
"            \n"
"            o = x;\n"
"            dir = reflect(dir,n);\n"
"            d = 1.e-2;\n"
"            \n"
"            N = 150;\n"
"            \n"
"            for(i = 0; i<N; ++i)\n"
"            {\n"
"                x = o + d * dir;\n"
"                scene(x,s);\n"
"                if(s.x < 1.e-4)break;\n"
"                d += s.x;\n"
"            }\n"
"                \n"
"            if(s.x < 1.e-4)\n"
"            {\n"
"                normal(x,n, 5.e-4);\n"
"                vec3 l = normalize(x+.1*n);\n"
"\n"
"                if(s.y == 1.)\n"
"                {\n"
"                    c1 = mix(vec3(0.81,0.15,0.18),vec3(0.62,0.27,0.35),clamp(1.-x.z/.3,0.,1.));\n"
"                    c1 = .3*c1\n"
"                        + .4*c1 * abs(dot(l,n))\n"
"                        + .9 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                }\n"
"                else if(s.y == 2.)\n"
"                {\n"
"                    c1 = .7*c.xxx;\n"
"                    c1 = .5*c1\n"
"                        + .4*c1 * abs(dot(l,n))\n"
"                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                    \n"
"                }\n"
"                else if(s.y == 3.)\n"
"                {\n"
"                    c1 = .7*c.xxx;\n"
"                    c1 = .5*c1\n"
"                        + .4*c1 * abs(dot(l,n))\n"
"                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                }\n"
"                else if(s.y == 4.)\n"
"                {\n"
"                    c1 = c.xxx;\n"
"                    c1 = .5*c1\n"
"                        + .4*c1 * abs(dot(l,n))\n"
"                        + .8 * c1 * abs(pow(dot(reflect(-l,n),dir),3.));\n"
"                }\n"
"                \n"
"                c1 = clamp(c1, 0., 1.);\n"
"                \n"
"                col = mix(col, c1, .2);\n"
"            }\n"
"            \n"
"            col = clamp(col, 0., 1.);\n"
"        }\n"
"     	\n"
"    }\n"
"    //col = mix(col,vec3(0.20,0.01,0.14),smoothstep(0.,1.,iTime-10.));\n"
"    \n"
"    fragColor = vec4(clamp(col,0.,1.),1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\0";
void Loadscale()
{
    int scale_size = strlen(scale_source);
    scale_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(scale_handle, 1, (GLchar **)&scale_source, &scale_size);
    glCompileShader(scale_handle);
#ifdef DEBUG
    printf("---> scale Shader:\n");
    debug(scale_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddsmoothvoronoi()
{
    int dsmoothvoronoi_size = strlen(dsmoothvoronoi_source);
    dsmoothvoronoi_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dsmoothvoronoi_handle, 1, (GLchar **)&dsmoothvoronoi_source, &dsmoothvoronoi_size);
    glCompileShader(dsmoothvoronoi_handle);
#ifdef DEBUG
    printf("---> dsmoothvoronoi Shader:\n");
    debug(dsmoothvoronoi_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrand()
{
    int rand_size = strlen(rand_source);
    rand_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rand_handle, 1, (GLchar **)&rand_source, &rand_size);
    glCompileShader(rand_handle);
#ifdef DEBUG
    printf("---> rand Shader:\n");
    debug(rand_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadhash31()
{
    int hash31_size = strlen(hash31_source);
    hash31_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(hash31_handle, 1, (GLchar **)&hash31_source, &hash31_size);
    glCompileShader(hash31_handle);
#ifdef DEBUG
    printf("---> hash31 Shader:\n");
    debug(hash31_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadlfnoise()
{
    int lfnoise_size = strlen(lfnoise_source);
    lfnoise_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(lfnoise_handle, 1, (GLchar **)&lfnoise_source, &lfnoise_size);
    glCompileShader(lfnoise_handle);
#ifdef DEBUG
    printf("---> lfnoise Shader:\n");
    debug(lfnoise_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadmfnoise()
{
    int mfnoise_size = strlen(mfnoise_source);
    mfnoise_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(mfnoise_handle, 1, (GLchar **)&mfnoise_source, &mfnoise_size);
    glCompileShader(mfnoise_handle);
#ifdef DEBUG
    printf("---> mfnoise Shader:\n");
    debug(mfnoise_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddbox()
{
    int dbox_size = strlen(dbox_source);
    dbox_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dbox_handle, 1, (GLchar **)&dbox_source, &dbox_size);
    glCompileShader(dbox_handle);
#ifdef DEBUG
    printf("---> dbox Shader:\n");
    debug(dbox_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddlinesegment3()
{
    int dlinesegment3_size = strlen(dlinesegment3_source);
    dlinesegment3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dlinesegment3_handle, 1, (GLchar **)&dlinesegment3_source, &dlinesegment3_size);
    glCompileShader(dlinesegment3_handle);
#ifdef DEBUG
    printf("---> dlinesegment3 Shader:\n");
    debug(dlinesegment3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadstroke()
{
    int stroke_size = strlen(stroke_source);
    stroke_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(stroke_handle, 1, (GLchar **)&stroke_source, &stroke_size);
    glCompileShader(stroke_handle);
#ifdef DEBUG
    printf("---> stroke Shader:\n");
    debug(stroke_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadzextrude()
{
    int zextrude_size = strlen(zextrude_source);
    zextrude_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(zextrude_handle, 1, (GLchar **)&zextrude_source, &zextrude_size);
    glCompileShader(zextrude_handle);
#ifdef DEBUG
    printf("---> zextrude Shader:\n");
    debug(zextrude_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadadd()
{
    int add_size = strlen(add_source);
    add_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(add_handle, 1, (GLchar **)&add_source, &add_size);
    glCompileShader(add_handle);
#ifdef DEBUG
    printf("---> add Shader:\n");
    debug(add_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadsmoothmin()
{
    int smoothmin_size = strlen(smoothmin_source);
    smoothmin_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(smoothmin_handle, 1, (GLchar **)&smoothmin_source, &smoothmin_size);
    glCompileShader(smoothmin_handle);
#ifdef DEBUG
    printf("---> smoothmin Shader:\n");
    debug(smoothmin_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddspline3()
{
    int dspline3_size = strlen(dspline3_source);
    dspline3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dspline3_handle, 1, (GLchar **)&dspline3_source, &dspline3_size);
    glCompileShader(dspline3_handle);
#ifdef DEBUG
    printf("---> dspline3 Shader:\n");
    debug(dspline3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddvoronoi()
{
    int dvoronoi_size = strlen(dvoronoi_source);
    dvoronoi_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dvoronoi_handle, 1, (GLchar **)&dvoronoi_source, &dvoronoi_size);
    glCompileShader(dvoronoi_handle);
#ifdef DEBUG
    printf("---> dvoronoi Shader:\n");
    debug(dvoronoi_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadnormal()
{
    int normal_size = strlen(normal_source);
    normal_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(normal_handle, 1, (GLchar **)&normal_source, &normal_size);
    glCompileShader(normal_handle);
#ifdef DEBUG
    printf("---> normal Shader:\n");
    debug(normal_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddbox3()
{
    int dbox3_size = strlen(dbox3_source);
    dbox3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dbox3_handle, 1, (GLchar **)&dbox3_source, &dbox3_size);
    glCompileShader(dbox3_handle);
#ifdef DEBUG
    printf("---> dbox3 Shader:\n");
    debug(dbox3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrot3()
{
    int rot3_size = strlen(rot3_source);
    rot3_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rot3_handle, 1, (GLchar **)&rot3_source, &rot3_size);
    glCompileShader(rot3_handle);
#ifdef DEBUG
    printf("---> rot3 Shader:\n");
    debug(rot3_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddtriangle()
{
    int dtriangle_size = strlen(dtriangle_source);
    dtriangle_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dtriangle_handle, 1, (GLchar **)&dtriangle_source, &dtriangle_size);
    glCompileShader(dtriangle_handle);
#ifdef DEBUG
    printf("---> dtriangle Shader:\n");
    debug(dtriangle_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddlinesegment()
{
    int dlinesegment_size = strlen(dlinesegment_source);
    dlinesegment_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dlinesegment_handle, 1, (GLchar **)&dlinesegment_source, &dlinesegment_size);
    glCompileShader(dlinesegment_handle);
#ifdef DEBUG
    printf("---> dlinesegment Shader:\n");
    debug(dlinesegment_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddpolygon()
{
    int dpolygon_size = strlen(dpolygon_source);
    dpolygon_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dpolygon_handle, 1, (GLchar **)&dpolygon_source, &dpolygon_size);
    glCompileShader(dpolygon_handle);
#ifdef DEBUG
    printf("---> dpolygon Shader:\n");
    debug(dpolygon_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrot()
{
    int rot_size = strlen(rot_source);
    rot_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rot_handle, 1, (GLchar **)&rot_source, &rot_size);
    glCompileShader(rot_handle);
#ifdef DEBUG
    printf("---> rot Shader:\n");
    debug(rot_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddcircle()
{
    int dcircle_size = strlen(dcircle_source);
    dcircle_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dcircle_handle, 1, (GLchar **)&dcircle_source, &dcircle_size);
    glCompileShader(dcircle_handle);
#ifdef DEBUG
    printf("---> dcircle Shader:\n");
    debug(dcircle_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddschnappsgirls()
{
    int dschnappsgirls_size = strlen(dschnappsgirls_source);
    dschnappsgirls_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dschnappsgirls_handle, 1, (GLchar **)&dschnappsgirls_source, &dschnappsgirls_size);
    glCompileShader(dschnappsgirls_handle);
#ifdef DEBUG
    printf("---> dschnappsgirls Shader:\n");
    debug(dschnappsgirls_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddspacepigs()
{
    int dspacepigs_size = strlen(dspacepigs_source);
    dspacepigs_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dspacepigs_handle, 1, (GLchar **)&dspacepigs_source, &dspacepigs_size);
    glCompileShader(dspacepigs_handle);
#ifdef DEBUG
    printf("---> dspacepigs Shader:\n");
    debug(dspacepigs_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddkewlers()
{
    int dkewlers_size = strlen(dkewlers_source);
    dkewlers_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dkewlers_handle, 1, (GLchar **)&dkewlers_source, &dkewlers_size);
    glCompileShader(dkewlers_handle);
#ifdef DEBUG
    printf("---> dkewlers Shader:\n");
    debug(dkewlers_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddfarbrausch()
{
    int dfarbrausch_size = strlen(dfarbrausch_source);
    dfarbrausch_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dfarbrausch_handle, 1, (GLchar **)&dfarbrausch_source, &dfarbrausch_size);
    glCompileShader(dfarbrausch_handle);
#ifdef DEBUG
    printf("---> dfarbrausch Shader:\n");
    debug(dfarbrausch_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddhaujobb()
{
    int dhaujobb_size = strlen(dhaujobb_source);
    dhaujobb_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dhaujobb_handle, 1, (GLchar **)&dhaujobb_source, &dhaujobb_size);
    glCompileShader(dhaujobb_handle);
#ifdef DEBUG
    printf("---> dhaujobb Shader:\n");
    debug(dhaujobb_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddmercury()
{
    int dmercury_size = strlen(dmercury_source);
    dmercury_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dmercury_handle, 1, (GLchar **)&dmercury_source, &dmercury_size);
    glCompileShader(dmercury_handle);
#ifdef DEBUG
    printf("---> dmercury Shader:\n");
    debug(dmercury_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrshort()
{
    int rshort_size = strlen(rshort_source);
    rshort_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rshort_handle, 1, (GLchar **)&rshort_source, &rshort_size);
    glCompileShader(rshort_handle);
#ifdef DEBUG
    printf("---> rshort Shader:\n");
    debug(rshort_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadrfloat()
{
    int rfloat_size = strlen(rfloat_source);
    rfloat_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(rfloat_handle, 1, (GLchar **)&rfloat_source, &rfloat_size);
    glCompileShader(rfloat_handle);
#ifdef DEBUG
    printf("---> rfloat Shader:\n");
    debug(rfloat_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddrhomboid()
{
    int drhomboid_size = strlen(drhomboid_source);
    drhomboid_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(drhomboid_handle, 1, (GLchar **)&drhomboid_source, &drhomboid_size);
    glCompileShader(drhomboid_handle);
#ifdef DEBUG
    printf("---> drhomboid Shader:\n");
    debug(drhomboid_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddcirclesegment()
{
    int dcirclesegment_size = strlen(dcirclesegment_source);
    dcirclesegment_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dcirclesegment_handle, 1, (GLchar **)&dcirclesegment_source, &dcirclesegment_size);
    glCompileShader(dcirclesegment_handle);
#ifdef DEBUG
    printf("---> dcirclesegment Shader:\n");
    debug(dcirclesegment_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddglyph()
{
    int dglyph_size = strlen(dglyph_source);
    dglyph_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dglyph_handle, 1, (GLchar **)&dglyph_source, &dglyph_size);
    glCompileShader(dglyph_handle);
#ifdef DEBUG
    printf("---> dglyph Shader:\n");
    debug(dglyph_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddstring()
{
    int dstring_size = strlen(dstring_source);
    dstring_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dstring_handle, 1, (GLchar **)&dstring_source, &dstring_size);
    glCompileShader(dstring_handle);
#ifdef DEBUG
    printf("---> dstring Shader:\n");
    debug(dstring_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddfloat()
{
    int dfloat_size = strlen(dfloat_source);
    dfloat_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dfloat_handle, 1, (GLchar **)&dfloat_source, &dfloat_size);
    glCompileShader(dfloat_handle);
#ifdef DEBUG
    printf("---> dfloat Shader:\n");
    debug(dfloat_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddint()
{
    int dint_size = strlen(dint_source);
    dint_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dint_handle, 1, (GLchar **)&dint_source, &dint_size);
    glCompileShader(dint_handle);
#ifdef DEBUG
    printf("---> dint Shader:\n");
    debug(dint_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loaddtime()
{
    int dtime_size = strlen(dtime_source);
    dtime_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(dtime_handle, 1, (GLchar **)&dtime_source, &dtime_size);
    glCompileShader(dtime_handle);
#ifdef DEBUG
    printf("---> dtime Shader:\n");
    debug(dtime_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadwindow()
{
    int window_size = strlen(window_source);
    window_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(window_handle, 1, (GLchar **)&window_source, &window_size);
    glCompileShader(window_handle);
#ifdef DEBUG
    printf("---> window Shader:\n");
    debug(window_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadprogressbar()
{
    int progressbar_size = strlen(progressbar_source);
    progressbar_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(progressbar_handle, 1, (GLchar **)&progressbar_source, &progressbar_size);
    glCompileShader(progressbar_handle);
#ifdef DEBUG
    printf("---> progressbar Shader:\n");
    debug(progressbar_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}
void Loadhash13()
{
    int hash13_size = strlen(hash13_source);
    hash13_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(hash13_handle, 1, (GLchar **)&hash13_source, &hash13_size);
    glCompileShader(hash13_handle);
#ifdef DEBUG
    printf("---> hash13 Shader:\n");
    debug(hash13_handle);
    printf(">>>>\n");
#endif
    progress += .2/(float)nsymbols;
}

void LoadSymbols()
{
    Loadscale();
    updateBar();
    Loaddsmoothvoronoi();
    updateBar();
    Loadrand();
    updateBar();
    Loadhash31();
    updateBar();
    Loadlfnoise();
    updateBar();
    Loadmfnoise();
    updateBar();
    Loaddbox();
    updateBar();
    Loaddlinesegment3();
    updateBar();
    Loadstroke();
    updateBar();
    Loadzextrude();
    updateBar();
    Loadadd();
    updateBar();
    Loadsmoothmin();
    updateBar();
    Loaddspline3();
    updateBar();
    Loaddvoronoi();
    updateBar();
    Loadnormal();
    updateBar();
    Loaddbox3();
    updateBar();
    Loadrot3();
    updateBar();
    Loaddtriangle();
    updateBar();
    Loaddlinesegment();
    updateBar();
    Loaddpolygon();
    updateBar();
    Loadrot();
    updateBar();
    Loaddcircle();
    updateBar();
    Loaddschnappsgirls();
    updateBar();
    Loaddspacepigs();
    updateBar();
    Loaddkewlers();
    updateBar();
    Loaddfarbrausch();
    updateBar();
    Loaddhaujobb();
    updateBar();
    Loaddmercury();
    updateBar();
    Loadrshort();
    updateBar();
    Loadrfloat();
    updateBar();
    Loaddrhomboid();
    updateBar();
    Loaddcirclesegment();
    updateBar();
    Loaddglyph();
    updateBar();
    Loaddstring();
    updateBar();
    Loaddfloat();
    updateBar();
    Loaddint();
    updateBar();
    Loaddtime();
    updateBar();
    Loadwindow();
    updateBar();
    Loadprogressbar();
    updateBar();
    Loadhash13();
    updateBar();
}
int voronoidesign_program, voronoidesign_handle, groundboxes_program, groundboxes_handle, graffiti_program, graffiti_handle, greet_program, greet_handle, evoke_program, evoke_handle, canal_program, canal_handle, text_program, text_handle, post_program, post_handle, logo210_program, logo210_handle, transbubbles_program, transbubbles_handle, volclouds_program, volclouds_handle, chart_program, chart_handle;
int voronoidesign_iTime_location,voronoidesign_iResolution_location,voronoidesign_iFader0_location,voronoidesign_iFader1_location,voronoidesign_iFader2_location,voronoidesign_iFader3_location,voronoidesign_iFader4_location,voronoidesign_iFader5_location,voronoidesign_iFader6_location,voronoidesign_iFader7_location;
int groundboxes_iTime_location,groundboxes_iResolution_location,groundboxes_iFader0_location,groundboxes_iFader1_location,groundboxes_iFader2_location,groundboxes_iFader3_location,groundboxes_iFader4_location,groundboxes_iFader5_location,groundboxes_iFader6_location,groundboxes_iFader7_location;
int graffiti_iTime_location,graffiti_iResolution_location,graffiti_iFader0_location,graffiti_iFader1_location,graffiti_iFader2_location,graffiti_iFader3_location,graffiti_iFader4_location,graffiti_iFader5_location,graffiti_iFader6_location,graffiti_iFader7_location;
int greet_iTime_location,greet_iResolution_location,greet_iFader0_location,greet_iFader1_location,greet_iFader2_location,greet_iFader3_location,greet_iFader4_location,greet_iFader5_location,greet_iFader6_location,greet_iFader7_location;
int evoke_iTime_location,evoke_iResolution_location,evoke_iFader0_location,evoke_iFader1_location,evoke_iFader2_location,evoke_iFader3_location,evoke_iFader4_location,evoke_iFader5_location,evoke_iFader6_location,evoke_iFader7_location;
int canal_iTime_location,canal_iResolution_location,canal_iFader0_location,canal_iFader1_location,canal_iFader2_location,canal_iFader3_location,canal_iFader4_location,canal_iFader5_location,canal_iFader6_location,canal_iFader7_location;
int text_iFontWidth_location,text_iTime_location,text_iResolution_location,text_iChannel0_location,text_iFont_location,text_iFSAA_location;
int post_iFSAA_location,post_iResolution_location,post_iChannel0_location,post_iTime_location;
int logo210_iTime_location,logo210_iResolution_location,logo210_iFader0_location,logo210_iFader1_location,logo210_iFader2_location,logo210_iFader3_location,logo210_iFader4_location,logo210_iFader5_location,logo210_iFader6_location,logo210_iFader7_location;
int transbubbles_iTime_location,transbubbles_iResolution_location,transbubbles_iFader0_location,transbubbles_iFader1_location,transbubbles_iFader2_location,transbubbles_iFader3_location,transbubbles_iFader4_location,transbubbles_iFader5_location,transbubbles_iFader6_location,transbubbles_iFader7_location;
int volclouds_iTime_location,volclouds_iResolution_location,volclouds_iFader0_location,volclouds_iFader1_location,volclouds_iFader2_location,volclouds_iFader3_location,volclouds_iFader4_location,volclouds_iFader5_location,volclouds_iFader6_location,volclouds_iFader7_location;
int chart_iTime_location,chart_iResolution_location,chart_iFader0_location,chart_iFader1_location,chart_iFader2_location,chart_iFader3_location,chart_iFader4_location,chart_iFader5_location,chart_iFader6_location,chart_iFader7_location;
const int nprograms = 12;

void Loadvoronoidesign()
{
    int voronoidesign_size = strlen(voronoidesign_source);
    voronoidesign_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(voronoidesign_handle, 1, (GLchar **)&voronoidesign_source, &voronoidesign_size);
    glCompileShader(voronoidesign_handle);
#ifdef DEBUG
    printf("---> voronoidesign Shader:\n");
    debug(voronoidesign_handle);
    printf(">>>>\n");
#endif
    voronoidesign_program = glCreateProgram();
    glAttachShader(voronoidesign_program,voronoidesign_handle);
    glAttachShader(voronoidesign_program,scale_handle);
    glAttachShader(voronoidesign_program,dsmoothvoronoi_handle);
    glAttachShader(voronoidesign_program,rand_handle);
    glAttachShader(voronoidesign_program,hash31_handle);
    glAttachShader(voronoidesign_program,lfnoise_handle);
    glAttachShader(voronoidesign_program,mfnoise_handle);
    glAttachShader(voronoidesign_program,dbox_handle);
    glAttachShader(voronoidesign_program,dlinesegment3_handle);
    glAttachShader(voronoidesign_program,stroke_handle);
    glAttachShader(voronoidesign_program,zextrude_handle);
    glAttachShader(voronoidesign_program,add_handle);
    glAttachShader(voronoidesign_program,smoothmin_handle);
    glAttachShader(voronoidesign_program,dspline3_handle);
    glAttachShader(voronoidesign_program,dvoronoi_handle);
    glAttachShader(voronoidesign_program,normal_handle);
    glLinkProgram(voronoidesign_program);
#ifdef DEBUG
    printf("---> voronoidesign Program:\n");
    debugp(voronoidesign_program);
    printf(">>>>\n");
#endif
    glUseProgram(voronoidesign_program);
    voronoidesign_iTime_location = glGetUniformLocation(voronoidesign_program, "iTime");
    voronoidesign_iResolution_location = glGetUniformLocation(voronoidesign_program, "iResolution");
    voronoidesign_iFader0_location = glGetUniformLocation(voronoidesign_program, "iFader0");
    voronoidesign_iFader1_location = glGetUniformLocation(voronoidesign_program, "iFader1");
    voronoidesign_iFader2_location = glGetUniformLocation(voronoidesign_program, "iFader2");
    voronoidesign_iFader3_location = glGetUniformLocation(voronoidesign_program, "iFader3");
    voronoidesign_iFader4_location = glGetUniformLocation(voronoidesign_program, "iFader4");
    voronoidesign_iFader5_location = glGetUniformLocation(voronoidesign_program, "iFader5");
    voronoidesign_iFader6_location = glGetUniformLocation(voronoidesign_program, "iFader6");
    voronoidesign_iFader7_location = glGetUniformLocation(voronoidesign_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadgroundboxes()
{
    int groundboxes_size = strlen(groundboxes_source);
    groundboxes_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(groundboxes_handle, 1, (GLchar **)&groundboxes_source, &groundboxes_size);
    glCompileShader(groundboxes_handle);
#ifdef DEBUG
    printf("---> groundboxes Shader:\n");
    debug(groundboxes_handle);
    printf(">>>>\n");
#endif
    groundboxes_program = glCreateProgram();
    glAttachShader(groundboxes_program,groundboxes_handle);
    glAttachShader(groundboxes_program,rand_handle);
    glAttachShader(groundboxes_program,lfnoise_handle);
    glAttachShader(groundboxes_program,dspline3_handle);
    glAttachShader(groundboxes_program,dbox3_handle);
    glAttachShader(groundboxes_program,dlinesegment3_handle);
    glAttachShader(groundboxes_program,stroke_handle);
    glAttachShader(groundboxes_program,zextrude_handle);
    glAttachShader(groundboxes_program,scale_handle);
    glAttachShader(groundboxes_program,smoothmin_handle);
    glAttachShader(groundboxes_program,add_handle);
    glAttachShader(groundboxes_program,dvoronoi_handle);
    glAttachShader(groundboxes_program,dbox_handle);
    glAttachShader(groundboxes_program,rot3_handle);
    glAttachShader(groundboxes_program,normal_handle);
    glLinkProgram(groundboxes_program);
#ifdef DEBUG
    printf("---> groundboxes Program:\n");
    debugp(groundboxes_program);
    printf(">>>>\n");
#endif
    glUseProgram(groundboxes_program);
    groundboxes_iTime_location = glGetUniformLocation(groundboxes_program, "iTime");
    groundboxes_iResolution_location = glGetUniformLocation(groundboxes_program, "iResolution");
    groundboxes_iFader0_location = glGetUniformLocation(groundboxes_program, "iFader0");
    groundboxes_iFader1_location = glGetUniformLocation(groundboxes_program, "iFader1");
    groundboxes_iFader2_location = glGetUniformLocation(groundboxes_program, "iFader2");
    groundboxes_iFader3_location = glGetUniformLocation(groundboxes_program, "iFader3");
    groundboxes_iFader4_location = glGetUniformLocation(groundboxes_program, "iFader4");
    groundboxes_iFader5_location = glGetUniformLocation(groundboxes_program, "iFader5");
    groundboxes_iFader6_location = glGetUniformLocation(groundboxes_program, "iFader6");
    groundboxes_iFader7_location = glGetUniformLocation(groundboxes_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadgraffiti()
{
    int graffiti_size = strlen(graffiti_source);
    graffiti_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(graffiti_handle, 1, (GLchar **)&graffiti_source, &graffiti_size);
    glCompileShader(graffiti_handle);
#ifdef DEBUG
    printf("---> graffiti Shader:\n");
    debug(graffiti_handle);
    printf(">>>>\n");
#endif
    graffiti_program = glCreateProgram();
    glAttachShader(graffiti_program,graffiti_handle);
    glAttachShader(graffiti_program,rand_handle);
    glAttachShader(graffiti_program,lfnoise_handle);
    glAttachShader(graffiti_program,mfnoise_handle);
    glAttachShader(graffiti_program,dtriangle_handle);
    glAttachShader(graffiti_program,dbox_handle);
    glAttachShader(graffiti_program,dlinesegment_handle);
    glAttachShader(graffiti_program,stroke_handle);
    glAttachShader(graffiti_program,dvoronoi_handle);
    glAttachShader(graffiti_program,rot3_handle);
    glAttachShader(graffiti_program,scale_handle);
    glAttachShader(graffiti_program,zextrude_handle);
    glAttachShader(graffiti_program,add_handle);
    glAttachShader(graffiti_program,normal_handle);
    glLinkProgram(graffiti_program);
#ifdef DEBUG
    printf("---> graffiti Program:\n");
    debugp(graffiti_program);
    printf(">>>>\n");
#endif
    glUseProgram(graffiti_program);
    graffiti_iTime_location = glGetUniformLocation(graffiti_program, "iTime");
    graffiti_iResolution_location = glGetUniformLocation(graffiti_program, "iResolution");
    graffiti_iFader0_location = glGetUniformLocation(graffiti_program, "iFader0");
    graffiti_iFader1_location = glGetUniformLocation(graffiti_program, "iFader1");
    graffiti_iFader2_location = glGetUniformLocation(graffiti_program, "iFader2");
    graffiti_iFader3_location = glGetUniformLocation(graffiti_program, "iFader3");
    graffiti_iFader4_location = glGetUniformLocation(graffiti_program, "iFader4");
    graffiti_iFader5_location = glGetUniformLocation(graffiti_program, "iFader5");
    graffiti_iFader6_location = glGetUniformLocation(graffiti_program, "iFader6");
    graffiti_iFader7_location = glGetUniformLocation(graffiti_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadgreet()
{
    int greet_size = strlen(greet_source);
    greet_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(greet_handle, 1, (GLchar **)&greet_source, &greet_size);
    glCompileShader(greet_handle);
#ifdef DEBUG
    printf("---> greet Shader:\n");
    debug(greet_handle);
    printf(">>>>\n");
#endif
    greet_program = glCreateProgram();
    glAttachShader(greet_program,greet_handle);
    glAttachShader(greet_program,dpolygon_handle);
    glAttachShader(greet_program,rot_handle);
    glAttachShader(greet_program,dcircle_handle);
    glAttachShader(greet_program,dbox_handle);
    glAttachShader(greet_program,dlinesegment_handle);
    glAttachShader(greet_program,dtriangle_handle);
    glAttachShader(greet_program,scale_handle);
    glAttachShader(greet_program,rand_handle);
    glAttachShader(greet_program,hash31_handle);
    glAttachShader(greet_program,lfnoise_handle);
    glAttachShader(greet_program,dbox_handle);
    glAttachShader(greet_program,stroke_handle);
    glAttachShader(greet_program,zextrude_handle);
    glAttachShader(greet_program,add_handle);
    glAttachShader(greet_program,smoothmin_handle);
    glAttachShader(greet_program,dbox3_handle);
    glAttachShader(greet_program,dschnappsgirls_handle);
    glAttachShader(greet_program,dspacepigs_handle);
    glAttachShader(greet_program,dkewlers_handle);
    glAttachShader(greet_program,dfarbrausch_handle);
    glAttachShader(greet_program,dhaujobb_handle);
    glAttachShader(greet_program,dmercury_handle);
    glAttachShader(greet_program,normal_handle);
    glLinkProgram(greet_program);
#ifdef DEBUG
    printf("---> greet Program:\n");
    debugp(greet_program);
    printf(">>>>\n");
#endif
    glUseProgram(greet_program);
    greet_iTime_location = glGetUniformLocation(greet_program, "iTime");
    greet_iResolution_location = glGetUniformLocation(greet_program, "iResolution");
    greet_iFader0_location = glGetUniformLocation(greet_program, "iFader0");
    greet_iFader1_location = glGetUniformLocation(greet_program, "iFader1");
    greet_iFader2_location = glGetUniformLocation(greet_program, "iFader2");
    greet_iFader3_location = glGetUniformLocation(greet_program, "iFader3");
    greet_iFader4_location = glGetUniformLocation(greet_program, "iFader4");
    greet_iFader5_location = glGetUniformLocation(greet_program, "iFader5");
    greet_iFader6_location = glGetUniformLocation(greet_program, "iFader6");
    greet_iFader7_location = glGetUniformLocation(greet_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadevoke()
{
    int evoke_size = strlen(evoke_source);
    evoke_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(evoke_handle, 1, (GLchar **)&evoke_source, &evoke_size);
    glCompileShader(evoke_handle);
#ifdef DEBUG
    printf("---> evoke Shader:\n");
    debug(evoke_handle);
    printf(">>>>\n");
#endif
    evoke_program = glCreateProgram();
    glAttachShader(evoke_program,evoke_handle);
    glAttachShader(evoke_program,rand_handle);
    glAttachShader(evoke_program,dbox_handle);
    glAttachShader(evoke_program,stroke_handle);
    glAttachShader(evoke_program,dlinesegment_handle);
    glAttachShader(evoke_program,lfnoise_handle);
    glAttachShader(evoke_program,dvoronoi_handle);
    glAttachShader(evoke_program,scale_handle);
    glLinkProgram(evoke_program);
#ifdef DEBUG
    printf("---> evoke Program:\n");
    debugp(evoke_program);
    printf(">>>>\n");
#endif
    glUseProgram(evoke_program);
    evoke_iTime_location = glGetUniformLocation(evoke_program, "iTime");
    evoke_iResolution_location = glGetUniformLocation(evoke_program, "iResolution");
    evoke_iFader0_location = glGetUniformLocation(evoke_program, "iFader0");
    evoke_iFader1_location = glGetUniformLocation(evoke_program, "iFader1");
    evoke_iFader2_location = glGetUniformLocation(evoke_program, "iFader2");
    evoke_iFader3_location = glGetUniformLocation(evoke_program, "iFader3");
    evoke_iFader4_location = glGetUniformLocation(evoke_program, "iFader4");
    evoke_iFader5_location = glGetUniformLocation(evoke_program, "iFader5");
    evoke_iFader6_location = glGetUniformLocation(evoke_program, "iFader6");
    evoke_iFader7_location = glGetUniformLocation(evoke_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadcanal()
{
    int canal_size = strlen(canal_source);
    canal_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(canal_handle, 1, (GLchar **)&canal_source, &canal_size);
    glCompileShader(canal_handle);
#ifdef DEBUG
    printf("---> canal Shader:\n");
    debug(canal_handle);
    printf(">>>>\n");
#endif
    canal_program = glCreateProgram();
    glAttachShader(canal_program,canal_handle);
    glAttachShader(canal_program,scale_handle);
    glAttachShader(canal_program,rand_handle);
    glAttachShader(canal_program,hash31_handle);
    glAttachShader(canal_program,lfnoise_handle);
    glAttachShader(canal_program,stroke_handle);
    glAttachShader(canal_program,zextrude_handle);
    glAttachShader(canal_program,add_handle);
    glAttachShader(canal_program,smoothmin_handle);
    glAttachShader(canal_program,dsmoothvoronoi_handle);
    glAttachShader(canal_program,normal_handle);
    glLinkProgram(canal_program);
#ifdef DEBUG
    printf("---> canal Program:\n");
    debugp(canal_program);
    printf(">>>>\n");
#endif
    glUseProgram(canal_program);
    canal_iTime_location = glGetUniformLocation(canal_program, "iTime");
    canal_iResolution_location = glGetUniformLocation(canal_program, "iResolution");
    canal_iFader0_location = glGetUniformLocation(canal_program, "iFader0");
    canal_iFader1_location = glGetUniformLocation(canal_program, "iFader1");
    canal_iFader2_location = glGetUniformLocation(canal_program, "iFader2");
    canal_iFader3_location = glGetUniformLocation(canal_program, "iFader3");
    canal_iFader4_location = glGetUniformLocation(canal_program, "iFader4");
    canal_iFader5_location = glGetUniformLocation(canal_program, "iFader5");
    canal_iFader6_location = glGetUniformLocation(canal_program, "iFader6");
    canal_iFader7_location = glGetUniformLocation(canal_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadtext()
{
    int text_size = strlen(text_source);
    text_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(text_handle, 1, (GLchar **)&text_source, &text_size);
    glCompileShader(text_handle);
#ifdef DEBUG
    printf("---> text Shader:\n");
    debug(text_handle);
    printf(">>>>\n");
#endif
    text_program = glCreateProgram();
    glAttachShader(text_program,text_handle);
    glAttachShader(text_program,rand_handle);
    glAttachShader(text_program,lfnoise_handle);
    glAttachShader(text_program,rshort_handle);
    glAttachShader(text_program,rfloat_handle);
    glAttachShader(text_program,dbox_handle);
    glAttachShader(text_program,dcircle_handle);
    glAttachShader(text_program,dlinesegment_handle);
    glAttachShader(text_program,drhomboid_handle);
    glAttachShader(text_program,dcirclesegment_handle);
    glAttachShader(text_program,stroke_handle);
    glAttachShader(text_program,dglyph_handle);
    glAttachShader(text_program,dstring_handle);
    glAttachShader(text_program,dfloat_handle);
    glAttachShader(text_program,smoothmin_handle);
    glAttachShader(text_program,dint_handle);
    glAttachShader(text_program,dtime_handle);
    glAttachShader(text_program,window_handle);
    glAttachShader(text_program,progressbar_handle);
    glAttachShader(text_program,dvoronoi_handle);
    glLinkProgram(text_program);
#ifdef DEBUG
    printf("---> text Program:\n");
    debugp(text_program);
    printf(">>>>\n");
#endif
    glUseProgram(text_program);
    text_iFontWidth_location = glGetUniformLocation(text_program, "iFontWidth");
    text_iTime_location = glGetUniformLocation(text_program, "iTime");
    text_iResolution_location = glGetUniformLocation(text_program, "iResolution");
    text_iChannel0_location = glGetUniformLocation(text_program, "iChannel0");
    text_iFont_location = glGetUniformLocation(text_program, "iFont");
    text_iFSAA_location = glGetUniformLocation(text_program, "iFSAA");
    progress += .2/(float)nprograms;
}

void Loadpost()
{
    int post_size = strlen(post_source);
    post_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(post_handle, 1, (GLchar **)&post_source, &post_size);
    glCompileShader(post_handle);
#ifdef DEBUG
    printf("---> post Shader:\n");
    debug(post_handle);
    printf(">>>>\n");
#endif
    post_program = glCreateProgram();
    glAttachShader(post_program,post_handle);
    glAttachShader(post_program,rand_handle);
    glAttachShader(post_program,lfnoise_handle);
    glAttachShader(post_program,stroke_handle);
    glAttachShader(post_program,dvoronoi_handle);
    glAttachShader(post_program,rot3_handle);
    glAttachShader(post_program,dbox3_handle);
    glAttachShader(post_program,add_handle);
    glAttachShader(post_program,normal_handle);
    glLinkProgram(post_program);
#ifdef DEBUG
    printf("---> post Program:\n");
    debugp(post_program);
    printf(">>>>\n");
#endif
    glUseProgram(post_program);
    post_iFSAA_location = glGetUniformLocation(post_program, "iFSAA");
    post_iResolution_location = glGetUniformLocation(post_program, "iResolution");
    post_iChannel0_location = glGetUniformLocation(post_program, "iChannel0");
    post_iTime_location = glGetUniformLocation(post_program, "iTime");
    progress += .2/(float)nprograms;
}

void Loadlogo210()
{
    int logo210_size = strlen(logo210_source);
    logo210_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(logo210_handle, 1, (GLchar **)&logo210_source, &logo210_size);
    glCompileShader(logo210_handle);
#ifdef DEBUG
    printf("---> logo210 Shader:\n");
    debug(logo210_handle);
    printf(">>>>\n");
#endif
    logo210_program = glCreateProgram();
    glAttachShader(logo210_program,logo210_handle);
    glAttachShader(logo210_program,dbox3_handle);
    glAttachShader(logo210_program,rot3_handle);
    glAttachShader(logo210_program,stroke_handle);
    glAttachShader(logo210_program,add_handle);
    glLinkProgram(logo210_program);
#ifdef DEBUG
    printf("---> logo210 Program:\n");
    debugp(logo210_program);
    printf(">>>>\n");
#endif
    glUseProgram(logo210_program);
    logo210_iTime_location = glGetUniformLocation(logo210_program, "iTime");
    logo210_iResolution_location = glGetUniformLocation(logo210_program, "iResolution");
    logo210_iFader0_location = glGetUniformLocation(logo210_program, "iFader0");
    logo210_iFader1_location = glGetUniformLocation(logo210_program, "iFader1");
    logo210_iFader2_location = glGetUniformLocation(logo210_program, "iFader2");
    logo210_iFader3_location = glGetUniformLocation(logo210_program, "iFader3");
    logo210_iFader4_location = glGetUniformLocation(logo210_program, "iFader4");
    logo210_iFader5_location = glGetUniformLocation(logo210_program, "iFader5");
    logo210_iFader6_location = glGetUniformLocation(logo210_program, "iFader6");
    logo210_iFader7_location = glGetUniformLocation(logo210_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadtransbubbles()
{
    int transbubbles_size = strlen(transbubbles_source);
    transbubbles_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(transbubbles_handle, 1, (GLchar **)&transbubbles_source, &transbubbles_size);
    glCompileShader(transbubbles_handle);
#ifdef DEBUG
    printf("---> transbubbles Shader:\n");
    debug(transbubbles_handle);
    printf(">>>>\n");
#endif
    transbubbles_program = glCreateProgram();
    glAttachShader(transbubbles_program,transbubbles_handle);
    glAttachShader(transbubbles_program,scale_handle);
    glAttachShader(transbubbles_program,rand_handle);
    glAttachShader(transbubbles_program,lfnoise_handle);
    glAttachShader(transbubbles_program,dbox3_handle);
    glAttachShader(transbubbles_program,rot3_handle);
    glAttachShader(transbubbles_program,stroke_handle);
    glAttachShader(transbubbles_program,hash13_handle);
    glAttachShader(transbubbles_program,add_handle);
    glAttachShader(transbubbles_program,smoothmin_handle);
    glAttachShader(transbubbles_program,normal_handle);
    glLinkProgram(transbubbles_program);
#ifdef DEBUG
    printf("---> transbubbles Program:\n");
    debugp(transbubbles_program);
    printf(">>>>\n");
#endif
    glUseProgram(transbubbles_program);
    transbubbles_iTime_location = glGetUniformLocation(transbubbles_program, "iTime");
    transbubbles_iResolution_location = glGetUniformLocation(transbubbles_program, "iResolution");
    transbubbles_iFader0_location = glGetUniformLocation(transbubbles_program, "iFader0");
    transbubbles_iFader1_location = glGetUniformLocation(transbubbles_program, "iFader1");
    transbubbles_iFader2_location = glGetUniformLocation(transbubbles_program, "iFader2");
    transbubbles_iFader3_location = glGetUniformLocation(transbubbles_program, "iFader3");
    transbubbles_iFader4_location = glGetUniformLocation(transbubbles_program, "iFader4");
    transbubbles_iFader5_location = glGetUniformLocation(transbubbles_program, "iFader5");
    transbubbles_iFader6_location = glGetUniformLocation(transbubbles_program, "iFader6");
    transbubbles_iFader7_location = glGetUniformLocation(transbubbles_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadvolclouds()
{
    int volclouds_size = strlen(volclouds_source);
    volclouds_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(volclouds_handle, 1, (GLchar **)&volclouds_source, &volclouds_size);
    glCompileShader(volclouds_handle);
#ifdef DEBUG
    printf("---> volclouds Shader:\n");
    debug(volclouds_handle);
    printf(">>>>\n");
#endif
    volclouds_program = glCreateProgram();
    glAttachShader(volclouds_program,volclouds_handle);
    glAttachShader(volclouds_program,scale_handle);
    glLinkProgram(volclouds_program);
#ifdef DEBUG
    printf("---> volclouds Program:\n");
    debugp(volclouds_program);
    printf(">>>>\n");
#endif
    glUseProgram(volclouds_program);
    volclouds_iTime_location = glGetUniformLocation(volclouds_program, "iTime");
    volclouds_iResolution_location = glGetUniformLocation(volclouds_program, "iResolution");
    volclouds_iFader0_location = glGetUniformLocation(volclouds_program, "iFader0");
    volclouds_iFader1_location = glGetUniformLocation(volclouds_program, "iFader1");
    volclouds_iFader2_location = glGetUniformLocation(volclouds_program, "iFader2");
    volclouds_iFader3_location = glGetUniformLocation(volclouds_program, "iFader3");
    volclouds_iFader4_location = glGetUniformLocation(volclouds_program, "iFader4");
    volclouds_iFader5_location = glGetUniformLocation(volclouds_program, "iFader5");
    volclouds_iFader6_location = glGetUniformLocation(volclouds_program, "iFader6");
    volclouds_iFader7_location = glGetUniformLocation(volclouds_program, "iFader7");
    progress += .2/(float)nprograms;
}

void Loadchart()
{
    int chart_size = strlen(chart_source);
    chart_handle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(chart_handle, 1, (GLchar **)&chart_source, &chart_size);
    glCompileShader(chart_handle);
#ifdef DEBUG
    printf("---> chart Shader:\n");
    debug(chart_handle);
    printf(">>>>\n");
#endif
    chart_program = glCreateProgram();
    glAttachShader(chart_program,chart_handle);
    glLinkProgram(chart_program);
#ifdef DEBUG
    printf("---> chart Program:\n");
    debugp(chart_program);
    printf(">>>>\n");
#endif
    glUseProgram(chart_program);
    chart_iTime_location = glGetUniformLocation(chart_program, "iTime");
    chart_iResolution_location = glGetUniformLocation(chart_program, "iResolution");
    chart_iFader0_location = glGetUniformLocation(chart_program, "iFader0");
    chart_iFader1_location = glGetUniformLocation(chart_program, "iFader1");
    chart_iFader2_location = glGetUniformLocation(chart_program, "iFader2");
    chart_iFader3_location = glGetUniformLocation(chart_program, "iFader3");
    chart_iFader4_location = glGetUniformLocation(chart_program, "iFader4");
    chart_iFader5_location = glGetUniformLocation(chart_program, "iFader5");
    chart_iFader6_location = glGetUniformLocation(chart_program, "iFader6");
    chart_iFader7_location = glGetUniformLocation(chart_program, "iFader7");
    progress += .2/(float)nprograms;
}

void LoadPrograms()
{
    Loadvoronoidesign();
    updateBar();
    Loadgroundboxes();
    updateBar();
    Loadgraffiti();
    updateBar();
    Loadgreet();
    updateBar();
    Loadevoke();
    updateBar();
    Loadcanal();
    updateBar();
    Loadtext();
    updateBar();
    Loadpost();
    updateBar();
    Loadlogo210();
    updateBar();
    Loadtransbubbles();
    updateBar();
    Loadvolclouds();
    updateBar();
    Loadchart();
    updateBar();
}
#endif
