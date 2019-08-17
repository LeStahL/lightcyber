//Generated with Symbolize (c) 2019 Alexander Kraus <nr4@z10.info>.
#ifndef SYMBOLIZE_H
#define SYMBOLIZE_H

extern float progress;int scale_handle, dsmoothvoronoi_handle, rand_handle, hash31_handle, lfnoise_handle, mfnoise_handle, dbox_handle, dlinesegment3_handle, stroke_handle, zextrude_handle, add_handle, smoothmin_handle, dspline3_handle, dvoronoi_handle, normal_handle, dbox3_handle, rot3_handle, dtriangle_handle, dlinesegment_handle, dpolygon_handle, rot_handle, dcircle_handle, dschnappsgirls_handle, dspacepigs_handle, dkewlers_handle, dfarbrausch_handle, dhaujobb_handle, dmercury_handle, rshort_handle, rfloat_handle, drhomboid_handle, dcirclesegment_handle, dglyph_handle, dstring_handle, dfloat_handle, dint_handle, dtime_handle, window_handle, progressbar_handle, hash13_handle;
const int nsymbols = 40;
const char *scale_source = "#version 130\nuniform float iTime;void scale(out float s){if(iTime>=0.0&&iTime<4.705882){s=mod(iTime+.3-0.0,0.588225)-0.2941125;s=smoothstep(-0.04901875,0.,s)*(1.-smoothstep(0.,0.14705625,s));}if(iTime>=4.705882&&iTime<18.552036){s=mod(iTime+.3-4.705882,0.576975)-0.2884875;s=smoothstep(-0.04808125,0.,s)*(1.-smoothstep(0.,0.14424375,s));}if(iTime>=18.552036&&iTime<22.996481){s=mod(iTime+.3-18.552036,0.55555)-0.277775;s=smoothstep(-0.046295833333333335,0.,s)*(1.-smoothstep(0.,0.1388875,s));}if(iTime>=22.996481&&iTime<25.139338){s=mod(iTime+.3-22.996481,0.535675)-0.2678375;s=smoothstep(-0.04463958333333334,0.,s)*(1.-smoothstep(0.,0.13391875,s));}if(iTime>=25.139338&&iTime<27.208303){s=mod(iTime+.3-25.139338,0.517275)-0.2586375;s=smoothstep(-0.043106250000000006,0.,s)*(1.-smoothstep(0.,0.12931875,s));}if(iTime>=27.208303&&iTime<65.208303){s=mod(iTime+.3-27.208303,0.5)-0.25;s=smoothstep(-0.041666666666666664,0.,s)*(1.-smoothstep(0.,0.125,s));}if(iTime>=65.208303&&iTime<71.109943){s=mod(iTime+.3-65.208303,0.491825)-0.2459125;s=smoothstep(-0.04098541666666667,0.,s)*(1.-smoothstep(0.,0.12295625,s));}if(iTime>=71.109943&&iTime<82.722846){s=mod(iTime+.3-71.109943,0.48385)-0.241925;s=smoothstep(-0.040320833333333334,0.,s)*(1.-smoothstep(0.,0.1209625,s));}if(iTime>=82.722846&&iTime<86.722846){s=mod(iTime+.3-82.722846,0.5)-0.25;s=smoothstep(-0.041666666666666664,0.,s)*(1.-smoothstep(0.,0.125,s));}if(iTime>=86.722846&&iTime<94.998708){s=mod(iTime+.3-86.722846,0.517275)-0.2586375;s=smoothstep(-0.043106250000000006,0.,s)*(1.-smoothstep(0.,0.12931875,s));}if(iTime>=94.998708&&iTime<103.134301){s=mod(iTime+.3-94.998708,0.50845)-0.254225;s=smoothstep(-0.04237083333333333,0.,s)*(1.-smoothstep(0.,0.1271125,s));}if(iTime>=103.134301&&iTime<105.134301){s=mod(iTime+.3-103.134301,0.5)-0.25;s=smoothstep(-0.041666666666666664,0.,s)*(1.-smoothstep(0.,0.125,s));}if(iTime>=105.134301&&iTime<107.069785){s=mod(iTime+.3-105.134301,0.48385)-0.241925;s=smoothstep(-0.040320833333333334,0.,s)*(1.-smoothstep(0.,0.1209625,s));}if(iTime>=107.069785&&iTime<129.569785){s=mod(iTime+.3-107.069785,0.468775)-0.2343875;s=smoothstep(-0.03906458333333333,0.,s)*(1.-smoothstep(0.,0.11719375,s));}if(iTime>=129.569785&&iTime<131.505269){s=mod(iTime+.3-129.569785,0.48385)-0.241925;s=smoothstep(-0.040320833333333334,0.,s)*(1.-smoothstep(0.,0.1209625,s));}if(iTime>=131.505269&&iTime<141.341334){s=mod(iTime+.3-131.505269,0.491825)-0.2459125;s=smoothstep(-0.04098541666666667,0.,s)*(1.-smoothstep(0.,0.12295625,s));}if(iTime>=141.341334&&iTime<143.276818){s=mod(iTime+.3-141.341334,0.48385)-0.241925;s=smoothstep(-0.040320833333333334,0.,s)*(1.-smoothstep(0.,0.1209625,s));}if(iTime>=143.276818&&iTime<145.151818){s=mod(iTime+.3-143.276818,0.468775)-0.2343875;s=smoothstep(-0.03906458333333333,0.,s)*(1.-smoothstep(0.,0.11719375,s));}if(iTime>=145.151818&&iTime<186.97){s=mod(iTime+.3-145.151818,0.45455)-0.227275;s=smoothstep(-0.037879166666666665,0.,s)*(1.-smoothstep(0.,0.1136375,s));}}\0";
const char *dsmoothvoronoi_source = "#version 130\nuniform float iTime;uniform float iFader0;const vec3 c=vec3 (1.,0.,-1.);void lfnoise(in vec2 t,out float n);void smoothmin(in float a,in float b,in float k,out float dst);void rand(in vec2 x,out float d);void dsmoothvoronoi(in vec2 x,out float d,out vec2 z){float n;vec2 y=floor(x);float ret=1.;vec2 pf=c.yy,p;float df=10.;for(int i=-1;i<=1;i+=1)for(int j=-1;j<=1;j+=1){p=y+vec2 (float (i),float (j));float pa;rand(p,pa);p+=pa;d=length(x-p);if(d<df){df=d;pf=p;}}for(int i=-1;i<=1;i+=1)for(int j=-1;j<=1;j+=1){p=y+vec2 (float (i),float (j));float pa;rand(p,pa);p+=pa;vec2 o=p-pf;d=length(.5*o-dot(x-pf,o)/dot(o,o)*o);smoothmin(ret,d,.05,ret);}d=ret;z=pf;}\0";
const char *rand_source = "#version 130\nvoid rand(in vec2 x,out float n){x+=400.;n=fract(sin(dot(sign(x)*abs(x),vec2 (12.9898,78.233)))*43758.5453);}\0";
const char *hash31_source = "void hash31(in float p,out vec3 d){vec3 p3=fract(vec3 (p)*vec3 (.1031,.1030,.0973));p3+=dot(p3,p3.yzx+33.33);d=fract((p3.xxy+p3.yzz)*p3.zyx);}\0";
const char *lfnoise_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void rand(in vec2 x,out float d);void lfnoise(in vec2 t,out float n){vec2 i=floor(t);t=fract(t);t=smoothstep(c.yy,c.xx,t);vec2 v1,v2;rand(i,v1.x);rand(i+c.xy,v1.y);rand(i+c.yx,v2.x);rand(i+c.xx,v2.y);v1=c.zz+2.*mix(v1,v2,t.y);n=mix(v1.x,v1.y,t.x);}\0";
const char *mfnoise_source = "#version 130\nvoid lfnoise(in vec2 x,out float d);void mfnoise(in vec2 x,in float d,in float b,in float e,out float n){n=0.;float a=1.,nf=0.,buf;for(float f=d;f<b;f*=2.){lfnoise(f*x,buf);n+=a*buf;a*=e;nf+=1.;}n*=(1.-e)/(1.-pow(e,nf));}\0";
const char *dbox_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void dbox(in vec2 x,in vec2 b,out float d){vec2 da=abs(x)-b;d=length(max(da,c.yy))+min(max(da.x,da.y),0.0);}\0";
const char *dlinesegment3_source = "#version 130\nvoid dlinesegment3(in vec3 x,in vec3 p1,in vec3 p2,out float d){vec3 da=p2-p1;d=length(x-mix(p1,p2,clamp(dot(x-p1,da)/dot(da,da),0.,1.)));}\0";
const char *stroke_source = "void stroke(in float d0,in float s,out float d){d=abs(d0)-s;}\0";
const char *zextrude_source = "void zextrude(in float z,in float d2d,in float h,out float d){vec2 w=vec2 (-d2d,abs(z)-0.5*h);d=length(max(w,0.0));}\0";
const char *add_source = "void add(in vec2 sda,in vec2 sdb,out vec2 sdf){sdf=mix(sda,sdb,step(sdb.x,sda.x));}\0";
const char *smoothmin_source = "void smoothmin(in float a,in float b,in float k,out float dst){float h=max(k-abs(a-b),0.0)/k;dst=min(a,b)-h*h*h*k*(1.0/6.0);}\0";
const char *dspline3_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);const float pi=acos(-1.);float dist3(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t){t=clamp(t,0.,1.);return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);}void dspline3(in vec3 x,in vec3 p0,in vec3 p1,in vec3 p2,out float ds){vec3 E=x-p0,F=p2-2.*p1+p0,G=p1-p0,ai=vec3 (3.*dot(G,F),2.*dot(G,G)-dot(E,F),-dot(E,G))/dot(F,F);float tau=ai.x/3.,p=ai.y-tau*ai.x,q=-tau*(tau*tau+p)+ai.z,dis=q*q/4.+p*p*p/27.;if(dis>0.){vec2 ki=-.5*q*c.xx+sqrt(dis)*c.xz,ui=sign(ki)*pow(abs(ki),c.xx/3.);ds=dist3(p0,p1,p2,x,ui.x+ui.y-tau);return ;}float fac=sqrt(-4./3.*p),arg=acos(-.5*q*sqrt(-27./p/p/p))/3.;vec3 t=c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;ds=min(dist3(p0,p1,p2,x,t.x),min(dist3(p0,p1,p2,x,t.y),dist3(p0,p1,p2,x,t.z)));}\0";
const char *dvoronoi_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void rand(in vec2 x,out float d);void dvoronoi(in vec2 x,out float d,out vec2 z){vec2 y=floor(x);float ret=1.;vec2 pf=c.yy,p;float df=10.;for(int i=-1;i<=1;i+=1)for(int j=-1;j<=1;j+=1){p=y+vec2 (float (i),float (j));float pa;rand(p,pa);p+=pa;d=length(x-p);if(d<df){df=d;pf=p;}}for(int i=-1;i<=1;i+=1)for(int j=-1;j<=1;j+=1){p=y+vec2 (float (i),float (j));float pa;rand(p,pa);p+=pa;vec2 o=p-pf;d=length(.5*o-dot(x-pf,o)/dot(o,o)*o);ret=min(ret,d);}d=ret;z=pf;}\0";
const char *normal_source = "const vec3 c=vec3 (1.0,0.0,-1.0);void scene(in vec3 x,out vec2 s);void normal(in vec3 x,out vec3 n,in float dx){vec2 s,na;scene(x,s);scene(x+dx*c.xyy,na);n.x=na.x;scene(x+dx*c.yxy,na);n.y=na.x;scene(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);}\0";
const char *dbox3_source = "#version 130\nvoid dbox3(in vec3 x,in vec3 b,out float d){vec3 da=abs(x)-b;d=length(max(da,0.0))+min(max(da.x,max(da.y,da.z)),0.0);}\0";
const char *rot3_source = "const vec3 c=vec3 (1.,0.,-1.);void rot3(in vec3 p,out mat3 rot){rot=mat3 (c.xyyy,cos(p.x),sin(p.x),0.,-sin(p.x),cos(p.x))*mat3 (cos(p.y),0.,-sin(p.y),c.yxy,sin(p.y),0.,cos(p.y))*mat3 (cos(p.z),-sin(p.z),0.,sin(p.z),cos(p.z),c.yyyx);}\0";
const char *dtriangle_source = "void dtriangle(in vec2 p,in vec2 p0,in vec2 p1,in vec2 p2,out float dst){vec2 e0=p1-p0;vec2 e1=p2-p1;vec2 e2=p0-p2;vec2 v0=p-p0;vec2 v1=p-p1;vec2 v2=p-p2;vec2 pq0=v0-e0*clamp(dot(v0,e0)/dot(e0,e0),0.0,1.0);vec2 pq1=v1-e1*clamp(dot(v1,e1)/dot(e1,e1),0.0,1.0);vec2 pq2=v2-e2*clamp(dot(v2,e2)/dot(e2,e2),0.0,1.0);float s=sign(e0.x*e2.y-e0.y*e2.x);vec2 d=min(min(vec2 (dot(pq0,pq0),s*(v0.x*e0.y-v0.y*e0.x)),vec2 (dot(pq1,pq1),s*(v1.x*e1.y-v1.y*e1.x))),vec2 (dot(pq2,pq2),s*(v2.x*e2.y-v2.y*e2.x)));dst=-sqrt(d.x)*sign(d.y);}\0";
const char *dlinesegment_source = "#version 130\nvoid dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d){vec2 da=p2-p1;d=length(x-mix(p1,p2,clamp(dot(x-p1,da)/dot(da,da),0.,1.)));}\0";
const char *dpolygon_source = "#version 130\nconst float pi=acos(-1.);void dpolygon(in vec2 x,in float N,out float d){d=2.0*pi/N;float t=mod(acos(x.x/length(x)),d)-0.5*d;d=-0.5+length(x)*cos(t)/cos(0.5*d);}\0";
const char *rot_source = "#version 130\nvoid rot(in float phi,out mat2 m){vec2 cs=vec2 (cos(phi),sin(phi));m=mat2 (cs.x,-cs.y,cs.y,cs.x);}\0";
const char *dcircle_source = "#version 130\nvoid dcircle(in vec2 x,out float d){d=length(x)-1.0;}\0";
const char *dschnappsgirls_source = "#version 130\nvoid dtriangle(in vec2 x,in vec2 p0,in vec2 p1,in vec2 p2,out float d);void dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d);void stroke(in float d0,in float s,out float d);void dcircle(in vec2 x,out float d);void dpolygon(in vec2 x,in float N,out float d);void dschnappsgirls(in vec2 x,out float d){dpolygon(.5*x,6.0,d);float da,d0;dtriangle(x,vec2 (-.1,-.3),vec2 (.5,-.3),vec2 (.2,.6),d0);dlinesegment(x,vec2 (-.1,.325),vec2 (.5,.325),da);stroke(da,.06,da);d0=max(d0,-da);dcircle(7.*(x-vec2 (.2,.5)),da);d0=max(d0,-da+.5);d0=min(d0,da/7.);dlinesegment(x,vec2 (.125,-.3),vec2 (.125,-.6),da);stroke(da,.06,da);d0=min(d0,da);dlinesegment(x,vec2 (.275,-.3),vec2 (.275,-.6),da);stroke(da,.06,da);d0=min(d0,da);dlinesegment(x,vec2 (0.05,.25),vec2 (.35,.25),da);stroke(da,.085,da);d0=min(d0,da);dlinesegment(x,vec2 (.385,.25),vec2 (.5,-.1),da);stroke(da,.055,da);d0=min(d0,da);dlinesegment(x,vec2 (.017,.25),vec2 (-.1,-.1),da);stroke(da,.055,da);d0=min(d0,da);dtriangle(x,vec2 (-.6,.3),vec2 (-.4,.1),vec2 (-.2,.3),da);stroke(da,.0125,da);d0=min(d0,da);dlinesegment(x,vec2 (-.4,.15),vec2 (-.4,-.1),da);stroke(da,.0125,da);d0=min(d0,da);dtriangle(x,vec2 (-.5,-.15),vec2 (-.3,-.15),vec2 (-.4,-.1),da);d0=min(d0,da);dtriangle(x,vec2 (-.55,.25),vec2 (-.4,.1),vec2 (-.25,.25),da);d0=min(d0,da);dlinesegment(x,vec2 (-.4,.1),vec2 (-.2,.5),da);stroke(da,.01,da);d0=min(d0,da);dcircle(24.*(x-vec2 (-.3,.3)),da);d0=min(d0,da/24.);dcircle(24.*(x-vec2 (-.25,.4)),da);d0=min(d0,da/24.);d=max(d,-d0);}\0";
const char *dspacepigs_source = "#version 130\nvoid dpolygon(in vec2 x,in float N,out float d);void dcircle(in vec2 x,out float d);void dear(in vec2 x,out float d){d=abs(2.*x.y)-.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))-.5*min(-abs(x.x),.01);}void dspacepigs(in vec2 x,out float d){dpolygon(.5*x,6.0,d);float da,d0;dcircle(2.5*x,d0);d0/=2.5;dear(vec2 (2.,5.)*x-vec2 (.8,1.3),da);d0=min(d0,da/10.);dear(vec2 (2.,5.)*x+vec2 (.8,-1.3),da);d0=min(d0,da/10.);dcircle(6.*x-vec2 (0.,-.5),da);d0=max(d0,-da/6.);dcircle(24.*x-vec2 (-1.5,-2.),da);d0=min(d0,da/24.);dcircle(24.*x-vec2 (1.5,-2.),da);d0=min(d0,da/24.);dcircle(16.*x-vec2 (-3.5,2.5),da);d0=max(d0,-da/16.);dcircle(16.*x-vec2 (3.5,2.5),da);d0=max(d0,-da/16.);dcircle(24.*x-vec2 (-5.,3.5),da);d0=min(d0,da/24.);dcircle(24.*x-vec2 (5.,3.5),da);d0=min(d0,da/24.);d=max(d,-d0);}\0";
const char *dkewlers_source = "#version 130\nvoid dbox(in vec2 x,in vec2 b,out float d);void dpolygon(in vec2 x,in float N,out float d);void dkewlers(in vec2 x,out float d){dpolygon(.5*x,6.0,d);float da,d0;x*=1.2;dbox(x-vec2 (0.,-.3),vec2 (.6,.1),d0);dbox(x-vec2 (-.5,-.0),vec2 (.1,.25),da);d0=min(d0,da);dbox(x-vec2 (-.5+1./3.,.25),vec2 (.1,.5),da);d0=min(d0,da);dbox(x-vec2 (-.5+2./3.,-.0),vec2 (.1,.25),da);d0=min(d0,da);dbox(x-vec2 (.5,-.0),vec2 (.1,.25),da);d0=min(d0,da);d=max(d,-d0);}\0";
const char *dfarbrausch_source = "#version 130\nvoid dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d);void dpolygon(in vec2 x,in float N,out float d);void stroke(in float d0,in float s,out float d);void dfarbrausch(in vec2 x,out float d){dpolygon(.5*x,6.0,d);float da,d0;x+=vec2 (.1,0.);x*=1.2;dlinesegment(x,vec2 (-.65,.05),vec2 (-.5,.05),d0);dlinesegment(x,vec2 (-.5,.05),vec2 (-.2,-.49),da);d0=min(d0,da);dlinesegment(x,vec2 (-.2,-.49),vec2 (-.0,-.49),da);d0=min(d0,da);dlinesegment(x,vec2 (-.0,-.49),vec2 (-.27,.0),da);d0=min(d0,da);dlinesegment(x,vec2 (-.07,0.),vec2 (-.27,.0),da);d0=min(d0,da);dlinesegment(x,vec2 (.2,-.49),vec2 (-.07,.0),da);d0=min(d0,da);dlinesegment(x,vec2 (.4,-.49),vec2 (.13,.0),da);d0=min(d0,da);dlinesegment(x,vec2 (.4,-.49),vec2 (.2,-.49),da);d0=min(d0,da);dlinesegment(x,vec2 (.33,0.),vec2 (.13,.0),da);d0=min(d0,da);dlinesegment(x,vec2 (.33,0.),vec2 (.51,-.33),da);d0=min(d0,da);dlinesegment(x,vec2 (.6,-.15),vec2 (.51,-.33),da);d0=min(d0,da);dlinesegment(x,vec2 (.53,0.),vec2 (.6,-.15),da);d0=min(d0,da);dlinesegment(x,vec2 (.7,0.),vec2 (.53,.0),da);d0=min(d0,da);dlinesegment(x,vec2 (.7,0.),vec2 (.68,-.04),da);d0=min(d0,da);dpolygon(5.*(x+vec2 (.3,.65)),6.,da);d0=min(d0,da/5.);dpolygon(5.*(x+vec2 (-.5,.65)),6.,da);d0=min(d0,da/5.);stroke(d0,.035,d0);d=max(d,-d0);}\0";
const char *dhaujobb_source = "#version 130\nconst vec3 c=vec3 (1.0,0.0,-1.0);const float pi=acos(-1.);void dpolygon(in vec2 x,in float N,out float d);void rot(in float phi,out mat2 m);void dcircle(in vec2 x,out float d);void dbox(in vec2 x,in vec2 b,out float d);void dhaujobb(in vec2 x,out float d){dpolygon(.5*x,6.0,d);float da,d0;mat2 m;rot(.3,m);x=1.1*x*m;x.x*=1.1;x+=vec2 (-.05,.2);dbox(x+.35*c.xx,vec2 (.1,.05),d0);dbox(x+vec2 (.3,.25),vec2 (.05,.15),da);d0=min(d0,da);dbox(x+vec2 (.2,.15),vec2 (.1,.05),da);d0=min(d0,da);dbox(x+vec2 (.15,.05),vec2 (.05,.15),da);d0=min(d0,da);dbox(x-vec2 (.65,.35),vec2 (.05,.15),da);d0=min(d0,da);rot(.2,m);dbox(m*(x-vec2 (.25,.15)),vec2 (.45,.05),da);d0=min(d0,da);dbox(m*(x-vec2 (-.15,.35)),vec2 (.45,.05),da);d0=min(d0,da);rot(pi/8.,m);dbox(m*(x-vec2 (.0,.25)),vec2 (.1,.15),da);d0=min(d0,da);dbox(m*(x-vec2 (.1,-.0)),vec2 (.025,.1),da);d0=min(d0,da);rot(.3,m);dbox(m*(x-vec2 (.235,.535)),vec2 (.035,.15),da);d0=min(d0,da);dbox(m*(x-vec2 (.225,.7)),vec2 (.075,.025),da);d0=min(d0,da);rot(-.2,m);dbox(m*(x+vec2 (.585,-.2)),vec2 (.0375,.1),da);d0=min(d0,da);dcircle(6.*(x-vec2 (-.15,.58)),da);d0=min(d0,da/6.);d0-=.05*(abs(x.x)+abs(x.y)-.2);d=max(d,-d0);}\0";
const char *dmercury_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void dbox(in vec2 x,in vec2 b,out float d);void dpolygon(in vec2 x,in float N,out float d);void dmercury(in vec2 x,out float d){dpolygon(.5*x,6.0,d);float da;x+=.1*c.yx;dbox(x-.35*c.yx,vec2 (.4,.35),da);d=max(d,-da);dbox(x-.7*c.yx,vec2 (.2,.2),da);d=min(d,da);dbox(x-.25*c.yx,vec2 (.2,.05),da);d=min(d,da);dbox(x+.2*c.yx,vec2 (.1,.4),da);d=max(d,-da);dbox(x+.2*c.yx,vec2 (.4,.1),da);d=max(d,-da);}\0";
const char *rshort_source = "#version 130\nuniform float iFontWidth;uniform sampler2D iFont;void rshort(in float off,out float val){float hilo=mod(off,2.);off*=.5;vec2 ind=(vec2 (mod(off,iFontWidth),floor(off/iFontWidth))+.05)/iFontWidth;vec4 block=texture(iFont,ind);vec2 data=mix(block.rg,block.ba,hilo);val=round(dot(vec2 (255.,65280.),data));}\0";
const char *rfloat_source = "#version 130\nvoid rshort(in float off,out float val);void rfloat(in float off,out float val){float d;rshort(off,d);float sign=floor(d/32768.),exponent=floor(d/1024.-sign*32.),significand=d-sign*32768.-exponent*1024.;if(exponent==0.){val=mix(1.,-1.,sign)*5.960464477539063e-08*significand;}else {val=mix(1.,-1.,sign)*(1.+significand*9.765625e-4)*pow(2.,exponent-15.);}}\0";
const char *drhomboid_source = "#version 130\nvoid dbox(in vec2 x,in vec2 b,out float dst);void drhomboid(in vec2 x,in vec2 b,in float tilt,out float dst){x.x-=tilt/2./b.y*x.y;dbox(x,b,dst);}\0";
const char *dcirclesegment_source = "#version 130\nconst float pi=acos(-1.);void dcirclesegment(in vec2 x,in float R,in float p0,in float p1,out float d){float p=atan(x.y,x.x);vec2 philo=vec2 (max(p0,p1),min(p0,p1));if((p<philo.x&&p>philo.y)||(p+2.0*pi<philo.x&&p+2.0*pi>philo.y)||(p-2.0*pi<philo.x&&p-2.0*pi>philo.y))d=abs(length(x)-R);else d=min(length(x-vec2 (cos(p0),sin(p0))),length(x-vec2 (cos(p1),sin(p1))));}\0";
const char *dglyph_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void rfloat(in float off,out float val);void dbox(in vec2 x,in vec2 b,out float dst);void dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d);void dcircle(in vec2 x,out float d);void dcirclesegment(in vec2 x,in float r,in float p0,in float p1,out float d);void stroke(in float d0,in float s,out float d);void smoothmin(in float a,in float b,in float k,out float dst);void dglyph(in vec2 x,in float ordinal,in float size,out float dst){float dis;dbox(x,2.*size*c.xx,dis);if(dis>0.){dst=dis+.5*size;return ;}float nglyphs,offset=0;rfloat(1.,nglyphs);for(float i=0.;i<nglyphs;i+=1.){float ord;rfloat(2.+2.*i,ord);ord=floor(ord);if(ord==ordinal){rfloat(2.+2.*i+1.,offset);offset=floor(offset);break;}}if(offset==0.){dst=1.;return ;}float d=1.,da=1.;float nlines;rfloat(offset,nlines);nlines=floor(nlines);offset+=1.;for(float i=0.;i<nlines;i+=1.){float x1;rfloat(offset,x1);offset+=1.;float y1;rfloat(offset,y1);offset+=1.;float x2;rfloat(offset,x2);offset+=1.;float y2;rfloat(offset,y2);offset+=1.;dlinesegment(x,size*vec2 (x1,y1),size*vec2 (x2,y2),da);d=min(d,da);}stroke(d,.2*size,d);float nsmoothlines,db=1.;da=1.;rfloat(offset,nsmoothlines);nsmoothlines=floor(nsmoothlines);offset+=1.;for(float i=0.;i<nsmoothlines;i+=1.){float x1;rfloat(offset,x1);offset+=1.;float y1;rfloat(offset,y1);offset+=1.;float x2;rfloat(offset,x2);offset+=1.;float y2;rfloat(offset,y2);offset+=1.;dlinesegment(x,size*vec2 (x1,y1),size*vec2 (x2,y2),db);da=min(da,db);}stroke(da,.2*size,da);smoothmin(d,da,.1*size,d);dst=d;}\0";
const char *dstring_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void rfloat(in float off,out float val);void dbox(in vec2 x,in vec2 b,out float dst);void dglyph(in vec2 x,in float ordinal,in float size,out float dst);void dstring(in vec2 x,in float ordinal,in float size,out float dst){float stroff0;rfloat(0.,stroff0);stroff0=floor(stroff0);float nstrings;rfloat(stroff0,nstrings);nstrings=floor(nstrings);if(ordinal>=nstrings){dst=1.;return ;}float stroff;rfloat(stroff0+1.+2.*ordinal,stroff);stroff=floor(stroff);float len;rfloat(stroff0+2.+2.*ordinal,len);len=floor(len);vec2 dx=mod(x-size,2.*size)-size,ind=ceil((x-dx+size)/2./size);float bound;dbox(x-size*(len-3.)*c.xy,vec2 (size*len,1.*size),bound);if(bound>0.){dst=bound+.5*size;return ;}float da;rfloat(stroff+ind.x,da);da=floor(da);dglyph(dx,da,.7*size,dst);}\0";
const char *dfloat_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void dglyph(in vec2 x,in float ordinal,in float size,out float dst);void dfloat(in vec2 x,in float num,in float size,out float dst){float d=1.,index=0.;float sign=sign(num),exp=0.;if(sign<0.){float da;dglyph(x,45.,.7*size,da);d=min(d,da);index+=1.;num*=-1.;}for(exp=-15.;exp<15.;exp+=1.)if(floor(num*pow(10.,exp))!=0.)break;exp*=-1.;for(float i=exp;i>=max(exp-5.,-33);i-=1.){float po=pow(10.,i);float ca=floor(num/po);num-=ca*po;float da;dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,da);d=min(d,da);index+=1.;if(i==exp){dglyph(x-2.*index*size*c.xy,46.,.7*size,da);d=min(d,da);index+=1.;}}float db;dglyph(x+.7*size*c.xy-2.*index*size*c.xy,101.,.7*size,db);d=min(d,db);index+=1.;if(exp<0.){dglyph(x+.7*size*c.xy-2.*index*size*c.xy,45.,.7*size,db);d=min(d,db);index+=1.;exp*=-1.;}float ca=floor(exp/10.);dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,db);d=min(d,db);index+=1.;ca=floor(exp-10.*ca);dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,db);d=min(d,db);index+=1.;dst=d;}\0";
const char *dint_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void dglyph(in vec2 x,in float ordinal,in float size,out float dst);void dint(in vec2 x,in float num,in float size,in float ndigits,out float dst){float d=1.,index=0.;if(num==0.){index=ndigits;dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.,.7*size,dst);return ;}float sign=sign(num),exp=0.;if(sign<0.){float da;dglyph(x,45.,.7*size,da);d=min(d,da);index+=1.;num*=-1.;}for(exp=-15.;exp<15.;exp+=1.)if(floor(num*pow(10.,exp))!=0.)break;exp*=-1.;int hit=0;for(float i=ndigits;i>=0.;i-=1.){float po=pow(10.,i);float ca=floor(num/po);if(ca==0.){if(hit==0){index+=1.;continue;}}else hit=1;num-=ca*po;float da;dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,da);d=min(d,da);index+=1.;}dst=d;}\0";
const char *dtime_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void dglyph(in vec2 x,in float ordinal,in float size,out float dst);void dtime(in vec2 x,in float num,in float size,out float dst){float d=1.,index=0.;num=floor(num);float ca=floor(num/600.),da=1.;dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,da);d=min(d,da);index+=1.;num-=ca*600.;ca=floor(num/60.);dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,da);d=min(d,da);index+=1.;num-=ca*60.;dglyph(x+.7*size*c.xy-2.*index*size*c.xy,58.,.7*size,da);d=min(d,da);index+=1.;ca=floor(num/10.);dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,da);d=min(d,da);index+=1.;num-=ca*10.;ca=floor(num/1.);dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size,da);d=min(d,da);dst=d;}\0";
const char *window_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);uniform float iTime;void dhexagonpattern(in vec2 p,out float d,out vec2 ind);void stroke(in float d0,in float s,out float d);void lfnoise(in vec2 t,out float num);void box(in vec2 x,in vec2 b,out float dst);void drhomboid(in vec2 x,in vec2 b,in float tilt,out float dst);void window(in vec2 x,in vec2 size,in vec3 bg,in float title_index,out vec4 col){}\0";
const char *progressbar_source = "#version 130\nconst vec3 c=vec3 (1.,0.,-1.);void progressbar(in vec2 x,in float width,in float progress,out vec4 col){}\0";
const char *hash13_source = "void hash13(in vec3 p3,out float d){p3=fract(p3*.1031);p3+=dot(p3,p3.yzx+33.33);d=fract((p3.x+p3.y)*p3.z);}\0";
const char *voronoidesign_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;float nbeats;float iScale;const vec3 c=vec3 (1.0,0.0,-1.0);const float pi=acos(-1.);void scale(out float s);void dsmoothvoronoi(in vec2 x,out float d,out vec2 z);void rand(in vec2 x,out float n);void hash31(in float p,out vec3 d);void lfnoise(in vec2 t,out float n);void mfnoise(in vec2 x,in float d,in float b,in float e,out float n);void dbox(in vec2 x,in vec2 b,out float d);void dlinesegment3(in vec3 x,in vec3 p1,in vec3 p2,out float d);void stroke(in float d0,in float s,out float d);void zextrude(in float z,in float d2d,in float h,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);void smoothmin(in float a,in float b,in float k,out float dst);void dspline3(in vec3 x,in vec3 p0,in vec3 p1,in vec3 p2,out float ds);void dvoronoi(in vec2 x,out float d,out vec2 z);vec2 vind,vind2;float v,fn,r1,fb;void scene(in vec3 x,out vec2 sdf){x.y=mix(x.y,-x.y,step(150.,iTime));x.y+=mix(.2,-.2,step(150.,iTime))*iTime;dvoronoi(1.5*x.xy,v,vind);vec3 y=vec3 (vind/1.5-x.xy,x.z);float n,n2;lfnoise(c.xx-.3*iTime+vind*3.,n);lfnoise(5.*x.z*c.xx-iTime-vind*4.,n2);n2*=.2;mat2 RR=mat2 (cos(n2),sin(n2),-sin(n2),cos(n2));vec2 a=x.xy;x.xy=RR*x.xy;rand(vind,r1);float r2;rand(vind-1336.,r2);float phi=atan(y.y,y.x),dp=pi/24.,phii=mod(phi,dp)-.5*dp,pa=phi-phii,R1=.05,R2=mix(.4,.25,1.-r2);R2=mix(R1,R2,.5+.5*n2);float r0;rand(pa*c.xx,r0);r0=mix(r0,.5+.5*n,.5);dspline3(y,vec3 (1.4*R1*cos(pa),1.4*R1*sin(pa),-.5),vec3 (R1*cos(pa),R1*sin(pa),.1*r1),vec3 (mix(R1,R2,.5)*cos(pa),mix(R1,R2,.5)*sin(pa),.1*r1),sdf.x);float da;dspline3(y,vec3 (mix(R1,R2,.5)*cos(pa),mix(R1,R2,.5)*sin(pa),.1*r1),vec3 (R2*cos(pa),R2*sin(pa),.1*r1),vec3 (R2*cos(pa),R2*sin(pa),.1-.4*r0),da);sdf.x=min(sdf.x,da);stroke(sdf.x,.25*mix(.02,.05,.5+.5*n2),sdf.x);sdf.y=2.;add(sdf,vec2 (length(y-vec3 (R2*cos(pa),R2*sin(pa),.1-.4*r0))-.01,3.),sdf);float fa;lfnoise(4.*a,fa);dvoronoi(a,fn,vind2);fa=x.z+.4+.1*mix((v+fn),fa,.5);add(sdf,vec2 (fa,4.),sdf);smoothmin(sdf.x,fa,.1,sdf.x);}void normal(in vec3 x,out vec3 n,in float dx);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}float nan;void vs(in vec3 x,out vec2 sdf){vec2 vi;dsmoothvoronoi(3.*(x.xy+.2*iTime*c.yx),sdf.x,vi);sdf.x=x.z-.1-.2*sdf.x;}void mainImage(out vec4 fragColor,in vec2 fragCoord){vec2 uv=(fragCoord-.5*iResolution.xy)/iResolution.y,s;scale(iScale);uv*=2.;vec3 col=c.yyy,o=c.yzx,r=c.xyy,u=normalize(c.yxx),t=c.yyy,dir,n,x;int N=400,i;t=uv.x*r+uv.y*u;dir=normalize(t-o);float d=-(o.z-.3)/dir.z;for(i=0;i<N;++i){x=o+d*dir;vs(x,s);if(s.x<1.e-4)break;d+=s.x;}float v1,rar,dx=1.e-3;vec2 vi1,na;vec3 cv,l;x=o+d*dir;dsmoothvoronoi(3.*(x.xy+.2*iTime*c.yx),v1,vi1);rand(vi1,rar);cv=mix(c.yyy,vec3 (.23,.23,.23),rar);v1=abs(v1)-.01;cv=mix(cv,c.yyy,sm(v1));v1=abs(v1-.01)-.005;cv=mix(cv,c.xxx,sm(v1));vs(x,s);vs(x+dx*c.xyy,na);n.x=na.x;vs(x+dx*c.yxy,na);n.y=na.x;vs(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);l=normalize(x+c.yyx);cv=.2*cv+.2*cv*abs(dot(l,n))+.4*cv*pow(abs(dot(reflect(-l,n),dir)),3.);cv=mix(cv,1.5*vec3 (0.76,0.20,0.13),smoothstep(0.858,1.02,dot(n,c.yyx)));dir=refract(dir,n,.98);for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=min(s.x,3.e-2);}if(i<N){normal(x,n,1.e-4);if(s.y==3.){l=normalize(x+c.yyx);float r;col=mix(c.xxx,vec3 (0.76,0.20,0.13),.8);col=.2*col+.2*col*abs(dot(l,n))+.8*col*pow(abs(dot(reflect(-l,n),dir)),2.);}else if(s.y==4.){l=normalize(x+c.yyx);float r;rand(vind+vind2,r);col=mix(.023*c.xxx,vec3 (0.76,0.40,0.23),r);col=.2*col+.2*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);stroke(v,.01,v);stroke(fn,.01,fn);col=mix(col,c.yyy,sm(v));col=mix(col,col*col,sm(fn));}if(s.y==2.){l=normalize(x+c.yyx);float r;rand(vind,r);col=mix(mix(vec3 (0.76,0.20,0.23),vec3 (.18,.32,.13),r),vec3 (0.23,0.23,0.23),clamp((x.z)/r1/.1,0.,1.));col=.2*col+.2*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);col=mix(col,5.*col,.25*n.x*n.x);}}col*=3.6;col*=col;if(s.y!=3.){col=mix(length(col)/sqrt(3.)*c.xxx,col,.3);}col=mix(col,cv,.8);col=mix(col,c.yyy,smoothstep(1.,5.,d));col*=mix(1.,15.,mix(.28,.88,0.*iScale));col*=col;col=mix(col,.01*col,smoothstep(-.6,-1.,uv.y));fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *groundboxes_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.0,0.0,-1.0);float a=1.0;float iScale,nbeats;void rand(in vec2 x,out float n);void lfnoise(in vec2 t,out float n);void dspline3(in vec3 x,in vec3 p0,in vec3 p1,in vec3 p2,out float ds);void dbox3(in vec3 x,in vec3 b,out float d);void dlinesegment3(in vec3 x,in vec3 p1,in vec3 p2,out float d);void stroke(in float d0,in float s,out float d);void zextrude(in float z,in float d2d,in float h,out float d);void scale(out float s);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void smoothmin(in float a,in float b,in float k,out float dst);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);void dvoronoi(in vec2 x,out float d,out vec2 z);void dbox(in vec2 x,in vec2 b,out float d);void rot3(in vec3 p,out mat3 rot);vec2 ind=c.yy;void scene(in vec3 x,out vec2 sdf){float d,s=.1;sdf=c.xy;for(float size=.0;size<=.5;size+=.025){dbox3(x,size*c.xxx,d);stroke(d,.001,d);vec2 sda=vec2 (d,3.+size);float n;vec3 y=mod(x,.125*size)-.5*.125*size,yi=(x-y)/size;ind=yi.xy+yi.yz+yi.xz;lfnoise(3.6*ind+15.*size-1.1*iTime,n);if(n>-.3){dbox3(y,.25*size*c.xxx,d);sda.x=max(sda.x,-d);}add(sdf,sda,sdf);}}void normal(in vec3 x,out vec3 n,in float dx);void colorize(in vec2 x,out vec3 col){col=.5*c.xxx;}void mainImage(out vec4 fragColor,in vec2 fragCoord){a=iResolution.x/iResolution.y;mat3 RR;rot3(.2*iTime*vec3 (1.1,1.4,1.6),RR);scale(iScale);vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0),s;vec3 col=c.yyy,o=RR*c.yyx,r=normalize(c.xyy),u=normalize(c.yxy),t=c.yyy,dir,n,x,size=.301*c.xxx;int N=150,i=0;t=uv.x*r+uv.y*u;t=RR*t;dir=normalize(t-o);float d=0.,ra,ss,inside=0.,dd;for(ss=.5;ss>=.0;ss-=.025){size=ss*c.xxx+1.e-4;vec3 tlo=min((size-o)/dir,(-size-o)/dir);vec2 abxlo=abs(o.yz+tlo.x*dir.yz),abylo=abs(o.xz+tlo.y*dir.xz),abzlo=abs(o.xy+tlo.z*dir.xy);vec4 dn=100.*c.xyyy;dn=mix(dn,vec4 (tlo.x,c.xyy),float (all(lessThan(abxlo,size.yz)))*step(tlo.x,dn.x));dn=mix(dn,vec4 (tlo.y,c.yxy),float (all(lessThan(abylo,size.xz)))*step(tlo.y,dn.x));dn=mix(dn,vec4 (tlo.z,c.yyx),float (all(lessThan(abzlo,size.xy)))*step(tlo.z,dn.x));inside+=.05;if(ss==3.)dd=dn.r;d=dn.r;x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;normal(x,n,5.e-4);lfnoise(x.xy*vec2 (3.,8.)-iTime,ra);vec3 f;float dd=5.e-1;r=RR*c.xyy;f=RR*c.yzy;u=RR*c.yyx;vec3 dp=abs(vec3 (dot(n,r),dot(n,f),dot(n,u)));if(dp.y<dd&&dp.z<dd)n=r;else if(dp.x<dd&&dp.z<dd)n=f;else if(dp.x<dd&&dp.y<dd)n=u;s.y=-1.;vec3 l=normalize(x+.03*normalize(x-o));vec3 c1=col;if(s.y==1.){colorize(x.xy,c1);c1=.1*c1+1.*c1*abs(dot(l,n))+1.5*c1*abs(pow(dot(reflect(x-l,n),dir),2.));}else if(s.y==2.){c1=mix(vec3 (0.99,0.43,0.15),vec3 (0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));vec3 c1=mix(vec3 (0.99,0.43,0.15),vec3 (0.44,0.07,0.66),.5*sin(2.*iScale*ra*x));c1=mix(c1,c1,.5+.5*ra);c1=.3*c.xxx;c1=.1*c1+.4*c1*abs(dot(l,n))+.3*c1*abs(pow(dot(reflect(x-l,n),dir),3.));}else if(s.y>=3.){lfnoise(.1*ind+3.*s.y*c.xx-iTime,c1.x);lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime,c1.y);lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime,c1.z);float na;ind=x.xy+x.yz+x.zx;ind=mod(ind,.01)-.005;ind=x.xy+x.yz+x.zx-ind;rand(ind,na);na=0.;c1=.8+.2*c1;c1*=na;c1=mix(.1,mix(.2,.4,iScale),step(na,.05))*c1+mix(.1,.2,step(na,.95))*c1*abs(dot(l,n))+.5*c1*abs(pow(dot(reflect(x-l,n),dir),2.));}else if(s.y==-1.){lfnoise(.1*ind+3.*s.y*c.xx-iTime,c1.x);lfnoise(.1*ind+3.*s.y*c.xx+1337.-iTime,c1.y);lfnoise(.1*ind+3.*s.y*c.xx+2337.-iTime,c1.z);c1=.8+.2*c1;c1=.1*c1+1.*c1*abs(dot(l,n))+1.5*c1*abs(pow(dot(reflect(x-l,n),dir),3.));vec3 dc;vec3 zz=mod(x,.025)-.5*.025,zi=x-zz;rand(zi.xy+zi.yz+zi.zx,dc.x);rand(zi.xy+zi.yz+zi.zx+1337.,dc.y);rand(zi.xy+zi.yz+zi.zx+2337.,dc.z);float da;dbox3(zz,.01*c.xxx,da);stroke(da,.001,da);c1=mix(c1,1.2*c1*c1+.2*dc,sm(da));stroke(da-.002,.001,da);c1=mix(c1,1.6*c1*c1,sm(da));}else if(s.y==-2.){float s=.05;c1=.5*c.xxx;vec2 dd=mod(x.xz,s)-.5*s;stroke(dd.x,.005,dd.x);stroke(dd.y,.005,dd.y);c1=mix(c1,c.xxx,sm(min(dd.x,dd.y)));}col=mix(col,c1,mix(.2,.6,iScale));col=.9*col;}if(s.x>1.e-4&&uv.y<-2.e-4){d=-(o.y+.375)/dir.y;x=o+d*dir;scene(x,s);i=N;}else x=o+d*dir;if(s.y==-1.){vec3 c1=mix(vec3 (0.99,0.43,0.15),vec3 (0.44,0.07,0.66),.5+.5*sin(2.*iScale*ra*x));col=mix(col,c1,.3+.5*inside);}col=mix(col,c.yyy,smoothstep(2.,3.,d));col*=1.2;col*=col;col*=col*col;fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *graffiti_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.0,0.0,-1.0);float a=1.0;float iScale,nbeats;void rand(in vec2 x,out float n);void lfnoise(in vec2 t,out float n);void mfnoise(in vec2 x,in float d,in float b,in float e,out float n);void dtriangle(in vec2 p,in vec2 p0,in vec2 p1,in vec2 p2,out float dst);void dbox(in vec2 x,in vec2 b,out float d);void dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d);void stroke(in float d0,in float s,out float d);void dvoronoi(in vec2 x,out float d,out vec2 z);void rot3(in vec3 phi,out mat3 R);void scale(out float s);void graf(in vec2 x,out float d){x.y*=.7;float size=.4,n,da;vec2 y=vec2 (mod(x.x,size)-.5*size,x.y),yi=(x-y)/size,x1,x2,x3;dbox(y,vec2 (.75,.75)*size,d);rand(yi,n);x1=vec2 (-.5+.02+n*.96,.75)*size,x2=vec2 (.5-.02-n*.96,-.75)*size;x1.x=floor(5.*x1.x)/5.;x2.x=floor(5.*x2.x)/5.;x1.x=max(x1.x,-.4*size);x1.x=min(x1.x,.4*size);x2.x=max(x2.x,-.4*size);x2.x=min(x2.x,.4*size);dlinesegment(y,x1,x2,da);stroke(da,.02,da);d=max(d,-da);rand(yi+1337.,n);x1=vec2 (-.55+n,.75)*size*1.05,x2=vec2 (.5-n,.75)*size*1.05,x3=.75*(.8-n)*size*1.05*c.yx-.01*c.yx;x1=round(15.*x1)/15.;x2=round(15.*x2)/15.;x3=round(15.*x3)/15.;dtriangle(y,x1,x2,x3,da);d=max(d,-da);rand(yi+2337.,n);x1=vec2 (-.5+n,-.75)*size*1.05,x2=vec2 (.55-n,-.75)*size*1.05,x3=-.75*(.8-n)*size*1.05*c.yx+.01*c.yx;x1=round(15.*x1)/15.;x2=round(15.*x2)/15.;x3=round(15.*x3)/15.;dtriangle(y,x1,x2,x3,da);d=max(d,-da);}void zextrude(in float z,in float d2d,in float h,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);mat3 R;void scene(in vec3 x,out vec2 sdf){x=R*x;x.x+=.3*iTime;x*=2.;vec3 n;lfnoise(x.x*c.xx-iTime,n.x);lfnoise(2.*x.x*c.xx-iTime-1337.,n.y);lfnoise(x.x*c.xx+2.*iTime-2337.,n.z);x.yz+=.1*vec2 (cos(x.x),sin(x.x))*n.xy;mat3 RR;rot3(1.3*mix(.2,1.5,.5+.5*n.x)*n.z*c.xyy,RR);x=RR*x;x.z=abs(x.z);float d,da,db;graf(x.xy,d);stroke(d+mix(.01,.04,iScale),mix(.01,.04,iScale),da);float v;vec2 ind;dvoronoi(12.*x.xy,v,ind);zextrude(x.z,-d,.1-.1*v,d);sdf=vec2 (d,1.);float modsize=.025,y=mod(d-.3-.02*iTime,modsize)-.5*modsize,yi=(d-y)/modsize;float na;lfnoise(2.*yi*c.xx-.3*iTime,na);zextrude(x.z-.05*na,-y,mix(0.,.05+.05*na,iScale),d);stroke(d,.035,d);zextrude(x.z,-da,.25,da);add(sdf,vec2 (da,1.),sdf);lfnoise(5.*x.xy,da);mfnoise(x.xy,32.,422.,.45,db);da=.5*(db+da);sdf.x-=.001*da;stroke(da,.1,da);sdf.x-=.005*da;add(sdf,vec2 (d,1.),sdf);add(sdf,vec2 (x.z+.25,1.),sdf);float xa=mix(x.x+3.*a,x.x-3.*a,clamp(iTime/3.,0.,1.));xa=mix(xa,-xa,clamp((iTime-3.)/6.,0.,1.));sdf.x+=mix(0.,2.,smoothstep(-.5,.5,xa));}void normal(in vec3 x,out vec3 n,in float dx);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void colorize(in vec2 x,out vec3 col){x.x+=.3*iTime;x*=2.;float n;lfnoise(x.x*c.xx-iTime,n);x.y+=.3*cos(x.x)*n;float d;graf(x,d);col=mix(col,mix(mix(vec3 (0.85,0.87,0.89),c.xxx,step(50.,iTime)),mix(vec3 (0.04,0.18,0.24),vec3 (0.00,0.20,0.36),step(50.,iTime)),clamp(abs(x.y/2.),0.,1.)),sm(d-.2));col=mix(col,mix(vec3 (1.00,0.40,0.39),vec3 (0.00,0.67,0.91),step(50.,iTime)),sm(d));float da=d;stroke(d+mix(.01,.03,iScale),mix(.01,.04,iScale),d);col=mix(col,1.4*col,sm(d));stroke(d,.001,d);col=mix(col,1.3*col,sm(d));if(da<.02&&da>-.02){lfnoise(5.*x,da);mfnoise(x,32.,422.,.45,d);d=.5*(d+da);col=mix(col,vec3 (0.27,0.27,0.27),sm(d));stroke(d,.1,d);col=mix(col,1.5*col,sm(d));}col*=mix(1.,1.6,iScale);}void mainImage(out vec4 fragColor,in vec2 fragCoord){a=iResolution.x/iResolution.y;scale(iScale);vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0),s;if(iTime>71.)rot3(mix(pi/4.*c.xyy,7.*pi/4.*c.xxy,smoothstep(71.,86.,iTime)),R);else R=mat3 (1.);float sc2=0.,sc3=0.;vec3 col=c.yyy,o=mix(1.,.5,smoothstep(0.,5.,clamp(iTime-71.,0.,5.)))*c.yzx,r=c.xyy,u=normalize(c.yxx),t=c.yyy,dir,n,x;int N=250,i;t=uv.x*r+uv.y*u;dir=normalize(t-o);vec3 c1;float d=-(o.z-.35)/dir.z;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;if(x.z<-.15){i=N;break;}d+=min(s.x,5.e-2);}if(i<N){normal(x,n,1.e-2);if(s.y==1.){vec3 l=normalize(x+.5*c.yzx);colorize(x.xy,c1);c1=.1*c1+1.*c1*abs(dot(l,n))+1.5*c1*abs(pow(dot(reflect(x-l,n),dir),2.));}else if(s.y==2.){vec3 l=normalize(x+c.xzx);float r;lfnoise(x.xy,r);c1=mix(vec3 (0.99,0.43,0.15),vec3 (0.44,0.07,0.66),sin(2.*iScale*r*x));c1=.1*c1+.8*c1*abs(dot(l,n))+6.5*c1*abs(pow(dot(reflect(x-l,n),dir),3.));}col=c1;}col*=col*col;col=mix(col,c.yyy,clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));if(length(col)<.001){float v,ra,v2;vec2 ind,ind2;lfnoise(iTime*c.xx,ind2.x);lfnoise(iTime*c.xx-1337.,ind2.y);dvoronoi(12.*(uv-.03*ind2),v,ind);rand(ind,ra);stroke(-v,.05,v);v=-v;col=mix(col,.3*ra*mix(.5*vec3 (1.00,0.40,0.39),.05*c.xxx,clamp(tanh(1.5*length(uv)),0.,1.)),sm(v));col*=mix(13.,1.,smoothstep(0.,.5,clamp((iTime-6.),0.,1.)));}fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *greet_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const vec3 c=vec3 (1.0,0.0,-1.0);const float pi=acos(-1.);float iScale;void scale(out float s);void rand(in vec2 x,out float n);void hash31(in float p,out vec3 d);void lfnoise(in vec2 t,out float n);void dbox(in vec2 x,in vec2 b,out float d);void stroke(in float d0,in float s,out float d);void zextrude(in float z,in float d2d,in float h,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);void smoothmin(in float a,in float b,in float k,out float dst);void dbox3(in vec3 x,in vec3 b,out float d);void dschnappsgirls(in vec2 x,out float d);void dspacepigs(in vec2 x,out float d);void dkewlers(in vec2 x,out float d);void dfarbrausch(in vec2 x,out float d);void dhaujobb(in vec2 x,out float d);void dmercury(in vec2 x,out float d);vec2 ind,indc;void scene(in vec3 x,out vec2 sdf){float d;d=mix(mix(mix(mix(mix(mix(mix(0.,0.14173228346456693,smoothstep(0.,0.341334,iTime)),0.14173228346456693+.25,smoothstep(0.341334,2.276818,iTime)),0.14173228346456693+.5,smoothstep(2.276818,4.151818,iTime)),0.14173228346456693+.75,smoothstep(4.151818,4.151818+1*1.8182,iTime)),0.14173228346456693+1.,smoothstep(4.151818+1*1.8182,4.151818+2*1.8182,iTime)),0.14173228346456693+1.25,smoothstep(4.151818+2*1.8182,4.151818+3*1.8182,iTime)),0.,smoothstep(4.151818+3*1.8182,4.151818+5*1.8182,iTime));x.z-=d;dbox3(x,vec3 (.1,.1,1.e3),d);sdf=vec2 (-d,2.);float distortion;lfnoise(5.2e2*x.yz,distortion);float tsize=.005,dy=mod(x.y,tsize)-.5*tsize,yi=(x.y-dy)/tsize,zpar=x.z+mix(0.,.5*tsize,mod(yi,2.)),dz=mod(zpar,tsize)-.5*tsize,zi=(zpar-dz)/tsize;dbox3(vec3 (abs(x.x)-.1,dy,dz),vec3 (.0005+.00001*distortion,.39*tsize*c.xx),d);add(sdf,vec2 (d,3.),sdf);smoothmin(sdf.x,d,.001,sdf.x);ind=vec2 (yi,zi);tsize=.025;dz=mod(x.z,tsize)-.5*tsize;float dx=mod(x.x,tsize)-.5*tsize;zi=(x.z-dz)/tsize;float xi=(x.x-dx)/tsize;dbox3(vec3 (dx,abs(x.y)-.1,dz),vec3 (.48*tsize,.0005,.48*tsize),d);add(sdf,vec2 (d,4.),sdf);smoothmin(sdf.x,d,.002,sdf.x);indc=vec2 (xi,zi);tsize=.25;float tw=.0005;dz=mod(x.z-.5*tsize,tsize)-.5*tsize;zi=round((x.z-dz)/tsize);zi=mod(zi,6.);if(zi<.5)dmercury(20.*x.xy,d);else if(zi<1.5)dhaujobb(20.*x.xy,d);else if(zi<2.5)dfarbrausch(20.*x.xy,d);else if(zi<3.5)dkewlers(20.*x.xy,d);else if(zi<4.5)dspacepigs(20.*x.xy,d);else if(zi<5.5)dschnappsgirls(20.*x.xy,d);stroke(d/20.,tw,d);zextrude(dz,-d,.005,d);if(zi==0.)add(sdf,vec2 (d,6.),sdf);else add(sdf,vec2 (d,5.),sdf);}void normal(in vec3 x,out vec3 n,in float dx);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void mainImage(out vec4 fragColor,in vec2 fragCoord){vec2 uv=(fragCoord-.5*iResolution.xy)/iResolution.y,s;scale(iScale);float d0;vec3 o0,dir0;vec3 col=c.yyy,o=c.yyx,r=c.xyy,u=c.yxy,t=c.yyy,dir,n,x;int N=400,i;t=uv.x*r+uv.y*u;dir=normalize(t-o);float d=0.;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=min(s.x,3.e-2);}if(i<N){normal(x,n,5.e-4);vec3 l=x+.1*n;o0=o;d0=d;dir0=dir;if(s.y==2.){col=.23*c.xxx;col=.2*col+.2*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);}else if(s.y==3.){float r;rand(ind,r);col=.2*vec3 (0.02,0.11,0.24)+.2*vec3 (0.25,0.75,0.85)*abs(dot(l,n))+mix(.5,.1,r)*vec3 (0.45,0.69,0.76)*pow(abs(dot(reflect(-l,n),dir)),2.);d0=d;o0=o;dir0=dir;d=1.e-2;o=x;dir=reflect(dir,n);for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=min(s.x,3.e-1);}normal(x,n,5.e-4);vec3 l=x+.1*n;vec3 c1;if(s.y==2.){c1=.23*c.xxx;c1=.2*c1+.2*c1*abs(dot(l,n))+.6*c1*pow(abs(dot(reflect(-l,n),dir)),3.);}else if(s.y==3.){float r;rand(ind,r);c1=.2*vec3 (0.02,0.11,0.24)+.2*vec3 (0.25,0.75,0.85)*abs(dot(l,n))+mix(.5,.1,r)*vec3 (0.45,0.69,0.76)*pow(abs(dot(reflect(-l,n),dir)),2.);}else if(s.y==4.){float r;rand(indc,r);c1=.2*.2*c.xxx+.2*.5*mix(c.xxx,vec3 (0.02,0.11,0.24),step(0.,-x.y))*abs(dot(l,n))+mix(.5,.1,r)*.8*c.xxx*pow(abs(dot(reflect(-l,n),dir)),2.);}else if(s.y==5.){c1=.2*.2*c.xyy+.2*.5*mix(c.xyy,vec3 (0.24,0.11,0.024),step(0.,-x.y))*abs(dot(l,n))+.8*vec3 (.8,.3,.2)*pow(abs(dot(reflect(-l,n),dir)),2.);c1=mix(c1,c.xxx,.1);}else if(s.y==6.){c1=.2*.2*c.xyy+.2*.5*mix(c.xyy,vec3 (0.24,0.11,0.024),step(0.,-x.y))*abs(dot(l,n))+mix(.5,1.1,.5+.5*sin(6.*iTime))*mix(vec3 (.8,.3,.2),vec3 (1.,.3,.6),.5+.5*sin(6.*iTime))*pow(abs(dot(reflect(-l,n),dir)),2.);c1=mix(c1,c.xxx,.1);}col=mix(col,c1,.5);}else if(s.y==4.){float r;rand(indc,r);col=.2*.2*c.xxx+.2*.5*mix(c.xxx,vec3 (0.02,0.11,0.24),step(0.,-x.y))*abs(dot(l,n))+mix(.5,.1,r)*.8*c.xxx*pow(abs(dot(reflect(-l,n),dir)),2.);}else if(s.y==5.){col=.2*.2*c.xyy+.2*.5*mix(c.xyy,vec3 (0.24,0.11,0.024),step(0.,-x.y))*abs(dot(l,n))+.8*vec3 (.8,.3,.2)*pow(abs(dot(reflect(-l,n),dir)),2.);col=mix(col,c.xxx,.1);}else if(s.y==6.){col=.2*.2*c.xyy+.2*.5*mix(c.xyy,vec3 (0.24,0.11,0.024),step(0.,-x.y))*abs(dot(l,n))+mix(.5,1.1,.5+.5*sin(6.*iTime))*mix(vec3 (.8,.3,.2),vec3 (1.,.3,.6),.5+.5*sin(6.*iTime))*pow(abs(dot(reflect(-l,n),dir)),2.);col=mix(col,c.xxx,.1);}}float da;{float dd=mix(mix(mix(mix(mix(mix(mix(0.,0.14173228346456693,smoothstep(0.,0.341334,iTime)),0.14173228346456693+.25,smoothstep(0.341334,2.276818,iTime)),0.14173228346456693+.5,smoothstep(2.276818,4.151818,iTime)),0.14173228346456693+.75,smoothstep(4.151818,4.151818+1*1.8182,iTime)),0.14173228346456693+1.,smoothstep(4.151818+1*1.8182,4.151818+2*1.8182,iTime)),0.14173228346456693+1.25,smoothstep(4.151818+2*1.8182,4.151818+3*1.8182,iTime)),0.,smoothstep(4.151818+3*1.8182,4.151818+5*1.8182,iTime));da=-(o0.z-dd)/dir0.z;x=o0+da*dir0;float tsize=.25;float tw=.0005;float dz=mod(x.z-.5*tsize,tsize)-.5*tsize;float zi=round((x.z-dz)/tsize);zi=mod(zi,6.);if(da>0.){dmercury(20.*x.xy,da);stroke(da,tw,da);col=mix(col,mix(col,vec3 (1.,.3,.6),.5),sm(da/50.*(.6+.4*sin(6.*iTime))));}}col=mix(col,0.*.23*c.xxx*vec3 (0.76,0.20,0.13),smoothstep(1.,5.,d));col*=3.6;col*=col;col=clamp(col,0.,1.);col=mix(col,c.yyy,smoothstep(4.151818+3*1.8182,13.,iTime));fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *evoke_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.0,0.0,-1.0);float a=1.0;float iScale,nbeats;void rand(in vec2 x,out float n);void dbox(in vec2 x,in vec2 b,out float d);void stroke(in float d0,in float s,out float d);void dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d);void lfnoise(in vec2 t,out float n);void dvoronoi(in vec2 x,out float d,out vec2 z);void scale(out float s);void devoke(in vec2 x,out float d){x.x+=.225;x*=1.1;d=length(x+.35*c.xy)-.1;stroke(d,.06,d);float da;dbox(x+.1*c.xy,vec2 (.05,.25),da);d=min(d,da);x=2.*x-vec2 (.4,-.2);dbox(x-.35*c.yx,vec2 (.4,.35),da);d=min(d,da);dbox(x-.7*c.yx,vec2 (.2,.2),da);d=max(d,-da);dbox(x-.25*c.yx,vec2 (.2,.05),da);d=max(d,-da);dbox(x+.1*c.yx,vec2 (.1,.2),da);d=min(d,da);dbox(x+.2*c.yx,vec2 (.4,.1),da);d=min(d,da);x=.5*(x+vec2 (.4,-.2));dbox(x-.9*c.xy,vec2 (.05,.25),da);d=min(d,da);dbox(vec2 (x.x-.7,abs(x.y)-.2),vec2 (.2,.05),da);d=min(d,da);dbox(x-.7*c.xy,vec2 (.2,.05),da);d=min(d,da);dbox(vec2 (x.x-.95,x.y+.2),vec2 (.05,.05),da);d=min(d,da);stroke(d,.001,d);}void dstripe(in vec2 x,out float d){dlinesegment(x-a*mix(-.4*c.xy,.4*c.xy,clamp(iTime/6.,0.,1.)),-.5*c.yx,.5*c.yx,d);d-=.005;float dd;vec2 vi;dvoronoi(5.*x,dd,vi);vi=x-vi/5.;dd=abs(length(vi)-.002)-.001;d=min(d,dd);stroke(d,.001,d);d=mix(1.,d,clamp(iTime,0.,1.));}float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void mainImage(out vec4 fragColor,in vec2 fragCoord){a=iResolution.x/iResolution.y;scale(iScale);vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0),s;vec3 col=vec3 (0.20,0.01,0.14),o=c.yyx,r=c.xyy,u=c.yxy,t=c.yyy,dir,n,x;float d,i,ra;t=uv.x*r+uv.y*u;dir=normalize(t-o);for(i=1.4;i>=0.;i-=.01){lfnoise(102.*i*c.xx-mix(102.,0.,smoothstep(0.,1.,clamp(iTime-6.,0.,1.)))*iTime,ra);ra=.5+.5*ra;d=-(o.z-.2+i)/dir.z;x=o+d*dir;float da;dstripe(x.xy,da);devoke(x.xy,s.x);s.x=mix(da,s.x,smoothstep(0.,1.,clamp(iTime-5.5,0.,1.)));s.x-=.01*iScale;if(ra<.5){vec3 c1=mix(mix(vec3 (0.75,0.24,0.31),vec3 (1.00,0.87,0.57),smoothstep(1.25,1.4,1.4-i)),vec3 (0.20,0.01,0.14),i/1.4);if(iTime>6.)col=mix(col,c1,sm(s.x));col=mix(col,mix(1.1,1.,clamp(-iScale+smoothstep(0.,1.,clamp(iTime-6.,0.,1.)),0.,1.))*mix(col,vec3 (.7,.45,.3),mix(.02,.1,iScale)),sm(s.x/64.));}}col=mix(col,c.yyy,clamp((d-2.-(o.z-.2)/dir.z)/4.,0.,1.));vec3 c1=c.yyy;float v,v2;vec2 ind,ind2;lfnoise((iTime-12.)*c.xx,ind2.x);lfnoise((iTime-12.)*c.xx-1337.,ind2.y);dvoronoi(12.*(uv-.03*ind2),v,ind);rand(ind,ra);stroke(-v,.05,v);v=-v;c1=mix(c1,.3*ra*mix(.5*vec3 (1.00,0.40,0.39),.05*c.xxx,clamp(tanh(1.5*length(uv)),0.,1.)),sm(v));c1*=mix(1.,13.,smoothstep(0.,1.,clamp((iTime-11.),0.,1.)));col=mix(col,c1,smoothstep(0.,1.,clamp((iTime-11.),0.,1.)));fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *canal_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;float nbeats;float iScale;const vec3 c=vec3 (1.0,0.0,-1.0);const float pi=acos(-1.);void scale(out float s);void rand(in vec2 x,out float n);void hash31(in float p,out vec3 d);void lfnoise(in vec2 t,out float n);void stroke(in float d0,in float s,out float d);void zextrude(in float z,in float d2d,in float h,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);void smoothmin(in float a,in float b,in float k,out float dst);void dsmoothvoronoi(in vec2 x,out float d,out vec2 z);void rot3(in vec3 phi,out mat3 R);mat3 R;vec2 ind;void scene(in vec3 x,out vec2 sdf){x=R*x;x.z-=mix(1.3,-1.3,step(156.,iTime))*iTime;float dx,d,v;lfnoise(.5*x.z*c.xx,dx);x.xy-=.2*dx*c.xy;float phi=atan(x.y,x.x);dsmoothvoronoi(2.*vec2 (mod(phi+pi/4.,2.*pi),x.z),v,ind);stroke(v,.01,v);d=length(x.xy)-mix(1.,1.1,smoothstep(.0,.2,v));zextrude(length(x.xy)-1.0,d,.05,d);d-=.05;sdf=vec2 (d,1.);dsmoothvoronoi(8.*vec2 (mod(phi+pi/4.,2.*pi),x.z),v,ind);stroke(v,.02,v);d=length(x.xy)-mix(1.1,1.2,smoothstep(.0,.2,v));zextrude(length(x.xy)-1.11,d,.01,d);d-=.1;add(sdf,vec2 (d,1.),sdf);smoothmin(d,sdf.x,.1,sdf.x);float dy=0.;for(float i=1.;i<=15.;i+=1.){float f,a;vec2 dir;float n;lfnoise(i*c.xx-n,f);f=.5+.5*f;lfnoise(i*c.xx+1337.-n,a);a=.5+.5*a;lfnoise(i*c.xx+2337.-2.*n,dir.x);lfnoise(i*c.xx+3337.-3.*n,dir.y);dir=mix(c.yx-.2*c.xy,x.yx+.2*c.xy,2.*dir);dir=normalize(dir);float dya=pow(1.01,f)*a*sin(-2.e-3*2.*pi*pow(1.95,abs(f+.01*a))*(1.*f-.01*a)*iTime-2.e-4*2.*pi*pow(1.99,abs(i-.1*a))*dot(dir,vec2 (.5,4.)*(2.*(x.xz+1.3*(iTime))*c.yx)));dy+=2.*pow((dya+1.)/2.,4.)-1.;}dy=.4*dy;add(sdf,vec2 (x.y+.4+.001*dy,2.),sdf);}void normal(in vec3 x,out vec3 n,in float dx);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void mainImage(out vec4 fragColor,in vec2 fragCoord){vec2 uv=(fragCoord-.5*iResolution.xy)/iResolution.y,s;float phi=mix(0.,.3*iTime,smoothstep(5.,6.,iTime))*step(iTime,156.),co=cos(phi),si=sin(phi);float phi2=mix(mix(mix(0.,pi/4.,smoothstep(5.,6.,iTime)),-pi/4.,smoothstep(9.,10.,iTime)),pi/2.,smoothstep(13.,14.,iTime))*step(iTime,156.),co2=cos(phi2),si2=sin(phi2);rot3(mod(vec3 (phi,phi2,0.),pi),R);uv=mix(uv,mat2 (co,si,-si,co)*uv,step(156.,iTime));uv.y=mix(uv.y,-uv.y,step(156.,iTime));scale(iScale);float dx,dx2,d0;lfnoise(-.5*1.3*iTime*c.xx,dx);lfnoise(-.5*1.3*(iTime+1.e-3)*c.xx,dx2);vec3 col=c.yyy,o=c.yyx+.1*c.yxy+.2*dx*c.xyy,r=c.xyy,u=c.yxy,t=c.yyy+.2*dx2*c.xyy,dir,n,x;int N=300,i,a=0;t=uv.x*r+uv.y*u;dir=normalize(t-o);float d=.5/length(dir.xy);for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;if(length(x.xy-.4*dx*c.xy)>1.5){col=c.yyy;i=N;break;}d+=min(s.x,2.e-2);}if(i<N){normal(x,n,5.e-4);vec3 l=normalize(x+.5*n);if(s.y==1.){col=vec3 (0.25,0.25,0.25);vec3 c0=col;col=.2*col+.2*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);col=mix(col,1.2*c0,smoothstep(0.658,1.02,1.-abs(dot(n,c.yyz))));vec3 c2=mix(mix(vec3 (0.76,0.20,0.23),vec3 (0.07,0.64,0.29),step(166.,iTime)),vec3 (0.96,0.7,0.423),step(174.,iTime));col=mix(col,c2,smoothstep(0.658,1.02,abs(dot(n,c.yyz))));}else if(s.y==2.){col=.3*c.xxx;col=.2*col+.2*col*abs(dot(l,n))+.3*col*pow(abs(dot(reflect(-l,n),dir)),2.);N=50;o=x;dir=reflect(dir,n);d0=d;d=1.e-2;vec3 c1=c.yyy;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;if(length(x.xy)>1.5){c1=c.yyy;i=N;break;}d+=s.x;}if(i<N){normal(x,n,5.e-4);vec3 l=normalize(x+.5*n);if(s.y==1.){col=vec3 (0.25,0.25,0.25);vec3 c0=col;col=.2*col+.2*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);col=mix(col,1.2*c0,smoothstep(0.658,1.02,1.-abs(dot(n,c.yyz))));vec3 c2=mix(mix(vec3 (0.76,0.20,0.23),vec3 (0.07,0.64,0.29),step(166.,iTime)),vec3 (0.96,0.7,0.423),step(174.,iTime));col=mix(col,c2,smoothstep(0.658,1.02,abs(dot(n,c.yyz))));}}col=mix(col,c1,.35);col=mix(col,.1*c.yxx,.3);}}col=mix(col,c.yyy,smoothstep(2.,22.,d+d0));float nn;lfnoise(12.*(x.z-1.3*iTime)*c.xx,nn);col*=mix(1.1,2.6,mix(.5+.5*nn,1.,0.*iScale));col*=col;col=clamp(col,0.,1.);if(col==c.xxx)col=c.yyy;fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *text_source = "#version 130\nuniform float iFontWidth,iTime;uniform vec2 iResolution;uniform sampler2D iChannel0,iFont;uniform float iFSAA;out vec4 gl_FragColor;const vec3 c=vec3 (1.,0.,-1.);const float pi=acos(-1.);float a;void rand(in vec2 x,out float num);void lfnoise(in vec2 t,out float num);void rshort(in float off,out float val);void rfloat(in float off,out float val);void dbox(in vec2 x,in vec2 b,out float dst);void dcircle(in vec2 x,out float d);void dlinesegment(in vec2 x,in vec2 p1,in vec2 p2,out float d);void drhomboid(in vec2 x,in vec2 b,in float tilt,out float dst);void dcirclesegment(in vec2 x,in float r,in float p0,in float p1,out float d);void stroke(in float d0,in float s,out float d);void dglyph(in vec2 x,in float ordinal,in float size,out float dst);void dstring(in vec2 x,in float ordinal,in float size,out float dst);void dfloat(in vec2 x,in float num,in float size,out float dst);void smoothmin(in float a,in float b,in float k,out float dst);void dint(in vec2 x,in float num,in float size,in float ndigits,out float dst);void dtime(in vec2 x,in float num,in float size,out float dst);void window(in vec2 x,in vec2 size,in vec3 bg,in float title_index,out vec4 col);void progressbar(in vec2 x,in float width,in float progress,out vec4 col);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void dvoronoi(in vec2 x,out float d,out vec2 z);void colorize(in vec2 x,out vec3 col){vec3 c1;vec2 ind,xv,xi;float d,vs=16.,n,size=.1,xix=mod(x.x,size)-.5*size,xixj=(x.x-xix),ri,rim1,rip1,lines=8.,da,op,s;s=smoothstep(0.,.5,.5-abs(x.y));col=mix(1.e-4*c.xxx,vec3 (0.04,0.18,0.24),s);dvoronoi(vs*x,d,ind);xv=ind/vs-x;lfnoise(vec2 (3.,33.)*ind/vs-3.*iTime*c.xy,n);n=.5+.5*n;d=length(xv)-mix(.0,.35,n)/vs;col=mix(col,n*.5*vec3 (1.00,0.40,0.39),sm(d));d=abs(d-.005)-.002;col=mix(col,(1.-n)*vec3 (0.49,0.71,0.78),sm(d));for(float i=1.;i<9.;i+=1.){rand((9.-i)*c.xx,op);op=.5+.5*round(16.*op)/16.;x+=-.1+.2*op;xix=mod(x.x,size)-.5*size;xixj=(x.x-xix);lfnoise(2.e0*xixj*c.xx+14.*i,ri);lfnoise(2.e0*(xixj+size)*c.xx+14.*i,rip1);lfnoise(2.e0*(xixj-size)*c.xx+14.*i,rim1);float h=.2;ri=h*round(lines*ri)/lines;rip1=h*round(lines*rip1)/lines;rim1=h*round(lines*rim1)/lines;{dlinesegment(vec2 (xix,x.y),vec2 (-.5*size,mix(ri,rim1,.5)),vec2 (-.25*size,ri),d);dlinesegment(vec2 (xix,x.y),vec2 (-.25*size,ri),vec2 (.25*size,ri),da);d=min(d,da);dlinesegment(vec2 (xix,x.y),vec2 (.25*size,ri),vec2 (.5*size,mix(ri,rip1,.5)),da);d=min(d,da);stroke(d,.002+.002*op,d);col=mix(col,op*(1.-n)*vec3 (1.00,0.40,0.39),sm(d));lfnoise(8.*xixj*c.xx-3.*iTime*c.xy+14.*i,n);n=.5+.5*n;d=length(vec2 (xix,x.y-ri))-mix(.0,.35,n)/vs;c1=mix(vec3 (1.00,0.40,0.39),vec3 (0.85,0.87,0.89),n);col=mix(col,op*(1.-n)*c1,sm(d));stroke(d-.009,(1.-n)*.005,d);c1*=2.4;col=mix(col,op*(1.-n)*c1,sm(d));}x-=-.1+.2*op;}lfnoise(3.*x.xy-vec2 (1.,.1)*iTime,n);stroke(n,.3,n);col=mix(col,1.e-4*c.xxx,n);col=mix(col,.1*col,1.-s);col=mix(col,mix(col,vec3 (1.00,0.40,0.39),mix(.4,.8,.5+.5*x.y/.1)),sm(abs(x.y)-.1));col=mix(col,c.xxx,sm(abs(abs(x.y)-.11)-.001));col=mix(col,col*col,clamp(-x.y/.1,0.,1.));col*=col;}void mainImage(out vec4 fragColor,in vec2 fragCoord){a=iResolution.x/iResolution.y;vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0);float d;vec4 old=c.yyyy,new=c.yyyy;float bound=sqrt(iFSAA)-1.;for(float i=-.5*bound;i<=.5*bound;i+=1.)for(float j=-.5*bound;j<=.5*bound;j+=1.){old.gba+=texture(iChannel0,(fragCoord+vec2 (i,j)*3./max(bound,1.))/iResolution.xy).xyz;}old.gba/=iFSAA;new=old;if(uv.y>.4){float da;dstring((uv-.45*vec2 (-.55*a,1.+4.*.008)),9.,.004,d);dstring((uv-.45*vec2 (-.55*a,1.+2.*.008)),10.,.004,da);d=min(d,da);dstring((uv-.45*vec2 (-.55*a,1.)),11.,.004,da);d=min(d,da);dstring((uv-.45*vec2 (-.55*a,1.-2.*.008)),12.,.004,da);d=min(d,da);dstring((uv-.45*vec2 (-.55*a,1.-4.*.008)),13.,.004,da);d=min(d,da);new.gba=mix(new.gba,mix(new.gba,c.xxx,.5),sm(d));dstring((uv-.45*vec2 (-.85*a,1.)),3.,.02,d);stroke(d-.002,.001,d);new.gba=mix(new.gba,vec3 (1.00,0.40,0.39),sm(d));dtime((uv-.45*vec2 (.975*a,1.05)),iTime+11.,.01,d);new.gba=mix(new.gba,c.xxx,sm(d));dint(uv-.45*vec2 (.975*a,1.0),floor(1.e3*fract(iTime)),.01,4.,d);stroke(d-.001,.0005,d);new.gba=mix(new.gba,c.xxx,sm(d));}if(iTime<0.){new.gba=old.gba;float sc=smoothstep(0.,1.,clamp(iTime+3.,0.,1.))*(1.-smoothstep(0.,1.,clamp(iTime+1.,0.,1.)));dstring((uv-vec2 (-.085,-.3)),3.,.02,d);float da;dstring((uv-vec2 (-.08,-.35)),26.,.02,da);d=min(d,da);new.gba=mix(new.gba,mix(new.gba,c.yyy,sc),sm(d));}else if(iTime<6.){vec2 dx=(.25*a+.3*c.xy)*c.xy;if(iTime<3.){float ind=mix(100000.,2.,clamp(iTime/3.,0.,1)),da;dint(uv+dx*c.xy,ind,.02,6.,d);dstring(uv+dx-2.*9.*.02*c.xy,4.,.02,da);d=min(d,da);}else if(iTime<4.){dint(uv+dx,2.,.02,6.,d);float da;dstring(uv+dx-2.*9.*.02*c.xy,4.,.02,da);d=min(d,da);}else if(iTime<5.){dint(uv+dx+.04*c.yx,1.,.02,6.,d);float da;dint(uv+dx,2.,.02,6.,da);d=min(d,da);dstring(uv+dx-2.*9.*.02*c.xy+.04*c.yx,4.,.02,da);d=min(d,da);}else if(iTime<6.){dint(uv+dx+.08*c.yx,0.,.02,6.,d);float da;dint(uv+dx+.04*c.yx,1.,.02,6.,da);d=min(d,da);dint(uv+dx,2.,.02,6.,da);d=min(d,da);dstring(uv+dx-2.*9.*.02*c.xy+.08*c.yx,4.,.02,da);d=min(d,da);}new.gba=mix(new.gba,mix(new.gba,vec3 (1.00,0.87,0.57),.7),sm(d));stroke(d-.002,.001,d);new.gba=mix(new.gba,c.xxx,sm(d));}else if(iTime<12.&&iTime>7.){float da,db;dbox(vec2 (uv.x+.75,uv.y-.35),vec2 (.013,.035),da);stroke(da,.002,da);dglyph(vec2 (uv.x+.75,uv.y-.35-.02).yx*c.zx,101.,.01,db);da=min(da,db);dglyph(vec2 (uv.x+.75,uv.y-.35),118.,.01,db);da=min(da,db);dglyph(vec2 (uv.x+.75,uv.y-.35+.02).yx*c.zx,107.,.01,db);da=min(da,db);vec2 b=vec2 (uv.x+.75,uv.y-.35+.02)-.01*c.xx-.02*c.xy,b1=mod(b,.02)-.01,b1i=floor((b-b1)/.02);if(abs(b1i.y)<=1.&&b1i.x>=0.&&b1i.x<=10.){lfnoise(b1i-12.*iTime,db);db=97.+mod(floor(26.*(.5+.5*db)),26.);dglyph(b1,db,.008,db);da=min(da,db);}dlinesegment(vec2 (uv.x+.75,uv.y-.35+.06),-.015*c.xy,.25*c.xy,db);stroke(db,.001,db);da=min(da,db);dstring(vec2 (uv.x+.75,uv.y+.35),5.,.015,db);da=min(da,db);dstring(vec2 (uv.x-.2,uv.y+.35),6.,.015,db);float dc;dbox(vec2 (uv.x-.2-.12,uv.y+.35),vec2 (.165,.015),dc);db=max(dc,-db);da=min(da,db);dstring(vec2 (uv.x+.75,uv.y+.4),7.,.015,db);da=min(da,db);dstring(vec2 (uv.x-.2,uv.y+.4),8.,.015,db);dbox(vec2 (uv.x-.2-.12,uv.y+.4),vec2 (.165,.015),dc);db=max(dc,-db);da=min(da,db);new.gba=mix(new.gba,vec3 (0.75,0.24,0.30),sm(da));}else if(iTime<28.){float da=length(uv)-.45,db;dstring((uv+.3*c.xy),2.,.0415,db);db-=.001;da=max(da,-db);da=mix(1.,da,smoothstep(0.,.5,clamp(iTime-18.5,0.,1.))*(1.-smoothstep(0.,.5,clamp(iTime-22.,0.,1.))));new.gba=mix(new.gba,vec3 (1.00,0.40,0.39),sm(da));da=length(uv-.3*c.xx)-.2,db;dstring(2.*(uv+.075*c.xy-.3*c.xx),3.,.0415,db);db-=.001;da=mix(1.,da,smoothstep(0.,.5,clamp(iTime-19.5,0.,1.))*(1.-smoothstep(0.,.5,clamp(iTime-22.,0.,1.))));db=mix(1.,db,smoothstep(0.,.5,clamp(iTime-19.5,0.,1.))*(1.-smoothstep(0.,.5,clamp(iTime-22.,0.,1.))));new.gba=mix(new.gba,vec3 (1.00,0.40,0.39)*vec3 (1.00,0.40,0.39),sm(da));new.gba=mix(new.gba,c.yyy,sm(db));dstring((uv-vec2 (-.75,-.35)).yx*c.xz,18.,.045,db);dstring((uv-vec2 (-.65,-.35)).yx*c.xz,19.,.045,da);db=min(db,da);db=mix(1.,db,smoothstep(0.,.5,clamp(iTime-24.5,0.,1.))*(1.-smoothstep(0.,.5,clamp(iTime-28.,0.,1.))));new.gba=mix(new.gba,mix(new.gba,c.xxx,.8),sm(db));stroke(db-.005,.0005,db);new.gba=mix(new.gba,mix(new.gba,vec3 (1.00,0.40,0.39),.8),sm(db));da=length((uv-vec2 (-.6,-.325)).yx*c.xz)-.1;dstring((uv-vec2 (-.6,-.35)).yx*c.xz,20.,.015,db);da=max(da,-db);da=mix(1.,da,smoothstep(0.,.5,clamp(iTime-25.5,0.,1.))*(1.-smoothstep(0.,.5,clamp(iTime-28.,0.,1.))));new.gba=mix(new.gba,mix(new.gba,c.xxx,.6),sm(da));}else if(iTime<44.){float da,db;dstring((uv-vec2 (-.3,.3)),21.,.1,da);db=abs(mod(uv.x+uv.y,.3)-.15)-.075;vec3 c1=mix(mix(new.gba,.3*c.xxx,.5),c.xxx,sm(db));db=smoothstep(33.,34.,iTime);da=mix(1.,da,db);new.gba=mix(new.gba,c1,sm(da));stroke(da-.02,.001,da);new.gba=mix(new.gba,c.xxx,sm(da));dstring(uv-vec2 (.35,.34),22.,.05,da);dbox(uv-vec2 (.35,.34),vec2 (.15,.06),db);db=max(db,-da);float mx=smoothstep(34.,35.,iTime);db=mix(1.,db,mx);new.gba=mix(new.gba,c.xxx,sm(db));dstring(uv-vec2 (.25,.24),23.,.0277,da);mx=smoothstep(35.,36.,iTime);da=mix(1.,da,mx);new.gba=mix(new.gba,.8*c.xxy,sm(da));dstring(uv-vec2 (.25,.15),22.,.05,da);mx=smoothstep(36.,37.,iTime);da=mix(1.,da,mx);new.gba=mix(new.gba,.8*c.xxy,sm(da));dstring((uv-vec2 (.45,.05)).yx*c.zx,24.,.05,da);dbox((uv-vec2 (.45,-.1)),vec2 (.05,.3),db);db=max(db,-da);mx=smoothstep(37.,38.,iTime);db=mix(1.,db,mx);new.gba=mix(new.gba,c.xxx,sm(db));dstring((uv-vec2 (.6,.1)).yx*c.zx,25.,.1,da);db=smoothstep(38.,39.,iTime);da=mix(1.,da,db);new.gba=mix(new.gba,c1,sm(da));stroke(da-.02,.001,da);new.gba=mix(new.gba,c.xxx,sm(da));}else if(iTime<60.){float da,db;db=abs(mod(uv.x+uv.y,.3)-.15)-.075;vec3 c1=mix(mix(new.gba,vec3 (0.93,0.36,0.44),.5),c.xxx,sm(db));dstring((uv-vec2 (-.3,.3)),38.,.02,da);db=smoothstep(45.,46.,iTime)-smoothstep(50.,51.,iTime);da=mix(1.,da,db);stroke(da-.005,.0025,db);new.gba=mix(new.gba,c.yyy,sm(db));new.gba=mix(new.gba,c1,sm(da-.0025));db=da+.0025;new.gba=mix(new.gba,2.*c1,sm(db));dstring((uv-vec2 (-.3,.25)),39.,.02,da);db=smoothstep(46.,47.,iTime)-smoothstep(51.,52.,iTime);da=mix(1.,da,db);stroke(da-.005,.0025,db);new.gba=mix(new.gba,c.yyy,sm(db));new.gba=mix(new.gba,c1,sm(da-.0025));db=da+.0025;new.gba=mix(new.gba,2.*c1,sm(db));dstring((uv-vec2 (-.3,.2)),40.,.02,da);db=smoothstep(47.,48.,iTime)-smoothstep(52.,53.,iTime);da=mix(1.,da,db);stroke(da-.005,.0025,db);new.gba=mix(new.gba,c.yyy,sm(db));new.gba=mix(new.gba,c1,sm(da-.0025));db=da+.0025;new.gba=mix(new.gba,2.*c1,sm(db));dstring((uv-vec2 (-.3,.15)),41.,.02,da);db=smoothstep(48.,49.,iTime)-smoothstep(53.,54.,iTime);da=mix(1.,da,db);stroke(da-.005,.0025,db);new.gba=mix(new.gba,c.yyy,sm(db));new.gba=mix(new.gba,c1,sm(da-.0025));db=da+.0025;new.gba=mix(new.gba,2.*c1,sm(db));new.gba=clamp(new.gba,0.,1.);}else if(iTime<130.){float da,db;dbox(uv-vec2 (.05,.3),vec2 (1.6,.055),da);da=mix(1.,da,smoothstep(125.,126.,iTime));new.gba=mix(new.gba,mix(new.gba,c.xxx,.5),sm(da));dstring((uv-vec2 (-.4,.3)),28.,.05,da);lfnoise(55.*uv,db);stroke(db,0.535,db);vec3 c1=mix(mix(new.gba,c.yyy,.3),c.yyy,sm(db/50.));db=smoothstep(125.,126.,iTime);da=mix(1.,da,db);new.gba=mix(new.gba,c1,sm(da));stroke(da-.01,.001,da);new.gba=mix(new.gba,c.yyy,sm(da));}else {float da;dstring(uv-vec2 (-.55,0.),27.,.025,da);da=mix(1.,da,smoothstep(172.,172.5,iTime));new.gba=mix(new.gba,c.yyy,smoothstep(172.,172.5,iTime));new.gba=mix(new.gba,vec3 (.9,.2,.03),sm(da));stroke(da-.005,.001,da);new.gba=mix(new.gba,c.xxx,sm(da));}float dc;dbox(uv,.5*vec2 (a,1.),dc);stroke(dc,.005,dc);new.gba=mix(new.gba,c.yyy,sm(dc));fragColor=vec4 (new.gba,1.);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *post_source = "#version 130\nuniform float iFSAA;uniform vec2 iResolution;uniform sampler2D iChannel0;uniform float iTime;out vec4 gl_FragColor;const float pi=acos(-1.);const vec3 c=vec3 (1.,0.,-1.);float a=1.0;float lscale,rscale;float size;float nbeats;float iScale;void rand(in vec2 x,out float n);void lfnoise(in vec2 t,out float n);void stroke(in float d0,in float s,out float d);void dvoronoi(in vec2 x,out float d,out vec2 z);float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}float dot2(in vec3 v){return dot(v,v);}void dtriangle3(in vec3 p,in vec3 v1,in vec3 v2,in vec3 v3,out float dst){vec3 v21=v2-v1;vec3 p1=p-v1;vec3 v32=v3-v2;vec3 p2=p-v2;vec3 v13=v1-v3;vec3 p3=p-v3;vec3 nor=cross(v21,v13);dst=sqrt((sign(dot(cross(v21,nor),p1))+sign(dot(cross(v32,nor),p2))+sign(dot(cross(v13,nor),p3))<2.0)?min(min(dot2(v21*clamp(dot(v21,p1)/dot2(v21),0.0,1.0)-p1),dot2(v32*clamp(dot(v32,p2)/dot2(v32),0.0,1.0)-p2)),dot2(v13*clamp(dot(v13,p3)/dot2(v13),0.0,1.0)-p3)):dot(nor,p1)*dot(nor,p1)/dot2(nor));}void rot3(in vec3 p,out mat3 rot);void dbox3(in vec3 x,in vec3 b,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);mat3 R;void scene(in vec3 x,out vec2 sdf){float d;dbox3(x,.2*c.xxx,sdf.x);sdf.y=1.;dbox3(x-.1*c.xyy,vec3 (.02,.3,.12),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x-.05*c.xyy-.1*c.yyx,vec3 (.07,.3,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x,vec3 (.02,.3,.1),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x+.05*c.xyy+.1*c.yyx,vec3 (.07,.3,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x+.1*c.xyy-.1*c.yyx,vec3 (.02,.3,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x+.04*c.yyx,vec3 (.3,.02,.08),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x-.1*c.yyx,vec3 (.3,.02,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));vec3 y=vec3 (x.x,abs(x.y),x.z);dbox3(y-.05*c.yxy,vec3 (.1,.03,.3),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(y-.1*c.yxy-.06*c.xyy,vec3 (.08,.021,.3),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));vec3 z=vec3 (abs(x.x),x.yz);dbox3(z-.119*c.xyy,vec3 (.021,.08,.3),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));}void scene2(in vec3 x,out vec2 sdf){float v=0.;vec2 vi=c.yy;dvoronoi(x.xy/size,v,vi);vec3 y=vec3 (x.xy-vi*size,x.z);vec2 yi=vi*size;float n=0.;lfnoise(4.*(yi-.5*iTime),n);lfnoise(12.*vec2 (n,1.)*yi-(.8+.2*n)*c.xy,n);n*=iScale;sdf=vec2 (length(y-.05*n*c.yyx)-mix(.05,1.,length(texture(iChannel0,yi/vec2 (a,1.)).rgb)/sqrt(3.))*size,1.);}void normal2(in vec3 x,out vec3 n,in float dx){vec2 s,na;scene2(x,s);scene2(x+dx*c.xyy,na);n.x=na.x;scene2(x+dx*c.yxy,na);n.y=na.x;scene2(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);}void scene3(in vec3 x,out vec2 sdf){vec3 y=vec3 (mod(x.xy,2.*size)-size,x.z);vec2 yi=x.xy-y.xy;float ss=mix(.0,.05,size/.01);vec2 p0=.8*size*c.xx,p1=.8*size*c.zx,p2=.8*size*c.xz;vec2 ind;float y0,y1,y2;lfnoise(4.e1*(yi+p0-.5e-4*iTime),y0);lfnoise(12.e1*vec2 (y0,1.)*(yi+p0)-1.e-4*(.8+.2*y0)*iTime*c.xy,y0);lfnoise(4.e1*(yi+p1-.5e-4*iTime),y1);lfnoise(12.e1*vec2 (y1,1.)*(yi+p1)-1.e-4*(.8+.2*y1)*iTime*c.xy,y1);lfnoise(4.e1*(yi+p2-.5e-4*iTime),y2);lfnoise(12.e1*vec2 (y2,1.)*(yi+p2)-1.e-4*(.8+.2*y2)*iTime*c.xy,y2);y0*=ss;y1*=ss;y2*=ss;dtriangle3(y,vec3 (p0,y0),vec3 (p1,y1),vec3 (p2,y2),sdf.x);float d;vec2 p3=.8*size*c.zz,p4=.8*size*c.xz,p5=.8*size*c.zx;float y3,y4,y5;lfnoise(4.e1*(yi+p3-.5e-4*iTime),y3);lfnoise(12.e1*vec2 (y3,1.)*(yi+p3)-1.e-4*(.8+.2*y3)*iTime*c.xy,y3);lfnoise(4.e1*(yi+p4-.5e-4*iTime),y4);lfnoise(12.e1*vec2 (y4,1.)*(yi+p4)-1.e-4*(.8+.2*y4)*iTime*c.xy,y4);lfnoise(4.e1*(yi+p5-.5e-4*iTime),y5);lfnoise(12.e1*vec2 (y5,1.)*(yi+p5)-1.e-4*(.8+.2*y5)*iTime*c.xy,y5);y3*=ss;y4*=ss;y5*=ss;dtriangle3(y,vec3 (p3,y3),vec3 (p4,y4),vec3 (p5,y5),d);sdf.x=min(sdf.x,d);stroke(sdf.x,.1*size,sdf.x);sdf.y=1.;}void normal3(in vec3 x,out vec3 n,in float dx){vec2 s,na;scene3(x,s);scene3(x+dx*c.xyy,na);n.x=na.x;scene3(x+dx*c.yxy,na);n.y=na.x;scene3(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);}void normal(in vec3 x,out vec3 n,in float dx);void mainImage(out vec4 fragColor,in vec2 fragCoord_){vec2 fragCoord=fragCoord_;if(iTime<159.5456&&iTime>155.9092){vec2 n;lfnoise(22.*fragCoord/iResolution-3.*iTime,n.x);lfnoise(22.*fragCoord/iResolution-3.*iTime-1337.,n.y);fragCoord+=22.*n;}else if(iTime<165.&&iTime>163.182){vec2 n;lfnoise(22.*fragCoord/iResolution-3.*iTime,n.x);lfnoise(22.*fragCoord/iResolution-3.*iTime-1337.,n.y);fragCoord+=22.*n;}float a=iResolution.x/iResolution.y;vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0);nbeats=mod(iTime,60./29.);iScale=nbeats-30./29.;iScale=smoothstep(-5./29.,0.,iScale)*(1.-smoothstep(0.,15./29.,iScale));vec3 col=texture(iChannel0,fragCoord/iResolution).rgb;float delta=0.;rot3(vec3 (-2.*pi/8.,2.*pi/8.,2.*pi/4.)-iTime*vec3 (1.1,1.3,1.5),R);float d;vec2 s;vec3 o,r,u,t,ssize,dir,x,n;vec2 uv2=10.*(uv-vec2 (-.45*a,.45));o=R*c.yyx;r=c.xyy;u=c.yxy;t=c.yyy;int N=250,i;t=uv2.x*r+uv2.y*u;t=R*t;dir=normalize(t-o);ssize=.2*c.xxx;vec3 tlo=min((ssize-o)/dir,(-ssize-o)/dir);vec2 abxlo=abs(o.yz+tlo.x*dir.yz),abylo=abs(o.xz+tlo.y*dir.xz),abzlo=abs(o.xy+tlo.z*dir.xy);vec4 dn=100.*c.xyyy;dn=mix(dn,vec4 (tlo.x,c.xyy),float (all(lessThan(abxlo,ssize.yz)))*step(tlo.x,dn.x));dn=mix(dn,vec4 (tlo.y,c.yxy),float (all(lessThan(abylo,ssize.xz)))*step(tlo.y,dn.x));dn=mix(dn,vec4 (tlo.z,c.yyx),float (all(lessThan(abzlo,ssize.xy)))*step(tlo.z,dn.x));uv=(fragCoord)/iResolution.xy*vec2 (a,1.);d=dn.r;float nan;lfnoise(iTime*c.xx,nan);nan+=.5;if(nan>0.)d=3.;if(d<=2.){x=o+d*dir;scene(x,s);if(s.x>1.e-4){for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=s.x;}}if(i<N){normal(x,n,5.e-4);if(s.y==1.){vec3 l=normalize(x+c.zzx*vec3 (1.3,.9,1.2));col=vec3 (0.81,0.15,0.18);col=.3*col+.4*col*abs(dot(l,n))+.6*col*abs(pow(dot(reflect(-l,n),dir),2.));}else if(s.y==2.){vec3 l=normalize(x+c.zzx*vec3 (1.3,.9,1.2));col=.7*c.xxx;col=.5*col+.4*col*abs(dot(l,n))+.8*col*abs(pow(dot(reflect(-l,n),dir),2.));}}if(iTime<0.)col=texture(iChannel0,fragCoord/iResolution).rgb;}else {iScale=nbeats-30./29.;iScale=smoothstep(-5./29.,0.,iScale)*(1.-smoothstep(0./29.,35./29.,iScale));lscale=smoothstep(0.,.5,clamp((iTime-10.),0.,1.))*(1.-smoothstep(0.,.5,clamp((iTime-18.),0.,1.)));rscale=smoothstep(167.,167.5,iTime)-smoothstep(172.,172.5,iTime);size=mix(.005,.01,rscale);size=mix(0.,size,max(rscale,lscale));if(lscale>0.){col=c.yyy;o=c.yyx+.5*vec3 (cos(iTime),sin(iTime),0.);r=c.xyy;u=c.yxy;t=c.yyy;dir=c.yyy;n=c.yyy;x=c.yyy;N=200;t=uv.x*r+uv.y*u;dir=normalize(t-o);d=-(o.z-.05-.5*size)/dir.z;for(i=0;i<N;++i){x=o+d*dir;scene2(x,s);if(s.x<1.e-4)break;if(x.z<-.05-.5*size){col=c.yyy;i=N;break;}d+=min(s.x,1.e-3);}if(i<N){normal2(x,n,5.e-4);vec3 l=normalize(x+.5*n);if(s.y==1.){float v;vec2 vi;dvoronoi(x.xy/size,v,vi);vec3 y=vec3 (x.xy-vi*size,x.z);vec2 yi=vi*size;float bound=sqrt(iFSAA)-1.;for(float i=-.5*bound;i<=.5*bound;i+=1.)for(float j=-.5*bound;j<=.5*bound;j+=1.){col+=texture(iChannel0,yi/vec2 (a,1.)+vec2 (i,j)*3./max(bound,1.)/iResolution.xy).xyz;}col/=iFSAA;col=.4*col+.9*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);}}else col=c.yyy;}else if(rscale>0.){col=c.yyy;o=c.yyx+.5*vec3 (-1.,-1.,0.);r=c.xyy;u=c.yxy;t=c.yyy;dir=c.yyy;n=c.yyy;x=c.yyy;N=300;t=uv.x*r+uv.y*u;dir=normalize(t-o);d=-(o.z-.05-.5*size)/dir.z;for(i=0;i<N;++i){x=o+d*dir;scene3(x,s);if(s.x<1.e-4)break;if(x.z<-.05-.5*size){col=c.yyy;i=N;break;}d+=min(s.x,1.e-3);}if(i<N){normal3(x,n,5.e-4);vec3 l=normalize(x+.5*n);if(s.y==1.){vec3 y=vec3 (mod(x.xy,size)-.5*size,x.z);vec2 yi=x.xy-y.xy;col=texture(iChannel0,yi/vec2 (a,1.)).rgb;col=.4*col+.9*col*abs(dot(l,n))+.6*col*pow(abs(dot(reflect(-l,n),dir)),3.);}}else col=c.yyy;}}col+=vec3 (0.,0.05,0.1)*sin(uv.y*1050.+5.*iTime);fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *logo210_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.,0.,-1.);float a=1.0;float nbeats,iScale;void dbox3(in vec3 x,in vec3 b,out float d);void rot3(in vec3 p,out mat3 rot);void stroke(in float d0,in float s,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);void dbox210(in vec3 x,in float size,out vec2 sdf){x/=size;float d=1.;dbox3(x,.2*c.xxx,sdf.x);sdf.y=1.;dbox3(x-.1*c.xyy,vec3 (.02,.3,.12),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x-.05*c.xyy-.1*c.yyx,vec3 (.07,.3,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x,vec3 (.02,.3,.1),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x+.05*c.xyy+.1*c.yyx,vec3 (.07,.3,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x+.1*c.xyy-.1*c.yyx,vec3 (.02,.3,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x+.04*c.yyx,vec3 (.3,.02,.08),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(x-.1*c.yyx,vec3 (.3,.02,.02),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));vec3 y=vec3 (x.x,abs(x.y),x.z);dbox3(y-.05*c.yxy,vec3 (.1,.03,.3),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));dbox3(y-.1*c.yxy-.06*c.xyy,vec3 (.08,.021,.3),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));vec3 z=vec3 (abs(x.x),x.yz);dbox3(z-.119*c.xyy,vec3 (.021,.08,.3),d);sdf.x=max(-d,sdf.x);sdf.y=mix(sdf.y,2.,step(d,sdf.x));sdf.x*=size;}mat3 R;void scene(in vec3 x,out vec2 sdf){sdf=c.xy;float d,da;rot3(vec3 (-pi/2.,0.,pi/2.),R);x=R*x;vec2 sda=c.xy;dbox210(x+.1*c.xyy,.5,sdf);rot3(vec3 (pi/2.,0.,pi/2.),R);x=R*x;dbox210(x,5.,sda);add(sdf,sda,sdf);rot3(vec3 (pi/2.,-pi/2.,pi/2.),R);x=R*x;dbox210(x-2.*c.yxy,50.,sda);add(sdf,sda,sdf);stroke(sdf.x,.001,sdf.x);dbox3(x,100.*c.xxx,sda.x);sda.y=2.;add(sdf,sda*c.zx,sdf);}void normal(in vec3 x,out vec3 n,in float dx){vec2 s,na;scene(x,s);scene(x+dx*c.xyy,na);n.x=na.x;scene(x+dx*c.yxy,na);n.y=na.x;scene(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);}void mainImage(out vec4 fragColor,in vec2 fragCoord){float a=iResolution.x/iResolution.y;vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0);vec3 col=c.xxx;float d=0.;vec2 s;vec3 o,t,dir,x,n;mat3 Ra;rot3(mix(c.yyy,vec3 (-5.7884463,2.4242211,0.3463173),clamp((iTime-6.)/1.5,0.,1.)),Ra);o=Ra*mix(mix(mix(c.yyy-.1*c.yxy,c.yyx,clamp(iTime/2.,0.,1.)),10.*c.yyx,clamp((iTime-2.)/2.,0.,1.)),100.*c.yyx,clamp((iTime-4.)/2.,0.,1.));t=c.yyy;int N=650,i;dir=Ra*normalize(vec3 (uv,-1.));for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=s.x;}if(s.x<1.e-4){normal(x,n,5.e-4);vec3 l=normalize(x+.1*n);if(s.y==1.){col=vec3 (0.81,0.15,0.18);col=.3*col+.4*col*abs(dot(l,n))+.6*col*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==2.){col=.7*c.xxx;col=.5*col+.4*col*abs(dot(l,n))+.8*col*abs(pow(dot(reflect(-l,n),dir),3.));vec3 c1=c.yyy;o=x;dir=reflect(dir,n);d=1.e-1;N=150;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=s.x;}if(s.x<1.e-4){normal(x,n,5.e-4);vec3 l=normalize(x+.1*n);if(s.y==1.){c1=vec3 (0.81,0.15,0.18);c1=.3*c1+.4*c1*abs(dot(l,n))+.6*c1*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==2.){c1=.7*c.xxx;c1=.5*c1+.4*c1*abs(dot(l,n))+.8*c1*abs(pow(dot(reflect(-l,n),dir),3.));}c1=clamp(c1,0.,1.);col=mix(col,c1,.2);}col=clamp(col,0.,1.);}}col=mix(col,vec3 (0.20,0.01,0.14),smoothstep(0.,1.,iTime-10.));fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *transbubbles_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.,0.,-1.);float a=1.0;float nbeats,iScale;void scale(out float s);void rand(in vec2 x,out float n);void lfnoise(in vec2 t,out float n);void dbox3(in vec3 x,in vec3 b,out float d);void rot3(in vec3 p,out mat3 rot);void stroke(in float d0,in float s,out float d);void hash13(in vec3 p3,out float d);void add(in vec2 sda,in vec2 sdb,out vec2 sdf);void smoothmin(in float a,in float b,in float k,out float dst);void dvoronoi3(in vec3 x,out float d,out vec3 z){vec3 y=floor(x);float ret=1.;vec3 pf=c.yyy,p;float df=10.;for(int i=-1;i<=1;i+=1)for(int j=-1;j<=1;j+=1){for(int k=-1;k<=1;k+=1){p=y+vec3 (float (i),float (j),float (k));float pa;hash13(p,pa);p+=pa;d=length(x-p);if(d<df){df=d;pf=p;}}}for(int i=-1;i<=1;i+=1)for(int j=-1;j<=1;j+=1){for(int k=-1;k<=1;k+=1){p=y+vec3 (float (i),float (j),float (k));float pa;hash13(p,pa);p+=pa;vec3 o=p-pf;d=length(.5*o-dot(x-pf,o)/dot(o,o)*o);ret=min(ret,d);}}d=ret;z=pf;}mat3 R;vec3 ind;void scene(in vec3 x,out vec2 sdf){x=R*x;x-=mix(0.,.1*iTime,step(150.,iTime));sdf=c.xy;float bsize=2.;float d,da;dvoronoi3(bsize*x,d,ind);vec3 y=x-ind/bsize;float n;hash13(ind,n);add(sdf,vec2 (abs(length(y)-.3)-.001,2.),sdf);sdf.x=max(sdf.x,-length(x)+.5);}void normal(in vec3 x,out vec3 n,in float dx);void mainImage(out vec4 fragColor,in vec2 fragCoord){float a=iResolution.x/iResolution.y;vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0);rot3(.1*vec3 (1.1,1.3,1.5)*iTime,R);vec3 col=c.yyy;scale(iScale);float d=0.;vec2 s;vec3 o,t,dir,x,n;float mx=clamp((iTime-5.),0.,1.),my=clamp(iTime-10.,0.,1.);float nai;lfnoise(2.*iTime*c.xx,nai);nai=.5*(nai);o=c.yyx;t=c.yyy;int N=250,i;dir=normalize(vec3 (uv,-1.));for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=s.x;}if(s.x<1.e-4){normal(x,n,1.e-4);vec3 l=normalize(x+.1*n);if(s.y==2.){col=mix(vec3 (0.86,0.21,0.13),vec3 (0.02,0.46,0.44),mx);col=.1*col+.1*col*abs(dot(l,n))+.5*col*abs(pow(dot(reflect(-l,n),dir),2.));vec3 c1=c.yyy;for(float fraction=0.;fraction<=3.;fraction+=1.){o=x;vec3 ddir=refract(dir,n,.95+.05*nai);vec3 dddir=refract(dir,n,.1);dir=refract(dir,n,.5+.05*nai);dir=mix(ddir,dir,mx);dir=mix(dir,dddir,my);dir=mix(dir,ddir,step(150.,iTime));d=2.e-2;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=s.x;}if(s.x<1.e-4){normal(x,n,1.e-4);vec3 l=normalize(x+.1*n);if(s.y==2.){c1=(fraction==0.)?vec3 (0.86,0.21,0.13):(fraction==1.)?vec3 (0.85,0.80,0.62):(fraction==2.)?vec3 (0.22,0.25,0.25):(fraction==3.)?vec3 (0.16,0.17,0.17):vec3 (0.12,0.12,0.13);vec3 c2=(fraction==0.)?vec3 (0.99,0.33,0.05):(fraction==1.)?vec3 (0.94,0.94,0.94):(fraction==2.)?vec3 (0.75,0.82,0.88):(fraction==3.)?vec3 (0.25,0.34,0.39):vec3 (0.17,0.22,0.27);vec3 c3=(fraction==0.)?vec3 (1.00,0.55,0.03):(fraction==1.)?vec3 (0.84,0.20,0.18):(fraction==2.)?vec3 (0.13,0.55,0.57):(fraction==3.)?vec3 (0.29,0.22,0.30):vec3 (0.00,0.00,0.00);vec3 c4=(fraction==0.)?vec3 (0.09,0.00,0.13):(fraction==1.)?vec3 (0.87,0.38,0.61):(fraction==2.)?vec3 (0.97,0.45,0.50):(fraction==3.)?vec3 (0.97,0.70,0.72):vec3 (0.04,0.04,0.05);vec3 c5=(fraction==0.)?vec3 (0.94,0.24,0.27):(fraction==1.)?vec3 (0.48,0.49,0.55):(fraction==2.)?vec3 (0.81,0.76,0.90):(fraction==3.)?vec3 (0.44,0.33,0.60):vec3 (0.34,0.24,0.50);vec3 c6=(fraction==0.)?vec3 (0.11,0.32,0.23):(fraction==1.)?vec3 (0.25,0.89,0.79):(fraction==2.)?vec3 (0.31,0.59,0.56):(fraction==3.)?vec3 (0.00,0.34,0.39):vec3 (0.00,0.15,0.17);c1=mix(c1,c2,mx);c1=mix(c1,c3,my);c1=mix(c1,c4,step(150.,iTime));c1=mix(c1,c5,step(163.5,iTime));c1=mix(c1,c6,step(170.,iTime));c1=.1*c1+.4*c1*abs(dot(l,n))+mix(mix(5.8,3.79,mx),4.753,my)*c1*abs(pow(dot(reflect(-l,n),dir),2.));c1=mix(c1,2.*c1,smoothstep(mix(1.,.6,iScale),1.02,1.-abs(dot(n,c.xyy))));c1=mix(c1,2.*c1,smoothstep(mix(1.,.6,iScale),1.02,abs(dot(n,c.zyy))));}col=mix(col,c1,.15);}col=clamp(col,0.,1.);}}}col*=col;fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *volclouds_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.,0.,-1.);float a=1.0;float nbeats,iScale;void scale(out float s);void hash13(in vec3 p3,out float d){p3=fract(p3*.1031);p3+=dot(p3,p3.yzx+33.33);d=fract((p3.x+p3.y)*p3.z);}void lfnoise3(in vec3 t,out float num){t-=vec3 (11.,13.,5.);vec3 i=floor(t);t=fract(t);t=smoothstep(c.yyy,c.xxx,t);vec2 v1,v2,v3,v4;hash13(i,v1.x);hash13(i+c.xyy,v1.y);hash13(i+c.yxy,v2.x);hash13(i+c.xxy,v2.y);hash13(i+c.yyx,v3.x);hash13(i+c.xyx,v3.y);hash13(i+c.yxx,v4.x);hash13(i+c.xxx,v4.y);v1=c.zz+2.*mix(v1,v2,t.y);v3=c.zz+2.*mix(v3,v4,t.y);v2.x=mix(v1.x,v1.y,t.x);v2.y=mix(v3.x,v3.y,t.x);num=mix(v2.x,v2.y,t.z);}void mfnoise3(in vec3 x,in float d,in float b,in float e,out float n){n=0.;float a=1.,nf=0.,buf;for(float f=d;f<b;f*=2.){lfnoise3(f*x-vec3 (11.,13.,5.),buf);n+=a*buf;a*=e;nf+=1.;}n*=(1.-e)/(1.-pow(e,nf));}void stroke(in float d0,in float s,out float d){d=abs(d0)-s;}void rot3(in vec3 p,out mat3 rot){rot=mat3 (c.xyyy,cos(p.x),sin(p.x),0.,-sin(p.x),cos(p.x))*mat3 (cos(p.y),0.,-sin(p.y),c.yxy,sin(p.y),0.,cos(p.y))*mat3 (cos(p.z),-sin(p.z),0.,sin(p.z),cos(p.z),c.yyyx);}mat3 R;vec3 ind;void scene(in vec3 x,out vec2 sdf){x=R*x;float n;mfnoise3(x-.1*iTime*c.yyx,10.,400.,.25,n);n=.5+.5*n;sdf=vec2 (-n,2.);}void normal(in vec3 x,out vec3 n,in float dx){vec2 s,na;scene(x,s);scene(x+dx*c.xyy,na);n.x=na.x;scene(x+dx*c.yxy,na);n.y=na.x;scene(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);}void palette2(in float scale,out vec3 col){const int N=5;const vec3 colors[N]=vec3 [N](vec3 (0.82,0.27,0.13),vec3 (0.85,0.77,0.68),vec3 (0.65,0.59,0.55),vec3 (0.45,0.29,0.24),vec3 (0.85,0.27,0.15));float index=floor(scale*float (N)),remainder=scale*float (N)-index;col=mix(colors[int (index)],colors[int (index)+1],remainder);}void palette3(in float scale,out vec3 col){const int N=5;const vec3 colors[N]=vec3 [N](vec3 (0.86,0.21,0.13),vec3 (0.85,0.80,0.62),vec3 (0.22,0.25,0.25),vec3 (0.16,0.17,0.17),vec3 (0.12,0.12,0.13));float index=floor(scale*float (N)),remainder=scale*float (N)-index;col=mix(colors[int (index)],colors[int (index)+1],remainder);}void palette1(in float scale,out vec3 col){const int N=5;vec3 colors[N];if(iTime<150.)colors=vec3 [N](vec3 (0.99,0.33,0.05),vec3 (0.94,0.94,0.94),vec3 (0.75,0.82,0.88),vec3 (0.25,0.34,0.39),vec3 (0.17,0.22,0.27));else if(iTime<160.)colors=vec3 [N](vec3 (0.82,0.27,0.13),vec3 (0.85,0.77,0.68),vec3 (0.65,0.59,0.55),vec3 (0.45,0.29,0.24),vec3 (0.85,0.27,0.15));else if(iTime<165.)colors=vec3 [N](vec3 (0.11,0.32,0.23),vec3 (0.25,0.89,0.79),vec3 (0.31,0.59,0.56),vec3 (0.00,0.34,0.39),vec3 (0.00,0.15,0.17));else colors=vec3 [N](vec3 (0.86,0.21,0.13),vec3 (0.85,0.80,0.62),vec3 (0.22,0.25,0.25),vec3 (0.16,0.17,0.17),vec3 (0.12,0.12,0.13));float index=floor(scale*float (N)),remainder=scale*float (N)-index;col=mix(colors[int (index)],colors[int (index)+1],remainder);}void mainImage(out vec4 fragColor,in vec2 fragCoord){rot3(mix(mix(mix(0.,pi/2.,clamp((iTime/5.),0.,1.)),0.,clamp((iTime-5.)/5.,0.,1.)),pi/2.,clamp((iTime-10.)/5.,0.,1.))*c.xyx,R);scale(iScale);float a=iResolution.x/iResolution.y;vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0);vec3 col=c.yyy;float d=0.;vec2 s;vec3 o,t,dir,x,n;o=c.yyx;t=c.yyy;int N=40,i;dir=normalize(vec3 (uv,-1.));for(i=0;i<N;++i){d+=.5/float (N);x=o+d*dir;scene(x,s);normal(x,n,5.e-4);vec3 l=normalize(x+.1*n);vec3 c1;palette1(-s.x,c1);if(iTime<150.){vec3 c2;palette2(-s.x,c2);c1=mix(c1,c2,smoothstep(5.,6.,iTime));palette3(-s.x,c2);c1=mix(c1,c2,smoothstep(10.,11.,iTime));}c1=.1*c1+.1*c1*abs(dot(l,n))+1.354*c1*abs(pow(dot(reflect(-l,n),dir),2.));c1=mix(c1,2.*c1,smoothstep(mix(1.,.6,iScale),1.02,1.-abs(dot(n,c.xyy))));c1=mix(c1,2.*c1,smoothstep(mix(1.,.6,iScale),1.02,abs(dot(n,c.zyy))));c1=clamp(c1,0.,1.);col=mix(col,c1,d*d);}col*=col;fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
const char *chart_source = "#version 130\nuniform float iTime;uniform vec2 iResolution;const float pi=acos(-1.);const vec3 c=vec3 (1.,0.,-1.);float a=1.0;float nbeats,iScale;void rand(in vec2 x,out float n){x+=400.;n=fract(sin(dot(sign(x)*abs(x),vec2 (12.9898,78.233)))*43758.5453);}void lfnoise(in vec2 t,out float n){vec2 i=floor(t);t=fract(t);t=smoothstep(c.yy,c.xx,t);vec2 v1,v2;rand(i,v1.x);rand(i+c.xy,v1.y);rand(i+c.yx,v2.x);rand(i+c.xx,v2.y);v1=c.zz+2.*mix(v1,v2,t.y);n=mix(v1.x,v1.y,t.x);}void mfnoise(in vec2 x,in float d,in float b,in float e,out float n){n=0.;float a=1.,nf=0.,buf;for(float f=d;f<b;f*=2.){lfnoise(f*x,buf);n+=a*buf;a*=e;nf+=1.;}n*=(1.-e)/(1.-pow(e,nf));}void dbox3(in vec3 x,in vec3 b,out float d){vec3 da=abs(x)-b;d=length(max(da,0.0))+min(max(da.x,max(da.y,da.z)),0.0);}void rot3(in vec3 p,out mat3 rot){rot=mat3 (c.xyyy,cos(p.x),sin(p.x),0.,-sin(p.x),cos(p.x))*mat3 (cos(p.y),0.,-sin(p.y),c.yxy,sin(p.y),0.,cos(p.y))*mat3 (cos(p.z),-sin(p.z),0.,sin(p.z),cos(p.z),c.yyyx);}void stroke(in float d0,in float s,out float d){d=abs(d0)-s;}void add(in vec2 sda,in vec2 sdb,out vec2 sdf){sdf=sda.x<sdb.x?sda:sdb;}void rotAB(in vec3 a,in vec3 b,out mat3 R){a=normalize(a);b=normalize(b);vec3 v=cross(a,b);float co=dot(a,b);R=mat3 (0.,v.z,-v.y,-v.z,0.,v.x,v.y,-v.x,0.);R=R*R/(1.+co)+R;R+=mat3 (1.);}mat3 R;float ind;void scene(in vec3 x,out vec2 sdf){sdf=3.*c.yx;float d,da;dbox3(x,2.*c.xxx,sdf.x);sdf.x*=-1.;sdf.x=min(sdf.x,x.z);vec3 y=x;x=vec3 (x.x,mod(x.y,.2)-.1,x.z);vec3 xi=(y-x)/.2;ind=xi.y;float dx;lfnoise(ind*c.xx,dx);float n;mfnoise((x.x+.3*iTime)*c.xx-xi.y,12.,120.,.45,n);n*=clamp(.2-.01*xi.y,0.,1.);vec2 sda;dbox3(x-1.4*dx*c.xyy,vec3 (.5,.01,.3*(.5+.5*dx)-n),sda.x);sda.y=1.;add(sdf,sda,sdf);dbox3(x-(.3*(.5+.5*dx)-n)*c.yyx-1.4*dx*c.xyy,vec3 (.5,.03,.01),sda.x);sda.y=4.;add(sdf,sda,sdf);}void normal(in vec3 x,out vec3 n,in float dx){vec2 s,na;scene(x,s);scene(x+dx*c.xyy,na);n.x=na.x;scene(x+dx*c.yxy,na);n.y=na.x;scene(x+dx*c.yyx,na);n.z=na.x;n=normalize(n-s.x);}float sm(float d){return smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d);}void palette1(in float scale,out vec3 col){const int N=8;const vec3 colors[N]=vec3 [N](vec3 (0.18,0.30,0.65),vec3 (0.00,0.69,0.80),vec3 (0.45,0.72,0.50),vec3 (0.93,0.36,0.44),vec3 (1.00,0.51,0.45),vec3 (0.78,0.57,0.78),vec3 (0.54,0.33,0.60),vec3 (0.98,0.78,0.38));float index=floor(scale*float (N)),remainder=scale*float (N)-index;col=mix(colors[int (index)],colors[int (index)+1],remainder);}void mainImage(out vec4 fragColor,in vec2 fragCoord){float a=iResolution.x/iResolution.y;vec2 uv=fragCoord/iResolution.yy-0.5*vec2 (a,1.0);vec3 col=c.yyy,o=c.yzx,r=c.xyy,u=normalize(c.yxx),t=c.yyy,dir,n,x;int N=400,i;t=uv.x*r+uv.y*u;dir=normalize(t-o);float d=0.;vec2 s;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=min(s.x,7.e-3);}if(s.x<1.e-4){normal(x,n,5.e-4);vec3 l=normalize(x+.1*n);if(s.y==1.){ind=.5+.5*clamp(ind/10.,-1.,1.);palette1(ind,col);vec3 c1;palette1(ind+.1,c1);col=mix(col,c1,clamp(1.-x.z/.3,0.,1.));col=.3*col+.4*col*abs(dot(l,n))+.9*col*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==4.){col=c.xxx;col=.5*col+.4*col*abs(dot(l,n))+.8*col*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==2.||s.y==3.){col=c.xxx;vec2 ma=abs(mod(x.xy,.05)-.025)-.0015;col=mix(col,.5*c.xxx,sm(min(ma.x,ma.y)));col=.5*col+.4*col*abs(dot(l,n))+.8*col*abs(pow(dot(reflect(-l,n),dir),3.));vec3 c1=c.yyy;o=x;dir=reflect(dir,n);d=1.e-2;N=100;for(i=0;i<N;++i){x=o+d*dir;scene(x,s);if(s.x<1.e-4)break;d+=s.x;}if(s.x<1.e-4){normal(x,n,5.e-4);vec3 l=normalize(x+.1*n);if(s.y==1.){ind=.5+.5*clamp(ind/10.,-1.,1.);palette1(ind,c1);vec3 c2;palette1(ind+.1,c2);c1=mix(c1,c1,clamp(1.-x.z/.3,0.,1.));c1=.3*c1+.4*c1*abs(dot(l,n))+.9*c1*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==2.){c1=.7*c.xxx;c1=.5*c1+.4*c1*abs(dot(l,n))+.8*c1*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==3.){c1=.7*c.xxx;c1=.5*c1+.4*c1*abs(dot(l,n))+.8*c1*abs(pow(dot(reflect(-l,n),dir),3.));}else if(s.y==4.){c1=c.xxx;c1=.5*c1+.4*c1*abs(dot(l,n))+.8*c1*abs(pow(dot(reflect(-l,n),dir),3.));}c1=clamp(c1,0.,1.);col=mix(col,c1,.2);}col=clamp(col,0.,1.);}}fragColor=vec4 (clamp(col,0.,1.),1.0);}void main(){mainImage(gl_FragColor,gl_FragCoord.xy);}\0";
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
int voronoidesign_iTime_location,voronoidesign_iResolution_location;
int groundboxes_iTime_location,groundboxes_iResolution_location;
int graffiti_iTime_location,graffiti_iResolution_location;
int greet_iTime_location,greet_iResolution_location;
int evoke_iTime_location,evoke_iResolution_location;
int canal_iTime_location,canal_iResolution_location;
int text_iFontWidth_location,text_iTime_location,text_iResolution_location,text_iChannel0_location,text_iFont_location,text_iFSAA_location;
int post_iFSAA_location,post_iResolution_location,post_iChannel0_location,post_iTime_location;
int logo210_iTime_location,logo210_iResolution_location;
int transbubbles_iTime_location,transbubbles_iResolution_location;
int volclouds_iTime_location,volclouds_iResolution_location;
int chart_iTime_location,chart_iResolution_location;
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
    glAttachShader(canal_program,rot3_handle);
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
