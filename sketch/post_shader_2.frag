// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

float size = .01;

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void dvoronoi(in vec2 x, out float d, out vec2 z)
{
    vec2 y = floor(x);
       float ret = 1.;
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    d = ret;
    z = pf;
}

float dot2( in vec3 v ) { return dot(v,v); }

// Adapted from https://www.shadertoy.com/view/4sXXRN
void dtriangle3(in vec3 p,  in vec3 v1, in vec3 v2, in vec3 v3, out float dst)
{
    vec3 v21 = v2 - v1; vec3 p1 = p - v1;
    vec3 v32 = v3 - v2; vec3 p2 = p - v2;
    vec3 v13 = v1 - v3; vec3 p3 = p - v3;
    vec3 nor = cross( v21, v13 );

    dst = sqrt( (sign(dot(cross(v21,nor),p1)) + 
                  sign(dot(cross(v32,nor),p2)) + 
                  sign(dot(cross(v13,nor),p3))<2.0) 
                  ?
                  min( min( 
                  dot2(v21*clamp(dot(v21,p1)/dot2(v21),0.0,1.0)-p1), 
                  dot2(v32*clamp(dot(v32,p2)/dot2(v32),0.0,1.0)-p2) ), 
                  dot2(v13*clamp(dot(v13,p3)/dot2(v13),0.0,1.0)-p3) )
                  :
                  dot(nor,p1)*dot(nor,p1)/dot2(nor) );
}

float sm(float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

void scene(in vec3 x, out vec2 sdf)
{
    vec3 y = vec3(mod(x.xy,size)-.5*size, x.z);
    vec2 yi = x.xy-y.xy;
	float ss = .05;
    
    vec2 p0 = .4*size*c.xx,
        p1 = .4*size*c.zx,
        p2 = .4*size*c.xz;
    
    vec2 ind;
    
    float y0, y1, y2;
    lfnoise(4.e1*(yi+p0-.5e-2*iTime), y0);
    lfnoise(12.e1*vec2(y0,1.)*(yi+p0)-1.e-2*(.8+.2*y0)*iTime*c.xy, y0);
    lfnoise(4.e1*(yi+p1-.5e-2*iTime), y1);
    lfnoise(12.e1*vec2(y1,1.)*(yi+p1)-1.e-2*(.8+.2*y1)*iTime*c.xy, y1);
    lfnoise(4.e1*(yi+p2-.5e-2*iTime), y2);
    lfnoise(12.e1*vec2(y2,1.)*(yi+p2)-1.e-2*(.8+.2*y2)*iTime*c.xy, y2);
    y0 *= ss;
    y1 *= ss;
    y2 *= ss;
    
    dtriangle3(y, vec3(p0,y0), vec3(p1,y1), vec3(p2,y2), sdf.x);
    
    float d;
    vec2 p3 = .4*size*c.zz,
        p4 = .4*size*c.xz,
        p5 = .4*size*c.zx;
    
    float y3, y4, y5;
    lfnoise(4.e1*(yi+p3-.5e-2*iTime), y3);
    lfnoise(12.e1*vec2(y3,1.)*(yi+p3)-1.e-2*(.8+.2*y3)*iTime*c.xy, y3);
    lfnoise(4.e1*(yi+p4-.5e-2*iTime), y4);
    lfnoise(12.e1*vec2(y4,1.)*(yi+p4)-1.e-2*(.8+.2*y4)*iTime*c.xy, y4);
    lfnoise(4.e1*(yi+p5-.5e-2*iTime), y5);
    lfnoise(12.e1*vec2(y5,1.)*(yi+p5)-1.e-2*(.8+.2*y5)*iTime*c.xy, y5);
    y3 *= ss;
    y4 *= ss;
    y5 *= ss;
    
    dtriangle3(y, vec3(p3,y3), vec3(p4,y4), vec3(p5,y5), d);
    sdf.x = min(sdf.x, d);
    /*
    p3 = .01*size*c.zx;
	p4 = .49*size*c.zy;
	p5 = .49*size*c.yx;
    
    lfnoise(4.*(yi+p3-.5*iTime), y3);
    lfnoise(12.*vec2(y3,1.)*(yi+p3)-1.e-2*(.8+.2*y3)*iTime*c.xy, y3);
    lfnoise(4.*(yi+p4-.5*iTime), y4);
    lfnoise(12.*vec2(y4,1.)*(yi+p4)-1.e-2*(.8+.2*y4)*iTime*c.xy, y4);
    lfnoise(4.*(yi+p5-.5*iTime), y5);
    lfnoise(12.*vec2(y5,1.)*(yi+p5)-1.e-2*(.8+.2*y5)*iTime*c.xy, y5);
    y3 *= ss;
    y4 *= ss;
    y5 *= ss;
    
    dtriangle3(y, vec3(p3,y3), vec3(p4,y4), vec3(p5,y5), d);
    sdf.x = min(sdf.x, d);
    
    p3 = .01*size*c.xz;
	p4 = .49*size*c.xy;
	p5 = .49*size*c.yz;
    
    lfnoise(4.*(yi+p3-.5*iTime), y3);
    lfnoise(12.*vec2(y3,1.)*(yi+p3)-1.e-2*(.8+.2*y3)*iTime*c.xy, y3);
    lfnoise(4.*(yi+p4-.5*iTime), y4);
    lfnoise(12.*vec2(y4,1.)*(yi+p4)-1.e-2*(.8+.2*y4)*iTime*c.xy, y4);
    lfnoise(4.*(yi+p5-.5*iTime), y5);
    lfnoise(12.*vec2(y5,1.)*(yi+p5)-1.e-2*(.8+.2*y5)*iTime*c.xy, y5);
    y3 *= ss;
    y4 *= ss;
    y5 *= ss;
    
    dtriangle3(y, vec3(p3,y3), vec3(p4,y4), vec3(p5,y5), d);
    sdf.x = min(sdf.x, d);
    */
    stroke(sdf.x, .1*size, sdf.x);
    sdf.y = 1.;
    
    //sdf = vec2(length(y)-.5*size, 1.);
    
    
    /*
    float v;
    vec2 vi;
    dvoronoi(x.xy/size, v, vi);
    vec3 y = vec3(x.xy-vi*size, x.z);
    vec2 yi = vi*size;
    
    float n;
    lfnoise(4.*(yi-.5*iTime), n);
    lfnoise(12.*vec2(n,1.)*yi-1.e-2*(.8+.2*n)*iTime*c.xy, n);
    sdf = vec2(length(y-.05*n*c.yyx)-.5*size, 1.);
	*/
}

void normal(in vec3 x, out vec3 n, in float dx)
{
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = (fragCoord)/iResolution.xy*vec2(a,1.),
        s;
    
    vec3 col = c.yyy, 
        o = c.yyx+.5*vec3(cos(iTime), sin(iTime),0.),
        r = c.xyy,
        u = c.yxy,
        t = c.yyy, 
        dir,
        n,
        x;
    int N = 300,
        i;
    t = uv.x * r + uv.y * u;
    dir = normalize(t-o);

    float d = -(o.z-.05)/dir.z;
    
    for(i = 0; i<N; ++i)
    {
     	x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4)break;
        
        if(x.z<-.05-.5*size)
        {
            col = c.yyy;
            i = N;
            break;
        }
        d += min(s.x,1.e-3);
        //d += s.x;
    }
    
    if(i < N)
    {
        normal(x,n, 5.e-4);
        vec3 l = normalize(x+.5*n);
       
		if(s.y == 1.)
        {
            vec3 y = vec3(mod(x.xy,size)-.5*size, x.z);
    		vec2 yi = x.xy-y.xy;
            
            col = texture(iChannel0, yi/vec2(a,1.)).rgb;
            
            col = .4*col
                + .9*col * abs(dot(l,n))
                + .6*col * pow(abs(dot(reflect(-l,n),dir)),3.);
        }
    }
    
    //col = texture(iChannel0, uv/vec2(a,1.)).rgb;
    
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    
    fragColor = vec4(clamp(col,0.,1.),1.0);
}
